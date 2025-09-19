############################################################
# modules.R
# Modular functions for StaBiCut pipeline
############################################################

########################################
# 1. Winsorize Expression
########################################
# winsorize(): Limit influence of outliers in expression data
winsorize <- function(x, lower = 0.05, upper = 0.95) {
  qnt <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
  x[x < qnt[1]] <- qnt[1]
  x[x > qnt[2]] <- qnt[2]
  return(x)
}


########################################
# 2. Cutoff & Cox Regression
########################################
# compute_cutoff_and_cox(): Determine cutoff and compute Cox HR and p-value
compute_cutoff_and_cox <- function(df,gene,cut.type = "optimal", minprop = 0.25) {
  # Determine the optimal cut-off point
  cutoff <- NA
  cutoff_p <- NA
  if (cut.type == "optimal") {
    res.cut <- surv_cutpoint(df, time = "time", event = "event", variables ="expr", minprop = 0.25,progressbar = F) 
    df_cat <- surv_categorize(res.cut)
    colnames(df_cat)[3] <- "group"
    
    df$group <- df_cat[[3]]
    
    # Extract cutoff value
    cutoff <- res.cut$cutpoint$cutpoint
    
    # Calculate the P value corresponding to the cut-off
    surv_test <- survdiff(Surv(time, event) ~ group, data = df)
    cutoff_p <- 1 - pchisq(surv_test$chisq, df = 1)
    
    cat(gene, " cutoff = ", round(cutoff, 3), ", p = ", signif(cutoff_p, 3), "\n")  # Print gene names
    
  } else if (cut.type == "median") {
    cutoff <- median(df$expr, na.rm = TRUE)
    df$group <- ifelse(df$expr > cutoff, "high", "low")
    
    surv_test <- survdiff(Surv(time, event) ~ group, data = df)
    cutoff_p <- 1 - pchisq(surv_test$chisq, df = 1)
    
  } else if (cut.type == "quantile") {
    q25 <- quantile(df$expr, 0.25)
    q75 <- quantile(df$expr, 0.75)
    df <- df[df$expr <= q25 | df$expr >= q75, ]
    df$group <- ifelse(df$expr >= q75, "high", "low")
    cutoff <- c(q25 = q25, q75 = q75)
    
    surv_test <- survdiff(Surv(time, event) ~ group, data = df)
    cutoff_p <- 1 - pchisq(surv_test$chisq, df = 1)
    
  } else {
    stop("cut.type must be one of: optimal, median, quantile")
  }
  
  # HR values ​​calculated by Cox regression
  cox_model <- tryCatch({
    coxph(Surv(time, event) ~ group, data = df)
  }, error = function(e) {
    cat("  Cox model failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(cox_model)) next
  cox_summary <- summary(cox_model)
  
  hr <- cox_summary$coefficients[,"exp(coef)"]
  ci_low <- cox_summary$conf.int[,"lower .95"]
  ci_high <- cox_summary$conf.int[,"upper .95"]
  p_val <- cox_summary$coefficients[,"Pr(>|z|)"]
   
  return(list(
    cutoff = cutoff,
    cutoff_p = cutoff_p,
    hr = cox_summary$coefficients[,"exp(coef)"],
    ci_low = cox_summary$conf.int[,"lower .95"],
    ci_high = cox_summary$conf.int[,"upper .95"],
    p_val = cox_summary$coefficients[,"Pr(>|z|)"],
    df = df
  ))
}


########################################
# 3. Spline Analysis
########################################
# evaluate_spline_support(): Evaluate spline effect and generate spline plot
evaluate_spline_support <- function(df, gene, plot_dir) {
  fit_null <- coxph(Surv(time, event) ~ 1, data = df)
  fit_spline <- coxph(Surv(time, event) ~ ns(expr, df = 3), data = df)
  
  p_spline <- NA
  
  if (!is.null(fit_spline)) {
    term <- tryCatch(termplot(fit_spline, se = TRUE, plot = FALSE), error = function(e) NULL)
    
    if (!is.null(term)) {
      spline_df <- data.frame(expr = term$expr$x,
                              fit = term$expr$y,
                              upper = term$expr$y + 2 * term$expr$se,
                              lower = term$expr$y - 2 * term$expr$se)
      
      # The p-value is obtained by comparing the null and spline fitting results.
      anova_result <- tryCatch(anova(fit_null, fit_spline), error = function(e) NULL)
      if (!is.null(anova_result) && nrow(anova_result) >= 2 && "Chisq" %in% colnames(anova_result)) {
        chisq_val <- anova_result$Chisq[2]
        df_val <- anova_result$Df[2]
        if (!is.na(chisq_val) && !is.na(df_val)) {
          p_spline <- pchisq(chisq_val, df = df_val, lower.tail = FALSE)
        }
      }
      spline_support <- if (is.numeric(p_spline) && !is.na(p_spline) && p_spline < 0.05) TRUE else FALSE
      
      # Plot Output
      spline_plot <- ggplot(spline_df, aes(x = expr, y = fit)) +
        geom_line(color = "#2C7BB6", size = 1.2) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#ABD9E9", alpha = 0.5) +
        geom_rug(data = df, aes(x = expr), inherit.aes = FALSE, sides = "b", alpha = 0.4, color = "gray40") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
        labs(title = paste("Spline effect of", gene),
             subtitle = paste0("P = ", if (!is.na(p_spline)) signif(p_spline, 3) else "NA"),
             x = paste(gene, "expression (log2 TPM)"),
             y = "Log hazard ratio (partial effect)"
        ) + theme_minimal(base_size = 14)
      
      ggsave(filename = file.path(plot_dir, paste0(gene, "_spline_effect_ggplot.pdf")),
             plot = spline_plot, width = 7, height = 5)
    }
  }
  
  return(list(p_spline = p_spline, spline_support = spline_support))
}


########################################
# 4. Bootstrap Stability Analysis
########################################
# compute_bootstrap_cutoffs(): Generate bootstrap distributions of cutoffs
compute_bootstrap_cutoffs <- function(df, n_boot = 1000, minprop = 0.25) {
  
  replicate(n_boot, {
    idx <- sample(1:nrow(df), replace = TRUE)
    df_bs <- df[idx, ]
    # Use the same method as the original analysis (surv_cutpoint)
    tryCatch({
      bs_cut <- surv_cutpoint(df_bs, 
                              time = "time", 
                              event = "event", 
                              variables = "expr", 
                              minprop = 0.25)
      
      val <- as.numeric(bs_cut$cutpoint$cutpoint) # Returns the cutpoint of the expr column
      if (is.numeric(val)) as.numeric(val) else NA
    }, error = function(e) { NA })
  })
}

#is_unimodal(): Check whether a distribution is unimodal
is_unimodal <- function(x, lambda_thresh = 0.3) {
  x <- x[!is.na(x)]
  if (length(x) < 30) return(FALSE)
  mixmdl <- tryCatch({
    suppressMessages(suppressWarnings(normalmixEM(x, k = 2)))
  }, error = function(e) NULL)
  return(!is.null(mixmdl) && mixmdl$lambda[2] < lambda_thresh)
}

# evaluate_bootstrap_stability(): Assess bootstrap stability and provide rationale
evaluate_bootstrap_stability <- function(bootstrap_cutoffs, iqr_thresh = 1.5, lambda_thresh = 0.3) {
  boot_iqr <- IQR(bootstrap_cutoffs, na.rm = TRUE)
  if (length(na.omit(bootstrap_cutoffs)) < 10) {
    stable <- FALSE
  } else {
    stable <- boot_iqr < iqr_thresh && is_unimodal(bootstrap_cutoffs, lambda_thresh)
  }
  
  reason <- ifelse(boot_iqr >= iqr_thresh, "High variability (IQR >= 1.5)", 
                   ifelse(!is_unimodal(bootstrap_cutoffs), "Multimodal distribution", "Stable"))
  
  return(list(stable = stable, boot_iqr = boot_iqr, reason = reason))
}


########################################
# 5. Cutoff Position Evaluation
########################################
# analyze_cutoff_position(): Comprehensive cutoff scoring based on multiple criteria
analyze_cutoff_position <- function(exprset_full, expr_tumor, gene, cutoff,
                                    tumor_values, normal_values, bootstrap_cutoffs,
                                    p_spline, min_prop = 0.25) {
  expr_tumor <- as.numeric(expr_tumor)
  bootstrap_cutoffs <- as.numeric(bootstrap_cutoffs)
  
  # 1. cutoff region
  d <- tryCatch(density(expr_tumor, na.rm = TRUE), error = function(e) NULL)
  if (is.null(d)) {
    cutoff_region <- NA
    cutoff_distance_to_peak <- NA
  } else {
    x_peak <- d$x[which.max(d$y)]
    cutoff_distance_to_peak <- abs(cutoff - x_peak)
    cutoff_region <- ifelse(
      cutoff < quantile(expr_tumor, 0.1, na.rm = TRUE), "left_tail",
      ifelse(cutoff > quantile(expr_tumor, 0.9, na.rm = TRUE), "right_tail",
             ifelse(abs(cutoff - x_peak) < 0.25, "main_peak", "plateau"))
    )
  }
  
  # 2. small group
  group <- ifelse(expr_tumor > cutoff, "high", "low")
  group_sizes <- table(group)
  small_group <- any(group_sizes / sum(group_sizes) < min_prop)
  
  # 3. TN consistency
  tumor_median <- median(tumor_values, na.rm = TRUE)
  normal_median <- median(normal_values, na.rm = TRUE)
  tn_diff <- tumor_median - normal_median
  tn_consistent <- !is.na(tn_diff) && abs(tn_diff) > 0.2
  
  # 4. Bootstrap + unimodality
  boot_iqr <- IQR(bootstrap_cutoffs, na.rm = TRUE)
  bootstrap_stable <- boot_iqr < 1.5 && is_unimodal(bootstrap_cutoffs, lambda_thresh = 0.3)
  
  # 5. Spline support
  spline_support <-  is.numeric(p_spline) && !is.na(p_spline) && p_spline < 0.05
  
  # Score (Full score: 5 points)
  score <- sum(
    as.logical(cutoff_region == "main_peak"),
    as.logical(!small_group),
    as.logical(tn_consistent),
    as.logical(bootstrap_stable),
    as.logical(spline_support),
    na.rm = TRUE
  )
  
  return(list(
    cutoff_region = cutoff_region,
    cutoff_distance_to_peak = cutoff_distance_to_peak,
    small_group = small_group,
    tn_consistent = tn_consistent,
    bootstrap_stable = bootstrap_stable,
    spline_support = spline_support,
    score = score
  ))
}

########################################
# 6. Plotting Modules
########################################
# plot_bootstrap_density(): Bootstrap cutoff density plot
plot_bootstrap_density <- function(bootstrap_cutoffs, cutoff, gene, boot_eval, plot_dir) {
  df <- data.frame(cutoff = bootstrap_cutoffs)
  p <- ggplot(df, aes(x = cutoff)) +
    geom_density(na.rm = TRUE, fill = "#69b3a2", alpha = 0.4) +
    geom_vline(xintercept = cutoff, color =  "#F46D43", linetype = "dashed", size = 1.2) +
    labs(
      title = paste("Bootstrap Cutoff Distribution -", gene),
      x = "Cutoff",
      y = "Density",
      subtitle = paste0(
        ifelse(boot_eval$stable, "Stable (Unimodal)", "Unstable (Multimodal)"),
        "; IQR = ", round(boot_eval$boot_iqr, 3)
      )
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(color = ifelse(boot_eval$stable, "#66BD63","red"), hjust = 0.5))
  
  ggsave(file.path(plot_dir, paste0(gene, "_bootstrap_density.pdf")), p, width = 7, height = 5)
}

# plot_survival_curve(): Survival plot
plot_survival_curve <- function(df, gene, plot_dir) {
  df$group <- factor(df$group, levels = c("low", "high"))
  
  surv_fit <- survfit(Surv(time, event) ~ group, data = df)
  
  cat("surv_fit content")
  print(surv_fit)
  cat("summary of surv_fit")
  summary(surv_fit)
  cat(gene,": table(df$group)")
  table(df$group)
  
  
  p <- ggsurvplot(surv_fit, data = df, pval = FALSE, conf.int = TRUE,
                  palette = c("#2E9FDF", "#E7B800"),
                  legend.title = "Group", legend.labs = c("Low", "High"),
                  ggtheme = theme_bw())
  
  # HR and P value
  sdiff <- survdiff(Surv(time, event) ~ group, data = df)
  pval <- 1 - pchisq(sdiff$chisq, df = 1)
  hr <- (sdiff$obs[2] / sdiff$exp[2]) / (sdiff$obs[1] / sdiff$exp[1])
  se_loghr <- sqrt(1 / sdiff$exp[2] + 1 / sdiff$exp[1])
  ci_low <- exp(log(hr) - qnorm(0.975) * se_loghr)
  ci_high <- exp(log(hr) + qnorm(0.975) * se_loghr)
  ci_text <- paste0(sprintf("%.3f",hr)," [",sprintf("%.3f", ci_low),", ",sprintf("%.3f",ci_high),"]")
  
  p$plot <- p$plot +
    annotate("text", x = 0, y = 0.05,
             label = paste0("P = ", sprintf("%.3f", pval), "\nHR (95% CI) = ", ci_text),
             size = 5, color = "black", hjust = 0) +
    xlab("Month") +
    theme(text = element_text(size = 15))
  
  ggsave(file.path(plot_dir, paste0(gene, "_survival_curve.pdf")), p$plot, width = 6, height = 6)
}

# plot_expression_density(): Expression density plot (Tumor vs Normal)
plot_expression_density <- function(exprset_full, gene, cutoff, plot_dir) {
  sample_type <- ifelse(as.numeric(substr(colnames(exprset_full), 14, 15)) < 10, "Tumor", "Normal")
  expr_df <- data.frame(expr = as.numeric(exprset_full[gene, ]), group = sample_type)
  
  p <- ggplot(expr_df, aes(x = expr, fill = group)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = cutoff, color = "#F46D43", linetype = "dashed", size = 1.2) +
    scale_fill_manual(values = c("Tumor" = "#2C7BB6", "Normal" = "#ABD9E9")) +
    labs(
      title = paste0("Expression distribution of ", gene),
      subtitle = paste0("Cutoff = ", signif(cutoff, 3)),
      x = paste(gene, "expression (log2 TPM)"),
      y = "Density"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave(file.path(plot_dir, paste0(gene, "_expr_density_cutoff.pdf")), p, width = 7, height = 5)
}

# plot_cutoff_radar(): Scoring radar plot
plot_cutoff_radar <- function(top_genes, plot_dir) {
  if (nrow(top_genes) == 0) {
    message("[INFO] No genes with Cutoff_Score >= 1, radar plot not created.")
    return(NULL)
  }
  
  library(fmsb)
  
  # Save Cutoff_Region
  cutoff_region_vec <- top_genes$Cutoff_Region
  
  # 5 logical indicators for scoring:
  radar_data <- top_genes[, c("Bootstrap_Stable", "Spline_Support", "TN_Consistent", "Small_Group", "Cutoff_Region")]
  
  # Logical values were converted to binary (0/1), with missing values (NA) mapped to 0 (FALSE) prior to conversion.
  radar_data <- data.frame(lapply(radar_data, function(x) {
    x[is.na(x)] <- FALSE
    as.numeric(as.logical(x))
  }))
  
  # Group balance means Small_Group here. 
  # Inversion of Small_Group variable (TRUE -> 0, FALSE -> 1)
  radar_data$Small_Group <- 1 - radar_data$Small_Group
  
  # Cutoff_Region: 1 for main_peak, 0 for others.
  radar_data$Cutoff_Region <- as.numeric(cutoff_region_vec == "main_peak")
  
  # Insert maximum/minimum row
  radar_data <- rbind(rep(1, 5), rep(0, 5), radar_data)
  rownames(radar_data) <- c("Max", "Min", make.unique(as.character(top_genes$Gene)))
  
  radar_genes <- rownames(radar_data)[-c(1, 2)]
  
  # Plotting part
  pdf(file.path(plot_dir, "cutoff_score_radar_plot.pdf"), width = 7, height = 7)
  
  colorkey <- c("#4271D6", "#F768A1", "#A50026", "#F46D43", "#FDAE61", "#FEE090",
                "#ABD9E9", "#313695", "#D9EF8B", "#66BD63", "#66C2A5")
  
  radarchart(radar_data,
             axistype = 1, plwd = 2, plty = 1,
             cglcol = "grey", cglty = 1, axislabcol = "grey30", vlcex = 0.8,
             pcol = colorkey[seq_len(length(radar_genes))],
             title = "Cutoff reliability radar (score >= 1)")
  
  legend("bottomleft", legend = radar_genes,
         col = colorkey[seq_along(radar_genes)],
         lwd = ifelse(radar_genes == "ZG16", 3, 2),
         lty = 1, cex = ifelse(radar_genes == "ZG16", 1.0, 0.8))
  
  dev.off()
}

# plot_cutoff_score_barplot(): Multi-criteria scoring barplot
plot_cutoff_score_barplot <- function(top_genes, plot_dir) {
  if (nrow(top_genes) == 0) {
    message("[INFO] No genes to plot in barplot.")
    return(NULL)
  }
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  
  
  # 1. Data Preparation + Adding Boolean and Hierarchical Features
  top_genes <- top_genes %>%
    arrange(desc(Cutoff_Score)) %>%
    mutate(
      Gene = fct_reorder(Gene, Cutoff_Score),
      Survival_Sig = case_when(
        P < 0.05 ~ "P < 0.05",
        P < 0.1 ~ "0.05 ≤ P < 0.1",
        TRUE ~ "NS"
      ),
      Peak_Region = Cutoff_Region == "main_peak",
      Group_Balance = Small_Group == FALSE
    )
  
  # Constructing Logical Boolean Point Graph Data (Reliability)
  bool_features <- top_genes %>%
    select(Gene, Peak_Region, Group_Balance, TN_Consistent, Bootstrap_Stable, Spline_Support, Survival_Sig) %>%
    pivot_longer(cols = -c(Gene, Survival_Sig), names_to = "Feature", values_to = "Value") %>%
    mutate(
      # Specify Feature Order (determines legend and vertical arrangement in the plot).
      Feature = factor(Feature, levels = c("Peak_Region", "Group_Balance", "TN_Consistent", "Bootstrap_Stable", "Spline_Support")),
      # Map feature order to Y-axis coordinates (bottom at -0.4, incrementing upwards in reverse order)
      ypos = -0.4 - (length(levels(Feature)) - as.numeric(Feature)) * 0.4  # Reverse order
    )
  
  # 3.Main figure: bar chart + scoring label
  p <- ggplot(top_genes, aes(x = Gene, y = Cutoff_Score, fill = Cutoff_Score)) +
    geom_col(width = 0.45) +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 5),name = "Cutoff Score") + 
    geom_text(aes(label = round(Cutoff_Score, 1)), vjust = -1.2, color = "black", size = 5) +
    coord_cartesian(ylim = c(-2, 5)) + coord_flip() +
    theme_minimal(base_size = 14) +
    labs(title = "Reliability Features and Cutoff Score ",
         y = "Cutoff Score (Max = 5)", x = " ") +
    theme(legend.position = "right")
  
  # 4. Point layer: Logical Boolean indicator + survival significance
  p <- p + geom_point(
    data = bool_features,
    mapping = aes(x = Gene, y = ypos, shape = Feature, color = Survival_Sig, alpha = Value),
    size = 6, stroke = 1.5, inherit.aes = FALSE
  ) +
    scale_shape_manual(values = c(
      "Peak_Region" = 18,        # ◆
      "Group_Balance" = 25,      # ▼
      "TN_Consistent" = 16,      # ●
      "Bootstrap_Stable" = 17,   # ▲
      "Spline_Support" = 15      # ■
    )) +
    scale_color_manual(values = c(
      `P < 0.05` =   "#4777EFFF",
      `0.05 ≤ P < 0.1` = "#F37651FF",
      `NS` =  "#4D4D4D"
    )) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2)) +
    guides(
      fill = guide_colorbar(title = "Cutoff Score"), 
      shape = guide_legend(title = "Reliability Features"),
      color = guide_legend(title = "Survival Significance"),
      alpha = guide_legend(title = "Satisfied")
    )
  
  # Save figure
  ggsave(file.path(plot_dir, "cutoff_score_summary_plot.pdf"), plot = p, width = 10, height = 6)
  print(p)
}


########################################
# 7. Results Summarization
########################################
# summarize_results(): Aggregate results and filter top genes
summarize_results <- function(results, adjust_method = "BH", score_threshold = 1) {
  cat(">>> Entering summarize_results()\n")
  
  valid_results <- Filter(function(x) is.data.frame(x) && "Gene" %in% colnames(x), results)
  cat("Number of valid result entries:", length(valid_results), "\n")
  
  if (length(valid_results) == 0) {
    warning("No valid results to summarize.")
    return(list(results_df = data.frame(), top_genes = data.frame()))
  }
  
  results_df <- do.call(rbind, valid_results)
  results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)
  
  numeric_cols <- c("P", "Cutoff_P", "Cutoff_Score")
  for (col in numeric_cols) {
    if (col %in% names(results_df)) {
      results_df[[col]] <- as.numeric(results_df[[col]])
    }
  }
  
  results_df$P_adj <- p.adjust(results_df$P, method = adjust_method)
  results_df$Cutoff_P_adj <- p.adjust(results_df$Cutoff_P, method = adjust_method)
  results_df$Cutoff_Score[is.na(results_df$Cutoff_Score)] <- 0
  results_df <- results_df[order(results_df$P_adj), ]
  
  top_genes <- results_df %>%
    filter(Cutoff_Score >= score_threshold)
  
  cat("Returning list with results_df and top_genes\n")
  return(list(results_df = results_df, top_genes = top_genes))
}
