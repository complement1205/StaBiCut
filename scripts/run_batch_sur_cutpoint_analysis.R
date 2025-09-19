#' Stability-Informed Batch Survival Cutoff Analysis
#'
#' Performs cutoff optimization for multiple genes using expression and survival data,
#' integrating Cox regression, spline support, bootstrap stability, tumor–normal consistency,
#' and multi-criteria scoring with visualization.

#' It integrates multiple evaluation modules, including:
#'   1. Optimal cutoff detection (Cox regression)
#'   2. Spline-based model support
#'   3. Bootstrap cutoff distribution and stability
#'   4. Tumor–normal consistency check
#'   5. Multi-criteria scoring and visualization
#'
#' @param exprset   Expression matrix (genes x samples, log-transformed inside).
#' @param geneset   Vector of candidate gene symbols/IDs.
#' @param clin      Clinical data (must contain survival time and event).
#' @param n_boot    Number of bootstrap iterations (default = 1000).
#' @param adjust_method  Multiple testing correction method (default = "BH").
#' @param cut.type  Cutoff selection method. Options include "optimal" (default), "median", or "quantile".
#' @param save_plots Logical; save output plots to `plot_dir` (default = TRUE).
#' @param score_threshold Minimum score threshold for gene prioritization.
#' @param plot_dir  Output directory for generated plots (default = "./output").
#'
#' @return A data frame summarizing hazard ratio, confidence intervals, p-values,
#'         cutoff statistics, bootstrap stability metrics, tumor–normal consistency,
#'         spline support, and integrated cutoff scores for each gene.

#' @export 
run_batch_sur_cutpoint_analysis <- function(exprset, geneset, clin, n_boot = 1000,adjust_method = "BH",
                                            cut.type = "optimal",save_plots = TRUE,score_threshold = 1,
                                            plot_dir = "./output") {
  # 1. Initialization
  if (!dir.exists(plot_dir) & save_plots) dir.create(plot_dir, recursive = TRUE)
  results <- list()
  
  exprset <- log2(exprset + 1)
  
  # Save original expression matrix (including normal samples)
  exprset_full <- exprset
  
  is_tumor <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
  exprset <- exprset[, is_tumor]
  clin <- clin[is_tumor, ]
  
  # Add survival data (time and event) 
  clin$time <- ifelse(!is.na(clin$days_to_death),
                      clin$days_to_death,
                      clin$days_to_last_follow_up) / 30
  clin$event <- ifelse(clin$vital_status %in% c("Dead", "dead", "DEAD"), 1, 0)
  clin <- clin[!is.na(clin$time) & !is.na(clin$event), ]
  
  # Match samples between expression and clinical data
  sample_ids <- intersect(colnames(exprset), rownames(clin))
  exprset <- exprset[, sample_ids]
  clin <- clin[sample_ids, ]
  cat("Number of valid tumor samples:", nrow(clin), "\n")
  
  
  # 2. Iterate over candidate genes
  for (gene in geneset) {
    cat("Processing gene:", gene, "\n")
    expr <- exprset[gene, ]
    
    if (length(expr) != nrow(clin)) {
      cat(" Expression and clinical data have mismatched sample counts, skipping\n")
      next
    }
    
    df <- data.frame(expr = as.numeric(expr), time =as.numeric(clin$time), event = as.numeric(clin$event))
  
    
    # 3. Determine optimal cutoff (Cox-based)
    result <- compute_cutoff_and_cox(df, gene = gene, cut.type = "optimal", minprop = 0.25)
    
    df <- result$df
    cutoff <- result$cutoff
    cutoff_p <- result$cutoff_p
    cox <- list(
      hr = result$hr,
      ci_low = result$ci_low,
      ci_high = result$ci_high,
      p_val = result$p_val
    )
    
    
    # 4. Spline-based model support
    spline_result <- evaluate_spline_support(df=df, gene=gene,plot_dir = plot_dir)
    p_spline <- spline_result$p_spline
    spline_support <- spline_result$spline_support
    
    
    # 5. Bootstrap cutoff distribution
    bootstrap_cutoffs <- compute_bootstrap_cutoffs(df=df, n_boot = 1000, minprop = 0.25)
    cat(gene, ": NA bootstrap count = ", sum(is.na(bootstrap_cutoffs)), "\n")
    
    boot_mean <- mean(bootstrap_cutoffs, na.rm = TRUE)
    boot_median <- median(bootstrap_cutoffs, na.rm = TRUE)
    boot_ci <- quantile(bootstrap_cutoffs, probs = c(0.025, 0.975), na.rm = TRUE)
    
    # Assess bootstrap stability: a large IQR indicates an unstable cutoff.
    boot_iqr <- IQR(bootstrap_cutoffs, na.rm = TRUE)
    
    
    # 6. Bootstrap stability evaluation
    boot_eval <- evaluate_bootstrap_stability(bootstrap_cutoffs=bootstrap_cutoffs, iqr_thresh = 1.5, lambda_thresh = 0.3)
    
    stable <- boot_eval$stable
    Bootstrap_IQR <-  boot_eval$boot_iqr
    
    
    # 7. Tumor–normal consistency check
    expr_full <- exprset_full[gene, ]
    sample_type_full <- as.numeric(substr(colnames(expr_full), 14, 15))
    tumor_expr <- as.numeric(expr_full[sample_type_full < 10])
    normal_expr <- as.numeric(expr_full[sample_type_full >= 10])
    
    
    # 8. Integrated cutoff evaluation (multi-criteria scoring)
    cutoff_eval <- analyze_cutoff_position(exprset_full = exprset_full,
                                           expr_tumor = as.numeric(expr),
                                           cutoff = cutoff,
                                           tumor_values = tumor_expr,
                                           normal_values = normal_expr,
                                           bootstrap_cutoffs = bootstrap_cutoffs,
                                           p_spline = p_spline,
                                           min_prop = 0.25)
    
    
    # 9. Aggregate results for the gene
    results[[gene]] <- data.frame(
      Gene = gene,
      HR = cox$hr,
      CI_low = cox$ci_low,
      CI_high = cox$ci_high,
      P = cox$p_val,
      Cutoff = cutoff,
      p_spline = signif(p_spline, 3),
      Bootstrap_Mean = boot_mean,
      Bootstrap_Median = boot_median,
      Bootstrap_CI_low = boot_ci[1],
      Bootstrap_CI_high = boot_ci[2],
      Cutoff_P = cutoff_p,
      Bootstrap_IQR = boot_iqr,
      Stability_Class = ifelse(boot_eval$stable, "Stable", "Unstable"),
      Stability_Reason = boot_eval$reason,
      Cutoff_Region = cutoff_eval$cutoff_region,
      Small_Group = cutoff_eval$small_group,
      TN_Consistent = cutoff_eval$tn_consistent,
      Bootstrap_Stable = cutoff_eval$bootstrap_stable,
      Spline_Support = cutoff_eval$spline_support,
      Cutoff_Score = cutoff_eval$score,
      stringsAsFactors = FALSE
    )
  }
  
  
  # 10. Summarize across genes
  # Screen for candidate genes
  summary <- summarize_results(results, adjust_method = "BH", score_threshold = 1)
  results_df <- summary$results_df
  top_genes <- summary$top_genes
  
  
  # 11. Visualization
  plot_bootstrap_density(bootstrap_cutoffs, cutoff, gene, boot_eval, plot_dir)
  plot_survival_curve(df, gene, plot_dir)
  plot_cutoff_radar(top_genes, plot_dir)
  plot_cutoff_score_barplot(top_genes, plot_dir)   
  
  cat(">>> Ready to return results_df\n")
  return(results_df)
}

