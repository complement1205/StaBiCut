##########################################################
# StaBiCut - Stability-informed Survival Cutoff Analysis
# main.R
# 
# Entry point script to run the StaBiCut pipeline on example data.
# Loads dependencies, sources analysis functions, and executes workflow.
#
# Example objects provided in data/example_expr_clinical.Rdata:
#   - mrna_expr_tpm (expression matrix)
#   - clinicalSE   (clinical data)
#
# Author: Lu Chen

##########################################################


# Load required packages and scripts
source("requirements.R")
source("scripts/modules.R")
source("scripts/run_batch_sur_cutpoint_analysis.R")

# Load example data
load("data/example_expr_clinical.Rdata")
# Example objects provided: mrna_expr_tpm, clinicalSE

# Define gene set (can also be read from file: e.g., readLines("data/geneset.txt"))
geneset <- c("SPOCK2","PYCR1","CA4","CES1","ABCB1",
             "ZG16","TNXB","HMCN2","MEP1A","SLC37A2","CHGB")

# Ensure output directory exists
if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# Pipeline configuration
config <- list(
  exprset = mrna_expr_tpm,
  geneset = geneset,
  clin = clinicalSE,
  n_boot = 1000,
  adjust_method = "BH",
  cut.type = "optimal",
  save_plots = TRUE,
  score_threshold = 1,
  plot_dir = "output"
)


# Run analysis
all_results <- do.call(run_batch_sur_cutpoint_analysis, config)

# Save results
write.csv(all_results,
          file = "output/Cutoff_results.csv",
          row.names = FALSE)

cat("StaBiCut analysis completed. Results saved in output/Cutoff_results.csv\n")

