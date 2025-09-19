# Install and load required R packages for StaBiCut

packages <- c(
  "survival",
  "survminer",
  "splines",
  "ggplot2",
  "mixtools",
  "dplyr",
  "tidyr",
  "fmsb",
  "openxlsx",
  "forcats"
)

# Check and install missing packages
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

