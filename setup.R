# Run once to install the packages used for analysis, rendering, and numerical checks.
packages <- c("knitr", "maps", "sandwich")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Packages are ready. Place the RData file in data/raw/, then run run_analysis.R.")
