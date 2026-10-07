# Build the report from a fresh Rscript session or the open RStudio project.
main <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  script <- grep("^--file=", arguments, value = TRUE)
  root <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[1]))) else getwd()
  if (!file.exists(file.path(root, "analysis.Rmd"))) {
    stop("Open climate-seasonality-trends.Rproj or run this script from the project folder.")
  }
  previous <- setwd(root)
  on.exit(setwd(previous), add = TRUE)
  needed <- c("knitr", "maps")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Run setup.R first. Missing packages: ", paste(missing, collapse = ", "))
  if (!file.exists("data/raw/uscrn_daily.RData")) {
    stop("Place uscrn_daily.RData in data/raw/. See data/README.md.")
  }
  for (folder in c("figures", "results")) dir.create(folder, showWarnings = FALSE)
  knitr::knit("analysis.Rmd", output = "analysis.md", envir = new.env(parent = globalenv()), quiet = FALSE)
  # The source retains RStudio's output setting; the GitHub report starts at its title.
  report <- readLines("analysis.md", warn = FALSE)
  if (length(report) && report[1] == "---") {
    boundary <- which(report == "---")[2]
    if (!is.na(boundary)) writeLines(report[-seq_len(boundary)], "analysis.md")
  }
  message("Completed: analysis.md, figures/, and results/.")
}
main()
