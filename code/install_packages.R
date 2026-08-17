# Install the small set of non-base R packages used by the reproducibility code.

packages <- c("edfReader", "signal", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Required R packages are available: ", paste(packages, collapse = ", "))
