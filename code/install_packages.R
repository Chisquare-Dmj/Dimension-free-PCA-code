# Install the small set of non-base R packages used by the reproducibility code.
# Recorded computation-server versions are listed in software_versions.txt.
# This helper installs current CRAN releases when packages are missing; it is not a lockfile.

packages <- c("edfReader", "signal", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Required R packages are available: ", paste(packages, collapse = ", "))
