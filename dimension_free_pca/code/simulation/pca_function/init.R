# Central source hub. Keep dependency order explicit because later files use
# helpers defined by earlier files.

if (!exists("PROJECT_ROOT", inherits = FALSE)) {
  PROJECT_ROOT <- normalizePath(
    Sys.getenv("DIMENSION_FREE_PCA_ROOT", unset = getwd()),
    winslash = "/",
    mustWork = FALSE
  )
}
SIMULATION_CODE_ROOT <- file.path(PROJECT_ROOT, "code", "simulation", "pca_function")
source(file.path(SIMULATION_CODE_ROOT, "common_utils.R"))
source(file.path(SIMULATION_CODE_ROOT, "inference_utils.R"))
source(file.path(SIMULATION_CODE_ROOT, "model_builders.R"))
source(file.path(SIMULATION_CODE_ROOT, "experiment_phase.R"))
source(file.path(SIMULATION_CODE_ROOT, "experiment_functional.R"))
source(file.path(SIMULATION_CODE_ROOT, "experiment_functional_asymptotics.R"))
source(file.path(SIMULATION_CODE_ROOT, "experiment_large_domain.R"))
source(file.path(SIMULATION_CODE_ROOT, "experiment_robustness.R"))
source(file.path(SIMULATION_CODE_ROOT, "summarize_results.R"))
source(file.path(SIMULATION_CODE_ROOT, "plot_results.R"))
