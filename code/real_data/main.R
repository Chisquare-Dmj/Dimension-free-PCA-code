# Unified entry point for the EEGMMIDB and CAPTURE-24 applications.
#
# Examples:
#   Rscript code/real_data/main.R dataset=eeg action=all ncores=10
#   Rscript code/real_data/main.R dataset=capture24 action=all ncores=10
#   Rscript code/real_data/main.R dataset=sequence action=all ncores=10

find_real_project_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (!dir.exists(current)) current <- dirname(current)
  repeat {
    marker <- file.path(current, "code", "simulation", "main.R")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the project root from: ", start)
}

find_real_pipeline_root <- function() {
  configured <- Sys.getenv("DIMENSION_FREE_PCA_ROOT", unset = "")
  if (nzchar(configured) && file.exists(file.path(configured, "code", "simulation", "main.R"))) {
    return(normalizePath(configured, mustWork = TRUE))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (!identical(script_path, "-") && file.exists(script_path)) {
      return(find_real_project_root(script_path))
    }
  }
  find_real_project_root(getwd())
}

PROJECT_ROOT <- find_real_pipeline_root()
Sys.setenv(DIMENSION_FREE_PCA_ROOT = PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "code", "simulation", "main.R"))
source(file.path(PROJECT_ROOT, "code", "real_data", "common_real_data.R"))
source(file.path(PROJECT_ROOT, "code", "real_data", "eegmmidb.R"))
source(file.path(PROJECT_ROOT, "code", "real_data", "capture24.R"))

REAL_DATA_SELECTION <- "help"
REAL_DATA_ACTION <- "all"

REAL_DATA_CONFIG <- list(
  ncores = 10L,
  master_seed = MASTER_SEED,
  eeg = list(),
  capture24 = list()
)

real_data_selection <- function(selection) {
  selection <- tolower(as.character(selection))
  switch(
    selection,
    eeg = "eeg",
    eegmmidb = "eeg",
    capture24 = "capture24",
    capture = "capture24",
    `capture-24` = "capture24",
    all = c("eeg", "capture24"),
    sequence = c("eeg", "capture24"),
    help = "help",
    stop("Unknown real-data selection: ", selection)
  )
}

configured_real_data <- function(dataset, config) {
  defaults <- switch(
    dataset,
    eeg = default_eegmmidb_config(),
    capture24 = default_capture24_config()
  )
  global <- config[intersect(names(config), names(defaults))]
  global <- global[!vapply(global, is.list, logical(1))]
  local <- config[[dataset]]
  if (is.null(local)) local <- list()
  merge_config(merge_config(defaults, global), local)
}

regenerate_real_data_table <- function(dataset, root = PROJECT_ROOT, run_id = NULL) {
  manifest_dataset <- switch(
    dataset, eeg = "eegmmidb", capture24 = "capture24"
  )
  inference_path <- resolve_real_data_run(manifest_dataset, "inference", root, run_id)
  table_path <- resolve_real_data_run(manifest_dataset, "table", root, run_id)
  optional_separation_artifacts <- c("k0_sensitivity", "domain_expansion")
  for (artifact in optional_separation_artifacts) {
    path <- tryCatch(
      resolve_real_data_run(manifest_dataset, artifact, root, run_id),
      error = function(error) NULL
    )
    if (!is.null(path)) refresh_saved_separation_artifact(path)
  }
  inference <- refresh_saved_separation_artifact(inference_path)
  label <- switch(
    dataset,
    eeg = "EEGMMIDB R04",
    capture24 = "CAPTURE-24 wake-aligned ENMO"
  )
  write_real_inference_table(inference, table_path, label)
  if (dataset == "eeg") {
    comparison <- eeg_lateralization_comparison(inference, 2L)
    prefix <- sub("__spectral_inference[.]csv$", "", inference_path)
    artifacts <- list(
      fpca_comparison = paste0(prefix, "__pc2_lateralization_fpca_vs_proposed.csv")
    )
    write_csv_atomic(comparison, artifacts$fpca_comparison)
    write_csv_atomic(comparison, file.path(dirname(inference_path), "eeg_fpca_vs_proposed.csv"))
    write_csv_atomic(comparison,
                     file.path(dirname(inference_path), "eeg_pc2_lateralization_fpca_vs_proposed.csv"))
    augment_real_data_run_artifacts(manifest_dataset,
                                    resolved_real_data_run_id(manifest_dataset, root, run_id),
                                    artifacts, root)
  }
  if (dataset == "capture24") {
    pilot_path <- resolve_real_data_run(manifest_dataset, "pilot_summary", root, run_id)
    pilot_table_path <- resolve_real_data_run(manifest_dataset, "pilot_table", root, run_id)
    write_capture24_pilot_table(read_required_csv(pilot_path), pilot_table_path)
    refresh_capture24_supplementary_outputs(root, run_id)
  }
  invisible(table_path)
}

run_one_real_data <- function(dataset, action, config) {
  action <- match.arg(tolower(action), c("all", "download", "preprocess", "analyze", "summarize", "plot"))
  cfg <- configured_real_data(dataset, config)
  started <- Sys.time()
  pipeline_log(paste("REAL DATA", toupper(dataset), "|", toupper(action), "| START"))
  if (action %in% c("all", "download")) {
    download_started <- Sys.time()
    switch(
      dataset,
      eeg = download_eegmmidb(cfg), capture24 = download_capture24(cfg)
    )
    pipeline_log(paste("REAL DATA", toupper(dataset), "| DOWNLOAD COMPLETE | elapsed=", format_elapsed(download_started)))
  }
  if (action == "download") return(invisible(NULL))
  if (action == "preprocess") {
    switch(
      dataset,
      eeg = prepare_eegmmidb(cfg), capture24 = prepare_capture24(cfg)
    )
    pipeline_log(paste("REAL DATA", toupper(dataset), "| PREPROCESS COMPLETE | elapsed=", format_elapsed(started)))
    return(invisible(NULL))
  }
  if (action %in% c("all", "analyze")) {
    switch(
      dataset,
      eeg = run_eegmmidb_analysis(cfg), capture24 = run_capture24_analysis(cfg)
    )
  } else if (action == "summarize") {
    regenerate_real_data_table(dataset, cfg$root)
  } else if (action == "plot") {
    switch(dataset, eeg = plot_eegmmidb_results(cfg$root),
           capture24 = plot_capture24_results(cfg$root))
  }
  pipeline_log(paste("REAL DATA", toupper(dataset), "| COMPLETE | elapsed=", format_elapsed(started)))
  invisible(NULL)
}

run_real_data_pipeline <- function(selection = REAL_DATA_SELECTION, action = REAL_DATA_ACTION,
                                   config = REAL_DATA_CONFIG) {
  datasets <- real_data_selection(selection)
  if (identical(datasets, "help")) return(invisible(print_real_data_help()))
  for (dataset in datasets) run_one_real_data(dataset, action, config)
  invisible(NULL)
}

print_real_data_help <- function() {
  cat(
    "Dimension-Free PCA real-data pipeline\n\n",
    "Selections: dataset=eeg, capture24, all, or sequence\n",
    "Actions:    action=all, download, preprocess, analyze, summarize, or plot\n\n",
    "Examples:\n",
    "  Rscript code/real_data/main.R dataset=eeg action=all ncores=10\n",
    "  Rscript code/real_data/main.R dataset=capture24 action=all ncores=10\n",
    "  Rscript code/real_data/main.R dataset=sequence action=all ncores=10\n",
    "  Rscript code/real_data/main.R dataset=eeg action=plot\n",
    sep = ""
  )
}

if (sys.nframe() == 0L) {
  cli <- parse_cli_overrides(commandArgs(trailingOnly = TRUE))
  selection <- if (!is.null(cli$dataset)) cli$dataset else REAL_DATA_SELECTION
  action <- if (!is.null(cli$action)) cli$action else REAL_DATA_ACTION
  cli$dataset <- NULL; cli$action <- NULL
  run_real_data_pipeline(selection, action, merge_config(REAL_DATA_CONFIG, cli))
}
