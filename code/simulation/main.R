# Unified pipeline for the dimension-free PCA inference experiments.
#
# The selection, action, and configuration objects below are the normal user
# interface. Select "sequence" for a complete experiment-by-experiment run,
# "all" for the batch workflow, a numbered group ("1", "2", "3"), or one
# branch ("1a", "1b", "2a", "2b", "3a", "3b", "3c", "3d", "k0", "universality").

PIPELINE_MODE <- "help"
PIPELINE_EXPERIMENT <- "help"
PIPELINE_ACTION <- "all"

# The sequence runner completes every experiment bundle before starting the
# next one. Edit this vector only when a different execution order is needed.
PIPELINE_SEQUENCE <- c("1a", "1b", "2a", "2b", "3a", "3b", "3c", "3d", "k0", "universality")

PIPELINE_CONFIG <- list(
  # Number of parallel forked workers used for Monte Carlo replications.
  # Use 1 for serial execution and debugging. Large production runs should use
  # a value compatible with both the available CPU count and memory capacity.
  ncores = 100L,

  # Replication seeds equal master_seed + 100000 * experiment_id +
  # replication_id. Population-design seeds are separate and fixed.
  master_seed = 20260810L,

  # All intervals use this confidence level. The default 0.95 gives the
  # two-sided z_0.975 Wald interval and one-sided z_0.95 phase lower bound.
  confidence_level = 0.95,

  # Per-experiment overrides. Empty lists use the documented defaults. Add only
  # fields that should change, for example experiment_1a = list(replications=50L).
  experiment_1a = list(),
  experiment_1b = list(),
  experiment_2 = list(),
  experiment_2b = list(),
  experiment_3a = list(),
  experiment_3b = list(),
  experiment_3c = list(),
  experiment_3d = list(),
  robustness_k0 = list(),
  robustness_universality = list()
)

find_project_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (!dir.exists(current)) current <- dirname(current)
  repeat {
    marker <- file.path(current, "code", "simulation", "pca_function", "init.R")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the project root from: ", start)
}

find_script_root <- function() {
  configured_root <- Sys.getenv("DIMENSION_FREE_PCA_ROOT", unset = "")
  if (nzchar(configured_root) && file.exists(file.path(configured_root, "code", "simulation", "pca_function", "init.R"))) {
    return(normalizePath(configured_root, mustWork = TRUE))
  }
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (!identical(script_path, "-") && file.exists(script_path)) {
      return(find_project_root(script_path))
    }
  }
  find_project_root(getwd())
}

PROJECT_ROOT <- find_script_root()
Sys.setenv(DIMENSION_FREE_PCA_ROOT = PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "code", "simulation", "pca_function", "init.R"))
initialize_output_dirs(PROJECT_ROOT)

experiment_catalog <- function() {
  identity_config <- function(config) config
  list(
    `1a` = list(
      experiment = "experiment_1a_panel_phase", config_key = "experiment_1a",
      defaults = default_panel_phase_config, runner = run_panel_phase_experiment,
      prepare = identity_config,
      id_fields = c("n", "p", "T", "replications", "M", "K0", "target_Delta3", "distribution")
    ),
    `1b` = list(
      experiment = "experiment_1b_phase_margin", config_key = "experiment_1b",
      defaults = default_phase_margin_config, runner = run_phase_margin_experiment,
      prepare = prepare_phase_margin_config,
      id_fields = c("n", "p", "replications", "M", "K0", "Delta_grid_main", "Delta_grid_boundary", "distribution")
    ),
    `2a` = list(
      experiment = "experiment_2_functional_inference", config_key = "experiment_2",
      defaults = default_functional_inference_config, runner = run_functional_inference_experiment,
      prepare = function(config) { config$design_grid <- paste0(length(config$scenarios), "Sx", length(config$distributions), "D"); config },
      id_fields = c("n", "p", "N", "decay", "M", "K0", "core_replications", "other_replications", "design_grid")
    ),
    `2b` = list(
      experiment = "experiment_2b_functional_asymptotics", config_key = "experiment_2b",
      defaults = default_functional_asymptotic_config, runner = run_functional_asymptotic_experiment,
      prepare = prepare_functional_asymptotic_config,
      id_fields = c("n_values", "replications", "p_rule", "N", "decay", "M", "K0", "spike_indices", "design_grid")
    ),
    `3a` = list(
      experiment = "experiment_3a_simple_direction", config_key = "experiment_3a",
      defaults = default_large_domain_config, runner = run_simple_direction_experiment,
      prepare = function(config) { config$domain_rule <- "T-n-div-3"; config$basis_rule <- "N-2n"; config },
      id_fields = c("n_values", "replications", "domain_rule", "basis_rule", "decay", "M", "K0", "distribution")
    ),
    `3b` = list(
      experiment = "experiment_3b_repeated_eigenspace", config_key = "experiment_3b",
      defaults = default_large_domain_config, runner = run_repeated_eigenspace_experiment,
      prepare = function(config) { config$domain_rule <- "T-n-div-3"; config$basis_rule <- "N-2n"; config },
      id_fields = c("n_values", "replications", "domain_rule", "basis_rule", "decay", "M", "K0", "distribution")
    ),
    `3c` = list(
      experiment = "experiment_3c_grid_saturation", config_key = "experiment_3c",
      defaults = default_grid_saturation_config, runner = run_grid_saturation_experiment,
      prepare = identity_config,
      id_fields = c("n", "T", "N", "m_values", "replications", "decay", "M", "K0", "distribution")
    ),
    `3d` = list(
      experiment = "experiment_3d_eigengap_inference", config_key = "experiment_3d",
      defaults = default_eigengap_config, runner = run_eigengap_experiment,
      prepare = function(config) { config$domain_rule <- "T-n-div-3"; config$basis_rule <- "N-2n"; config },
      id_fields = c(
        "n_values", "replications", "delta_grid", "pair_alpha", "domain_rule",
        "basis_rule", "decay", "M", "K0", "distribution", "gap_definition",
        "fpca_equality_draws", "fpca_equality_seed"
      )
    ),
    k0 = list(
      experiment = "robustness_k0", config_key = "robustness_k0",
      defaults = default_k0_robustness_config, runner = run_k0_robustness_experiment,
      prepare = prepare_k0_robustness_config,
      id_fields = c(
        "n_values", "replications", "p_rule", "N", "M", "K0_values",
        "spike_indices", "scenario", "distribution_grid"
      )
    ),
    universality = list(
      experiment = "robustness_universality", config_key = "robustness_universality",
      defaults = default_universality_config, runner = run_universality_stress_experiment,
      prepare = function(config) { config$distribution_grid <- paste(config$distributions, collapse = "-"); config },
      id_fields = c("n", "p", "N", "replications", "M", "K0", "scenario", "distribution_grid")
    )
  )
}

expand_experiment_selection <- function(selection) {
  selection <- tolower(as.character(selection))
  aliases <- list(
    all = names(experiment_catalog()),
    sequence = PIPELINE_SEQUENCE,
    `1` = c("1a", "1b"), `1a` = "1a", `1b` = "1b",
    `2` = c("2a", "2b"), `2a` = "2a", `2b` = "2b",
    `3` = c("3a", "3b", "3c", "3d"),
    `3a` = "3a", `3b` = "3b", `3c` = "3c", `3d` = "3d",
    robustness = c("k0", "universality"), k0 = "k0",
    universality = "universality"
  )
  result <- aliases[[selection]]
  if (is.null(result)) stop("Unknown experiment selection: ", selection)
  result
}

pipeline_log <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), message))
  flush(stdout())
}

format_elapsed <- function(started_at) {
  seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  sprintf("%.1f seconds", seconds)
}

run_complete_experiment_sequence <- function(config = PIPELINE_CONFIG,
                                             sequence = PIPELINE_SEQUENCE) {
  catalog <- experiment_catalog()
  unknown <- setdiff(sequence, names(catalog))
  if (length(unknown)) stop("Unknown experiment(s) in PIPELINE_SEQUENCE: ", paste(unknown, collapse = ", "))
  root <- if (!is.null(config$root)) config$root else PROJECT_ROOT
  initialize_output_dirs(root)
  run_ids <- list()
  sequence_started_at <- Sys.time()
  pipeline_log(sprintf(
    "SEQUENCE START | experiments=%d | order=%s | ncores=%s",
    length(sequence), paste(sequence, collapse = " -> "), config$ncores
  ))
  for (index in seq_along(sequence)) {
    short_id <- sequence[index]
    specification <- catalog[[short_id]]
    experiment_started_at <- Sys.time()
    progress_label <- sprintf("EXPERIMENT %d/%d | selector=%s | name=%s",
                              index, length(sequence), short_id, specification$experiment)
    pipeline_log(paste(progress_label, "| START"))
    tryCatch({
      simulation_started_at <- Sys.time()
      pipeline_log(paste(progress_label, "| SIMULATION START"))
      specification$runner(configured_experiment(short_id, config))
      run_id <- resolved_run_id(specification$experiment, root)
      run_ids[[specification$experiment]] <- run_id
      pipeline_log(sprintf(
        "%s | SIMULATION COMPLETE | elapsed=%s | run_id=%s",
        progress_label, format_elapsed(simulation_started_at), run_id
      ))

      summary_started_at <- Sys.time()
      pipeline_log(paste(progress_label, "| SUMMARY AND TABLE START"))
      summary_output <- generate_experiment_summary(specification$experiment, root, run_id)
      pipeline_log(sprintf(
        "%s | SUMMARY AND TABLE COMPLETE | elapsed=%s | csv=%s | table=%s",
        progress_label, format_elapsed(summary_started_at), relative_to_root(summary_output$csv, root),
        relative_to_root(summary_output$table, root)
      ))

      figure_started_at <- Sys.time()
      pipeline_log(paste(progress_label, "| FIGURE START"))
      figure_output <- plot_experiment_results(specification$experiment, root, run_id)
      figure_files <- unlist(figure_output, recursive = TRUE, use.names = FALSE)
      figure_files <- unique(figure_files[file.exists(figure_files)])
      pipeline_log(sprintf(
        "%s | FIGURE COMPLETE | elapsed=%s | files=%s",
        progress_label, format_elapsed(figure_started_at),
        paste(vapply(figure_files, relative_to_root, character(1), root = root), collapse = ", ")
      ))
      pipeline_log(sprintf("%s | COMPLETE | elapsed=%s", progress_label, format_elapsed(experiment_started_at)))
    }, error = function(error) {
      pipeline_log(sprintf("%s | FAILED | %s", progress_label, conditionMessage(error)))
      stop(conditionMessage(error), call. = FALSE)
    })
  }
  experiments <- vapply(sequence, function(id) catalog[[id]]$experiment, character(1))
  if (all(c("experiment_1a_panel_phase", "experiment_1b_phase_margin") %in% experiments)) {
    figure_set_started_at <- Sys.time()
    pipeline_log("PHASE AND CONDITIONING FIGURE SET | START")
    plot_figure_2(root, run_ids)
    pipeline_log(sprintf("PHASE AND CONDITIONING FIGURE SET | COMPLETE | elapsed=%s", format_elapsed(figure_set_started_at)))
  }
  if (all(c("experiment_3a_simple_direction", "experiment_3b_repeated_eigenspace", "experiment_3c_grid_saturation") %in% experiments)) {
    figure_set_started_at <- Sys.time()
    pipeline_log("LARGE-DOMAIN FIGURE SET | START")
    plot_figure_4(root, run_ids)
    pipeline_log(sprintf("LARGE-DOMAIN FIGURE SET | COMPLETE | elapsed=%s", format_elapsed(figure_set_started_at)))
  }
  if (all(c("robustness_k0", "robustness_universality") %in% experiments)) {
    figure_set_started_at <- Sys.time()
    pipeline_log("ROBUSTNESS FIGURE SET | START")
    plot_robustness_results(root, run_ids)
    pipeline_log(sprintf("ROBUSTNESS FIGURE SET | COMPLETE | elapsed=%s", format_elapsed(figure_set_started_at)))
  }
  pipeline_log(sprintf("SEQUENCE COMPLETE | elapsed=%s", format_elapsed(sequence_started_at)))
  invisible(run_ids)
}

configured_experiment <- function(short_id, config) {
  specification <- experiment_catalog()[[short_id]]
  configured <- experiment_config(specification$defaults(), config, specification$config_key)
  requested_replications <- config[["replications"]]
  local_replications <- config[[specification$config_key]][["replications"]]
  if (!is.null(local_replications)) requested_replications <- local_replications
  if (short_id == "2a" && !is.null(requested_replications)) {
    configured$core_replications <- requested_replications
    configured$other_replications <- requested_replications
  }
  if (short_id %in% c("2b", "3a", "3b", "3d", "k0")) {
    if (length(configured$replications) == 1L) {
      configured$replications <- setNames(
        rep(configured$replications, length(configured$n_values)), configured$n_values
      )
    } else if (!is.null(names(configured$replications)) &&
               all(as.character(configured$n_values) %in% names(configured$replications))) {
      configured$replications <- configured$replications[as.character(configured$n_values)]
    }
  }
  requested_K0 <- config[["K0"]]
  local_K0 <- config[[specification$config_key]][["K0"]]
  if (!is.null(local_K0)) requested_K0 <- local_K0
  if (short_id == "k0" && !is.null(requested_K0)) configured$K0_values <- requested_K0
  if (short_id == "k0") {
    requested_distribution <- config[["distribution"]]
    local_distribution <- config[[specification$config_key]][["distribution"]]
    if (!is.null(local_distribution)) requested_distribution <- local_distribution
    if (!is.null(requested_distribution)) configured$distributions <- requested_distribution
  }
  specification$prepare(configured)
}

configured_run_id <- function(short_id, config) {
  specification <- experiment_catalog()[[short_id]]
  experiment_configured <- configured_experiment(short_id, config)
  make_run_id(specification$experiment, experiment_configured, specification$id_fields)
}

run_selected_experiments <- function(selection, action = "all", config = PIPELINE_CONFIG,
                                     exact_saved_config = FALSE) {
  action <- match.arg(tolower(action), c("all", "simulate", "summarize", "plot"))
  if (tolower(as.character(selection)) == "sequence") {
    if (action != "all") stop("experiment=sequence is a complete-run setting and requires action=all.")
    return(run_complete_experiment_sequence(config))
  }
  short_ids <- expand_experiment_selection(selection)
  catalog <- experiment_catalog()
  experiments <- vapply(short_ids, function(id) catalog[[id]]$experiment, character(1))
  root <- if (!is.null(config$root)) config$root else PROJECT_ROOT
  initialize_output_dirs(root)

  if (action %in% c("all", "simulate")) {
    for (id in short_ids) catalog[[id]]$runner(configured_experiment(id, config))
  }
  if (action == "simulate") return(invisible(NULL))

  run_ids <- if (!is.null(config$result_run_ids)) config$result_run_ids else list()
  if (action == "all") {
    for (experiment in experiments) run_ids[[experiment]] <- resolved_run_id(experiment, root)
  } else if (exact_saved_config) {
    for (id in short_ids) {
      experiment <- catalog[[id]]$experiment
      if (is.null(run_ids[[experiment]])) run_ids[[experiment]] <- configured_run_id(id, config)
    }
  }
  if (action %in% c("all", "summarize")) {
    generate_all_summaries(root, run_ids, experiments)
    fpca_experiments <- c(
      "experiment_2b_functional_asymptotics", "experiment_3a_simple_direction",
      "experiment_3b_repeated_eigenspace"
    )
    if (any(fpca_experiments %in% experiments)) {
      generate_fpca_comparison_summaries(root, run_ids)
    }
  }
  if (action %in% c("all", "plot")) {
    plot_selected_results(experiments, root, run_ids, include_comparison_exports = TRUE)
  }
  invisible(run_ids)
}

pipeline_modes <- function() {
  c(
    "help",
    "simulate_1a_panel_phase", "simulate_1b_phase_margin",
    "simulate_2_functional_inference", "simulate_2b_functional_asymptotics",
    "simulate_3a_simple_direction", "simulate_3b_repeated_eigenspace",
    "simulate_3c_grid_saturation", "simulate_3d_eigengap_inference",
    "simulate_robustness_k0", "simulate_robustness_universality",
    "simulate_all", "summarize_all",
    "plot_figure_2", "plot_figure_3", "plot_figure_4", "plot_robustness",
    "plot_all", "all", "smoke"
  )
}

print_pipeline_help <- function() {
  cat("Dimension-free PCA inference pipeline\n\n")
  cat("Recommended interface:\n")
  cat("  Rscript code/simulation/main.R experiment=1a\n")
  cat("  Rscript code/simulation/main.R experiment=1b n=200 p=200 replications=500\n")
  cat("  Rscript code/simulation/main.R experiment=3 action=plot\n")
  cat("  Rscript code/simulation/main.R experiment=all ncores=20\n")
  cat("  Rscript code/simulation/main.R experiment=sequence ncores=20\n\n")
  cat("Actions: all (simulate + summarize + plot), simulate, summarize, plot.\n")
  cat("Selections: sequence, all, 1, 1a, 1b, 2, 2a, 2b, 3, 3a, 3b, 3c, 3d, robustness, k0, universality.\n\n")
  cat("The sequence selection completes each experiment's full output bundle before starting the next one.\n\n")
  cat("Legacy mode interface remains available.\n")
  cat("Available modes:\n")
  cat(paste0("  - ", pipeline_modes(), collapse = "\n"), "\n")
  cat("\nSimulation modes save replicate-level CSV files. Plot modes only read saved results.\n")
}

run_all_simulations <- function(config = list()) {
  run_panel_phase_experiment(experiment_config(default_panel_phase_config(), config, "experiment_1a"))
  run_phase_margin_experiment(experiment_config(default_phase_margin_config(), config, "experiment_1b"))
  run_functional_inference_experiment(experiment_config(default_functional_inference_config(), config, "experiment_2"))
  run_functional_asymptotic_experiment(experiment_config(default_functional_asymptotic_config(), config, "experiment_2b"))
  run_simple_direction_experiment(experiment_config(default_large_domain_config(), config, "experiment_3a"))
  run_repeated_eigenspace_experiment(experiment_config(default_large_domain_config(), config, "experiment_3b"))
  run_grid_saturation_experiment(experiment_config(default_grid_saturation_config(), config, "experiment_3c"))
  run_eigengap_experiment(experiment_config(default_eigengap_config(), config, "experiment_3d"))
  run_k0_robustness_experiment(experiment_config(default_k0_robustness_config(), config, "robustness_k0"))
  run_universality_stress_experiment(experiment_config(default_universality_config(), config, "robustness_universality"))
  invisible(NULL)
}

smoke_config <- function(config = list()) {
  merge_config(config, list(
    root = tempfile("dimension-free-pca-smoke-", tmpdir = tempdir()),
    ncores = min(2L, parallel::detectCores()),
    experiment_1a = list(n = 30L, p = 4L, T = 10L, replications = 4L),
    experiment_1b = list(n = 30L, p = 30L, Delta_grid = c(0.30, 0.60), replications = 4L),
    experiment_2 = list(n = 20L, p = 10L, N = 12L, core_replications = 3L, other_replications = 3L),
    experiment_2b = list(
      n_values = c(30L, 300L),
      replications = setNames(c(3L, 3L), c("30", "300")),
      N = 12L
    ),
    experiment_3a = list(n_values = c(30L, 60L), replications = c(`30` = 3L, `60` = 3L)),
    experiment_3b = list(n_values = c(30L, 60L), replications = c(`30` = 3L, `60` = 3L)),
    experiment_3c = list(n = 30L, T = 10L, N = 20L, m_values = c(50L, 80L), replications = 3L),
    experiment_3d = list(
      n_values = c(30L, 60L), replications = c(`30` = 6L, `60` = 6L),
      delta_grid = c(0, 1)
    ),
    robustness_k0 = list(
      n_values = c(30L, 60L),
      replications = setNames(c(3L, 3L), c("30", "60")),
      N = 12L, K0_values = c(3L, 5L, 8L)
    ),
    robustness_universality = list(n = 20L, p = 10L, N = 12L, replications = 3L)
  ))
}

run_pipeline <- function(mode = PIPELINE_MODE, config = PIPELINE_CONFIG) {
  mode <- match.arg(mode, pipeline_modes())
  root <- if (!is.null(config$root)) config$root else PROJECT_ROOT
  selected_runs <- if (!is.null(config$result_run_ids)) config$result_run_ids else list()
  if (mode == "help") return(invisible(print_pipeline_help()))
  if (mode == "simulate_1a_panel_phase") return(run_panel_phase_experiment(experiment_config(default_panel_phase_config(), config, "experiment_1a")))
  if (mode == "simulate_1b_phase_margin") return(run_phase_margin_experiment(experiment_config(default_phase_margin_config(), config, "experiment_1b")))
  if (mode == "simulate_2_functional_inference") return(run_functional_inference_experiment(experiment_config(default_functional_inference_config(), config, "experiment_2")))
  if (mode == "simulate_2b_functional_asymptotics") return(run_functional_asymptotic_experiment(experiment_config(default_functional_asymptotic_config(), config, "experiment_2b")))
  if (mode == "simulate_3a_simple_direction") return(run_simple_direction_experiment(experiment_config(default_large_domain_config(), config, "experiment_3a")))
  if (mode == "simulate_3b_repeated_eigenspace") return(run_repeated_eigenspace_experiment(experiment_config(default_large_domain_config(), config, "experiment_3b")))
  if (mode == "simulate_3c_grid_saturation") return(run_grid_saturation_experiment(experiment_config(default_grid_saturation_config(), config, "experiment_3c")))
  if (mode == "simulate_3d_eigengap_inference") return(run_eigengap_experiment(experiment_config(default_eigengap_config(), config, "experiment_3d")))
  if (mode == "simulate_robustness_k0") return(run_k0_robustness_experiment(experiment_config(default_k0_robustness_config(), config, "robustness_k0")))
  if (mode == "simulate_robustness_universality") return(run_universality_stress_experiment(experiment_config(default_universality_config(), config, "robustness_universality")))
  if (mode == "simulate_all") return(run_all_simulations(config))
  if (mode == "summarize_all") return(generate_all_summaries(root, selected_runs))
  if (mode == "plot_figure_2") return(plot_figure_2(root, selected_runs))
  if (mode == "plot_figure_3") return(plot_figure_3(root, selected_runs))
  if (mode == "plot_figure_4") return(plot_figure_4(root, selected_runs))
  if (mode == "plot_robustness") return(plot_robustness_results(root, selected_runs))
  if (mode == "plot_all") return(plot_all_saved_results(root, selected_runs))
  if (mode == "all") {
    run_all_simulations(config)
    return(plot_all_saved_results(root, selected_runs))
  }
  if (mode == "smoke") {
    run_eigengap_smoke_checks()
    smoke <- smoke_config(config)
    on.exit(unlink(smoke$root, recursive = TRUE, force = TRUE), add = TRUE)
    run_all_simulations(smoke)
    return(plot_all_saved_results(smoke$root))
  }
}

parse_cli_overrides <- function(args) {
  output <- list()
  for (arg in args) {
    separator <- regexpr("=", arg, fixed = TRUE)[1]
    if (separator < 2L) stop("CLI arguments must use key=value syntax: ", arg)
    key <- substr(arg, 1L, separator - 1L)
    value <- substr(arg, separator + 1L, nchar(arg))
    if (key %in% c("mode", "experiment", "dataset", "action")) {
      output[[key]] <- gsub("^['\"]|['\"]$", "", value)
    } else {
      output[[key]] <- tryCatch(eval(parse(text = value), envir = baseenv()), error = function(error) value)
    }
  }
  aliases <- c(R = "replications", cores = "ncores", seed = "master_seed", m = "m_values", Delta = "Delta_grid", K0s = "K0_values")
  for (alias in intersect(names(aliases), names(output))) {
    output[[aliases[[alias]]]] <- output[[alias]]
    output[[alias]] <- NULL
  }
  output
}

if (sys.nframe() == 0L) {
  cli <- parse_cli_overrides(commandArgs(trailingOnly = TRUE))
  if (!is.null(cli$mode)) {
    mode <- cli$mode; cli$mode <- NULL
    run_pipeline(mode, merge_config(PIPELINE_CONFIG, cli))
  } else {
    selection <- if (!is.null(cli$experiment)) cli$experiment else PIPELINE_EXPERIMENT
    action <- if (!is.null(cli$action)) cli$action else PIPELINE_ACTION
    parameter_names <- setdiff(names(cli), c("experiment", "action"))
    cli$experiment <- NULL; cli$action <- NULL
    if (identical(tolower(selection), "help")) {
      print_pipeline_help()
    } else {
      selected_ids <- expand_experiment_selection(selection)
      cli_n <- cli[["n", exact = TRUE]]
      cli_n_values <- cli[["n_values", exact = TRUE]]
      if (any(selected_ids %in% c("2b", "3a", "3b", "3d", "k0")) &&
          !is.null(cli_n) && is.null(cli_n_values)) cli$n_values <- cli_n
      merged_config <- merge_config(PIPELINE_CONFIG, cli)
      catalog <- experiment_catalog()
      has_local_overrides <- any(vapply(selected_ids, function(id) length(merged_config[[catalog[[id]]$config_key]]) > 0L, logical(1)))
      scientific_parameters <- setdiff(parameter_names, c("root", "ncores"))
      run_selected_experiments(selection, action, merged_config,
                               exact_saved_config = length(scientific_parameters) > 0L || has_local_overrides)
    }
  }
}
