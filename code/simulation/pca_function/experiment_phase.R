# Experiments 1A and 1B: panel phase transition and phase-margin conditioning.

default_panel_phase_config <- function() {
  list(
    n = 300L, p = 20L, T = 30L, replications = 500L, M = 3L, K0 = 5L,
    scenarios = c("P1_independent", "P2_block_AR1"),
    distribution = "Gaussian", ncores = 1L,
    p2_leading_multipliers = c(5, 4), target_Delta3 = -0.10,
    population_seeds = c(P1_independent = 101L, P2_block_AR1 = 102L),
    master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

run_panel_phase_experiment <- function(config = list()) {
  cfg <- merge_config(default_panel_phase_config(), config)
  validate_inference_config(cfg, c("n", "p", "T"))
  if (!length(cfg$scenarios)) stop("scenarios must not be empty.")
  if (!all(cfg$scenarios %in% names(cfg$population_seeds))) stop("population_seeds must name every requested panel scenario.")
  if ("P2_block_AR1" %in% cfg$scenarios && cfg$T %% 10L != 0L) stop("Experiment 1A P2 requires T divisible by 10.")
  if (length(cfg$p2_leading_multipliers) != 2L || any(cfg$p2_leading_multipliers <= 1)) stop("p2_leading_multipliers must contain two values greater than one.")
  if (length(cfg$target_Delta3) != 1L || cfg$target_Delta3 >= 0) stop("target_Delta3 must be negative for a subcritical third component.")
  experiment <- "experiment_1a_panel_phase"
  run_id <- make_run_id(experiment, cfg, c("n", "p", "T", "replications", "M", "K0", "target_Delta3", "distribution"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_replications <- list()
  all_bulk <- list()
  all_truth <- list()

  for (scenario_id in seq_along(cfg$scenarios)) {
    scenario <- cfg$scenarios[scenario_id]
    population <- build_panel_population(
      scenario, cfg$p, cfg$T,
      seed = unname(cfg$population_seeds[[scenario]]),
      sample_size = cfg$n,
      p2_leading_multipliers = cfg$p2_leading_multipliers,
      target_Delta3 = cfg$target_Delta3
    )
    truth <- population_spike_truth(population$values, M = cfg$M, n = cfg$n, spike_indices = 1:3)
    truth$experiment <- experiment
    truth$run_id <- run_id
    truth$scenario <- scenario
    truth$score_distribution <- cfg$distribution
    truth$n <- cfg$n; truth$p <- cfg$p; truth$T <- cfg$T
    truth$M <- cfg$M; truth$K0 <- cfg$K0
    truth$population_seed <- unname(cfg$population_seeds[[scenario]])
    truth$calibrated_multiplier3 <- population$calibrated_multiplier3
    truth$population_phase_status <- ifelse(truth$true_Delta > 0, "supercritical", "subcritical")
    truth$true_r2[truth$true_Delta <= 0] <- NA_real_
    all_truth[[scenario]] <- truth

    worker <- function(replication) {
      set.seed(replication_seed(1L, replication, cfg$master_seed))
      X <- sample_panel_observations(cfg$n, population, cfg$distribution)
      fit <- gram_pca(X, vectors = FALSE, center = FALSE)
      estimates <- lapply(1:3, function(j) {
        if (j <= 2L) {
          estimate <- dimension_free_inference(
            fit$values, j, cfg$K0, cfg$n,
            cfg$confidence_level, cfg$certification_threshold
          )
          row <- attach_truth_and_coverage(estimate, truth[truth$spike_index == j, ])
          row$regular_inference_valid <- TRUE
        } else {
          row <- cbind(
            data.frame(spike_index = j, K0 = cfg$K0, hat_lambda = fit$values[j]),
            truth[truth$spike_index == j, setdiff(names(truth), "spike_index"), drop = FALSE]
          )
          row$regular_inference_valid <- FALSE
          row$cover_alpha <- NA_integer_; row$cover_Delta <- NA_integer_; row$cover_r2 <- NA_integer_
        }
        row$experiment <- experiment
        row$scenario <- scenario
        row$score_distribution <- cfg$distribution
        row$n <- cfg$n; row$p <- cfg$p; row$T <- cfg$T
        row$N <- NA_integer_; row$m <- NA_integer_
        row$replication <- replication
        row$sample_gap_to_lambda4 <- fit$values[j] - fit$values[4]
        row
      })
      bulk <- data.frame(
        experiment = experiment, run_id = run_id,
        scenario = scenario,
        score_distribution = cfg$distribution,
        n = cfg$n, p = cfg$p, T = cfg$T,
        replication = replication,
        bulk_rank = seq.int(4L, length(fit$values)),
        bulk_eigenvalue = fit$values[4:length(fit$values)]
      )
      list(spikes = bind_rows_base(estimates), bulk = bulk)
    }
    results <- parallel_map(seq_len(cfg$replications), worker, cfg$ncores)
    all_replications[[scenario]] <- add_run_metadata(
      bind_rows_base(lapply(results, `[[`, "spikes")), cfg, run_id,
      cfg$replications, unname(cfg$population_seeds[[scenario]])
    )
    all_bulk[[scenario]] <- add_run_metadata(
      bind_rows_base(lapply(results, `[[`, "bulk")), cfg, run_id,
      cfg$replications, unname(cfg$population_seeds[[scenario]])
    )
    progress_message(paste("Experiment 1A", scenario), cfg$replications, cfg$replications)
  }

  replication_df <- standardize_replicate_output(bind_rows_base(all_replications))
  bulk_df <- bind_rows_base(all_bulk)
  truth_df <- bind_rows_base(all_truth)
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  bulk_path <- run_artifact_path("replicate", run_id, "bulk-spectrum", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(replication_df, replicate_path)
  write_csv_atomic(bulk_df, bulk_path)
  write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, bulk = bulk_path, truth = truth_path, config = config_path), cfg$root)
  invisible(replication_df)
}

default_phase_margin_config <- function() {
  list(
    n = 300L, p = 300L, M = 1L, K0 = 1L,
    Delta_grid_main = seq(0.40, 0.90, by = 0.10),
    Delta_grid_boundary = seq(0.10, 0.30, by = 0.10),
    # Optional backward-compatible override for a custom combined grid.
    Delta_grid = NULL,
    replications = 2000L, distribution = "Gaussian", ncores = 1L,
    master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

prepare_phase_margin_config <- function(config) {
  if (!is.null(config$Delta_grid)) {
    requested_grid <- sort(unique(as.numeric(config$Delta_grid)))
    config$Delta_grid_main <- requested_grid[requested_grid >= 0.40]
    config$Delta_grid_boundary <- requested_grid[requested_grid < 0.40]
  }
  config$Delta_grid_main <- sort(unique(as.numeric(config$Delta_grid_main)))
  config$Delta_grid_boundary <- sort(unique(as.numeric(config$Delta_grid_boundary)))
  config$Delta_grid <- sort(unique(c(config$Delta_grid_boundary, config$Delta_grid_main)))
  config
}

run_phase_margin_experiment <- function(config = list()) {
  cfg <- prepare_phase_margin_config(merge_config(default_phase_margin_config(), config))
  validate_inference_config(cfg, c("n", "p"))
  if (cfg$M != 1L || cfg$K0 != 1L) stop("Experiment 1B baseline requires M=K0=1; over-deflation is studied separately.")
  if (!length(cfg$Delta_grid)) stop("Delta_grid must not be empty.")
  if (any(cfg$Delta_grid <= 0 | cfg$Delta_grid >= 1)) stop("Delta_grid values must lie strictly between 0 and 1.")
  experiment <- "experiment_1b_phase_margin"
  run_id <- make_run_id(
    experiment, cfg,
    c("n", "p", "replications", "M", "K0", "Delta_grid_main", "Delta_grid_boundary", "distribution")
  )
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list()
  all_truth <- list()

  for (scenario_id in seq_along(cfg$Delta_grid)) {
    Delta <- cfg$Delta_grid[scenario_id]
    population <- build_vector_calibration_population(cfg$n, cfg$p, Delta)
    truth <- population_spike_truth(population$values, cfg$M, cfg$n, 1L)
    truth$target_Delta <- Delta
    truth$phase_region <- if (Delta %in% cfg$Delta_grid_main) "main" else "boundary"
    all_truth[[scenario_id]] <- truth

    worker <- function(replication) {
      set.seed(replication_seed(2L, replication, cfg$master_seed))
      X <- sample_vector_observations(cfg$n, population, cfg$distribution)
      fit <- gram_pca(X, vectors = FALSE, center = FALSE)
      estimate <- dimension_free_inference(
        fit$values, 1L, cfg$K0, cfg$n,
        cfg$confidence_level, cfg$certification_threshold
      )
      row <- attach_truth_and_coverage(estimate, truth)
      row$experiment <- experiment
      row$scenario <- sprintf("Delta_%0.2f", Delta)
      row$score_distribution <- cfg$distribution
      row$n <- cfg$n; row$p <- cfg$p; row$T <- NA_integer_; row$N <- NA_integer_; row$m <- NA_integer_
      row$replication <- replication
      row$target_Delta <- Delta
      row$phase_region <- truth$phase_region
      row$relative_bias_naive <- fit$values[1] / truth$true_alpha - 1
      row$relative_bias_corrected <- estimate$hat_alpha / truth$true_alpha - 1
      row$ci_length_alpha <- estimate$ci_alpha_upper - estimate$ci_alpha_lower
      row$ci_length_Delta <- estimate$ci_Delta_upper - estimate$ci_Delta_lower
      row$ci_length_r2 <- estimate$ci_r2_upper - estimate$ci_r2_lower
      row
    }
    all_results[[scenario_id]] <- bind_rows_base(parallel_map(seq_len(cfg$replications), worker, cfg$ncores))
    progress_message(sprintf("Experiment 1B Delta=%.2f", Delta), cfg$replications, cfg$replications)
  }

  replication_df <- add_run_metadata(bind_rows_base(all_results), cfg, run_id, cfg$replications)
  replication_df <- standardize_replicate_output(replication_df)
  truth_df <- bind_rows_base(all_truth)
  truth_df$experiment <- experiment; truth_df$run_id <- run_id; truth_df$M <- cfg$M; truth_df$K0 <- cfg$K0
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(replication_df, replicate_path)
  write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(replication_df)
}
