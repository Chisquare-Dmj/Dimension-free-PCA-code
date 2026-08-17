# Experiment 2: complete inference calibration for high-dimensional functional data.

default_functional_inference_config <- function() {
  list(
    n = 300L, p = 100L, N = 200L, decay = 1.1, M = 3L, K0 = 5L,
    scenarios = c("F1_diagonal", "F2_block", "F3_Haar"),
    distributions = c("Gaussian", "t12", "Uniform"),
    core_combinations = c(
      "F1_diagonal|Gaussian", "F2_block|Gaussian", "F3_Haar|Gaussian",
      "F3_Haar|t12", "F3_Haar|Uniform"
    ),
    core_replications = 2000L, other_replications = 1000L,
    spike_indices = 1:2, ncores = 1L, population_seed = 201L,
    master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

run_functional_inference_experiment <- function(config = list()) {
  cfg <- merge_config(default_functional_inference_config(), config)
  validate_inference_config(cfg, c("n", "p", "N"))
  validate_positive_integer(cfg$core_replications, "core_replications")
  validate_positive_integer(cfg$other_replications, "other_replications")
  if ("F2_block" %in% cfg$scenarios && cfg$p %% 10L != 0L) stop("Experiment 2 F2 requires p divisible by 10.")
  if (!length(cfg$scenarios) || !length(cfg$distributions)) stop("scenarios and distributions must not be empty.")
  if (!length(cfg$spike_indices) || any(cfg$spike_indices < 1L | cfg$spike_indices > cfg$M)) stop("spike_indices must lie between 1 and M.")
  experiment <- "experiment_2_functional_inference"
  cfg$design_grid <- paste0(length(cfg$scenarios), "Sx", length(cfg$distributions), "D")
  run_id <- make_run_id(experiment, cfg, c("n", "p", "N", "decay", "M", "K0", "core_replications", "other_replications", "design_grid"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list()
  all_truth <- list()
  combination_id <- 0L

  for (scenario in cfg$scenarios) {
    population <- build_functional_population(scenario, cfg$p, cfg$N, cfg$decay, cfg$population_seed)
    truth <- population_spike_truth(population$values, cfg$M, cfg$n, cfg$spike_indices)
    derivatives <- lapply(cfg$spike_indices, function(j) {
      values <- population_derivatives(population$values, cfg$M, cfg$n, j)
      names(values) <- paste0("true_", names(values))
      as.data.frame(as.list(values))
    })
    derivatives <- bind_rows_base(derivatives)
    derivatives$spike_index <- cfg$spike_indices
    truth <- merge(truth, derivatives, by = "spike_index", suffixes = c("", "_population"), sort = FALSE)

    for (distribution in cfg$distributions) {
      combination_id <- combination_id + 1L
      combination <- paste(scenario, distribution, sep = "|")
      replications <- if (combination %in% cfg$core_combinations) cfg$core_replications else cfg$other_replications
      scenario_truth <- truth
      scenario_truth$experiment <- experiment
      scenario_truth$run_id <- run_id
      scenario_truth$scenario <- scenario
      scenario_truth$score_distribution <- distribution
      scenario_truth$n <- cfg$n; scenario_truth$p <- cfg$p; scenario_truth$N <- cfg$N
      all_truth[[combination]] <- scenario_truth
      target_directions <- lapply(cfg$spike_indices, function(j) functional_population_direction(population, j))

      worker <- function(replication) {
        set.seed(replication_seed(3L, replication, cfg$master_seed))
        X <- sample_functional_coefficients(cfg$n, population, distribution)
        fit <- gram_pca(X, vectors = TRUE, center = FALSE)
        rows <- lapply(seq_along(cfg$spike_indices), function(index) {
          j <- cfg$spike_indices[index]
          estimate <- dimension_free_inference(
            fit$values, j, cfg$K0, cfg$n,
            cfg$confidence_level, cfg$certification_threshold
          )
          row_truth <- scenario_truth[scenario_truth$spike_index == j, ]
          row <- attach_truth_and_coverage(estimate, row_truth)
          overlap <- feature_projection(fit, j, target_directions[[index]])
          row$experiment <- experiment
          row$scenario <- scenario
          row$score_distribution <- distribution
          row$n <- cfg$n; row$p <- cfg$p; row$T <- NA_integer_; row$N <- cfg$N; row$m <- NA_integer_
          row$replication <- replication
          row$actual_signal_overlap <- overlap^2
          row$theta_true <- NA_real_; row$theta_raw <- NA_real_; row$theta_corrected <- NA_real_
          row$Q_true <- NA_real_; row$Q_raw <- NA_real_; row$Q_corrected <- NA_real_
          row$individual_alignment_repeated <- NA_real_
          row
        })
        bind_rows_base(rows)
      }
      all_results[[combination]] <- add_run_metadata(
        bind_rows_base(parallel_map(seq_len(replications), worker, cfg$ncores)),
        cfg, run_id, replications, cfg$population_seed
      )
      progress_message(paste("Experiment 2", scenario, distribution), replications, replications)
    }
  }

  replication_df <- standardize_replicate_output(bind_rows_base(all_results))
  truth_df <- bind_rows_base(all_truth)
  truth_df$run_id <- run_id; truth_df$M <- cfg$M; truth_df$K0 <- cfg$K0
  truth_df$decay <- cfg$decay; truth_df$population_seed <- cfg$population_seed
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(replication_df, replicate_path)
  write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(replication_df)
}
