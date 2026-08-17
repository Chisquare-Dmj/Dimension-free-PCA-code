# Robustness experiments: K0 sensitivity and universality stress testing.

default_k0_robustness_config <- function() {
  list(
    n_values = c(150L, 300L, 600L),
    replications = setNames(c(2000L, 1500L, 1000L), c("150", "300", "600")),
    p_rule = FUNCTIONAL_ASYMPTOTIC_P_RULE,
    N = 200L, decay = 1.1, M = 3L,
    K0_values = c(3L, 5L, 8L),
    scenario = "F3_Haar", distributions = "Uniform",
    primary_distribution = "Uniform", spike_indices = 2L,
    population_seed = 201L, ncores = 1L, master_seed = MASTER_SEED,
    confidence_level = 0.95, certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

prepare_k0_robustness_config <- function(config) {
  requested_distribution <- config[["distribution"]]
  if (!is.null(requested_distribution)) config$distributions <- requested_distribution
  config$n_values <- sort(unique(as.integer(config$n_values)))
  config$K0_values <- sort(unique(as.integer(config$K0_values)))
  config$p_rule <- FUNCTIONAL_ASYMPTOTIC_P_RULE
  config$distribution_grid <- paste(config$distributions, collapse = "-")
  config
}

validate_k0_robustness_config <- function(config) {
  validate_positive_integer(config$N, "N")
  validate_positive_integer(config$M, "M")
  if (!length(config$n_values)) stop("n_values must not be empty.")
  invisible(lapply(config$n_values, validate_positive_integer, name = "n_values entry"))
  if (!length(config$K0_values)) stop("K0_values must not be empty.")
  invisible(lapply(config$K0_values, validate_positive_integer, name = "K0_values entry"))
  if (!config$M %in% config$K0_values) stop("K0_values must include the K0=M baseline.")
  if (any(config$K0_values < config$M | config$K0_values >= min(config$n_values))) {
    stop("Every K0_values entry must satisfy M <= K0 < min(n_values).")
  }
  if (is.null(names(config$replications)) ||
      !all(as.character(config$n_values) %in% names(config$replications))) {
    stop("replications must be named by every n value.")
  }
  invisible(lapply(
    config$replications[as.character(config$n_values)],
    validate_positive_integer, name = "replications entry"
  ))
  if (!length(config$spike_indices) ||
      any(config$spike_indices < 1L | config$spike_indices > config$M)) {
    stop("spike_indices must lie between 1 and M.")
  }
  if (!length(config$distributions)) stop("distributions must not be empty.")
  if (!config$primary_distribution %in% config$distributions) {
    stop("primary_distribution must be included in distributions.")
  }
  if (config$confidence_level <= 0 || config$confidence_level >= 1) {
    stop("confidence_level must lie strictly between 0 and 1.")
  }
  invisible(config)
}

run_k0_robustness_experiment <- function(config = list()) {
  cfg <- prepare_k0_robustness_config(
    merge_config(default_k0_robustness_config(), config)
  )
  validate_k0_robustness_config(cfg)
  experiment <- "robustness_k0"
  run_id <- make_run_id(
    experiment, cfg,
    c(
      "n_values", "replications", "p_rule", "N", "M", "K0_values",
      "spike_indices", "scenario", "distribution_grid"
    )
  )
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list()
  all_truth <- list()
  result_index <- 1L
  for (n_index in seq_along(cfg$n_values)) {
    n <- as.integer(cfg$n_values[n_index])
    p <- functional_asymptotic_p(n)
    replications <- as.integer(cfg$replications[[as.character(n)]])
    population <- build_functional_population(
      cfg$scenario, p, cfg$N, cfg$decay, cfg$population_seed
    )
    truth <- population_spike_truth(population$values, cfg$M, n, cfg$spike_indices)
    truth$experiment <- experiment
    truth$run_id <- run_id
    truth$scenario <- cfg$scenario
    truth$n <- n
    truth$p <- p
    truth$N <- cfg$N
    truth$M <- cfg$M
    truth$K0_values <- compact_config_value(cfg$K0_values)
    truth$distribution_grid <- cfg$distribution_grid
    all_truth[[n_index]] <- truth

    for (distribution_id in seq_along(cfg$distributions)) {
      distribution <- cfg$distributions[distribution_id]
      worker <- function(replication) {
        seed_id <- 70L + 10L * n_index + distribution_id
        set.seed(replication_seed(seed_id, replication, cfg$master_seed))
        X <- sample_functional_coefficients(n, population, distribution)
        fit <- gram_pca(X, vectors = FALSE, center = FALSE)
        rows <- list()
        row_index <- 1L
        for (K0 in cfg$K0_values) {
          for (j in cfg$spike_indices) {
            estimate <- dimension_free_inference(
              fit$values, j, K0, n,
              cfg$confidence_level, cfg$certification_threshold
            )
            row <- attach_truth_and_coverage(
              estimate, truth[truth$spike_index == j, ]
            )
            row$experiment <- experiment
            row$scenario <- cfg$scenario
            row$score_distribution <- distribution
            row$n <- n
            row$p <- p
            row$T <- NA_integer_
            row$N <- cfg$N
            row$m <- NA_integer_
            row$replication <- replication
            rows[[row_index]] <- row
            row_index <- row_index + 1L
          }
        }
        result <- bind_rows_base(rows)
        result$overdeflation_excess <- result$K0 - cfg$M
        result$overdeflation_error_alpha <- NA_real_
        result$overdeflation_error_Delta <- NA_real_
        result$overdeflation_error_r2 <- NA_real_
        for (j in cfg$spike_indices) {
          target <- result$spike_index == j
          reference <- result[target & result$K0 == cfg$M, ]
          result$overdeflation_error_alpha[target] <- result$hat_alpha[target] - reference$hat_alpha
          result$overdeflation_error_Delta[target] <- result$hat_Delta[target] - reference$hat_Delta
          result$overdeflation_error_r2[target] <- result$hat_r2[target] - reference$hat_r2
        }
        result
      }
      result <- bind_rows_base(parallel_map(seq_len(replications), worker, cfg$ncores))
      all_results[[result_index]] <- add_run_metadata(
        result, cfg, run_id, replications, cfg$population_seed
      )
      progress_message(
        paste("K0 robustness", distribution, "n=", n, "p=", p),
        replications, replications
      )
      result_index <- result_index + 1L
    }
  }
  out <- standardize_replicate_output(bind_rows_base(all_results))
  truth_df <- bind_rows_base(all_truth)
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(out, replicate_path)
  write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(out)
}

default_universality_config <- function() {
  list(
    n = 300L, p = 100L, N = 200L, decay = 1.1, M = 3L, K0 = 5L,
    replications = 2000L, scenario = "F1_diagonal",
    distributions = c("Gaussian", "t12", "Uniform"), spike_indices = 1:2,
    ncores = 1L, master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

run_universality_stress_experiment <- function(config = list()) {
  cfg <- merge_config(default_universality_config(), config)
  validate_inference_config(cfg, c("n", "p", "N"))
  if (!length(cfg$distributions)) stop("distributions must not be empty.")
  if (!length(cfg$spike_indices) || any(cfg$spike_indices < 1L | cfg$spike_indices > cfg$M)) stop("spike_indices must lie between 1 and M.")
  experiment <- "robustness_universality"
  cfg$distribution_grid <- paste(cfg$distributions, collapse = "-")
  run_id <- make_run_id(experiment, cfg, c("n", "p", "N", "replications", "M", "K0", "scenario", "distribution_grid"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  population <- build_functional_population(cfg$scenario, cfg$p, cfg$N, cfg$decay, seed = 0L)
  truth <- population_spike_truth(population$values, cfg$M, cfg$n, cfg$spike_indices)
  z <- qnorm(1 - (1 - cfg$confidence_level) / 2)
  all_results <- list()

  for (distribution_id in seq_along(cfg$distributions)) {
    distribution <- cfg$distributions[distribution_id]
    kappa <- score_fourth_cumulant(distribution)
    worker <- function(replication) {
      set.seed(replication_seed(8L, replication, cfg$master_seed))
      X <- sample_functional_coefficients(cfg$n, population, distribution)
      fit <- gram_pca(X, vectors = FALSE, center = FALSE)
      rows <- lapply(cfg$spike_indices, function(j) {
        estimate <- dimension_free_inference(fit$values, j, cfg$K0, cfg$n, cfg$confidence_level, cfg$certification_threshold)
        row_truth <- truth[truth$spike_index == j, ]
        vpsi_full <- kappa * row_truth$true_alpha^2 * row_truth$true_Delta^2 +
          2 * row_truth$true_alpha^2 * row_truth$true_Delta
        full_multiplier <- sqrt(1 + kappa * row_truth$true_Delta / 2)
        oracle_se_alpha <- estimate$se_alpha * full_multiplier
        oracle_se_Delta <- estimate$se_Delta * full_multiplier
        oracle_se_r2 <- estimate$se_r2 * full_multiplier
        row <- attach_truth_and_coverage(estimate, row_truth)
        row$experiment <- experiment; row$scenario <- cfg$scenario
        row$score_distribution <- distribution
        row$n <- cfg$n; row$p <- cfg$p; row$T <- NA_integer_; row$N <- cfg$N; row$m <- NA_integer_
        row$replication <- replication
        row$fourth_cumulant <- kappa; row$vpsi_oracle_full <- vpsi_full
        row$oracle_se_alpha <- oracle_se_alpha; row$oracle_se_Delta <- oracle_se_Delta; row$oracle_se_r2 <- oracle_se_r2
        row$oracle_cover_alpha <- as.integer(abs(estimate$hat_alpha - row_truth$true_alpha) <= z * oracle_se_alpha)
        row$oracle_cover_Delta <- as.integer(abs(estimate$hat_Delta - row_truth$true_Delta) <= z * oracle_se_Delta)
        row$oracle_cover_r2 <- as.integer(abs(estimate$hat_r2 - row_truth$true_r2) <= z * oracle_se_r2)
        row
      })
      bind_rows_base(rows)
    }
    all_results[[distribution]] <- bind_rows_base(parallel_map(seq_len(cfg$replications), worker, cfg$ncores))
    progress_message(paste("Universality", distribution), cfg$replications, cfg$replications)
  }
  out <- add_run_metadata(bind_rows_base(all_results), cfg, run_id, cfg$replications)
  out <- standardize_replicate_output(out)
  truth$experiment <- experiment; truth$run_id <- run_id; truth$M <- cfg$M; truth$K0 <- cfg$K0
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(out, replicate_path); write_csv_atomic(truth, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(out)
}
