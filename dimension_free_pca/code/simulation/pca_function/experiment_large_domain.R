# Experiments 3A--3C: direction, repeated eigenspace, and dense-grid approximation.

default_large_domain_config <- function() {
  list(
    n_values = c(150L, 300L, 600L),
    replications = c(`150` = 2000L, `300` = 1500L, `600` = 1000L),
    decay = 1.1, M = 3L, K0 = 5L, distribution = "Gaussian",
    ncores = 1L, master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

validate_large_domain_config <- function(config, label) {
  validate_positive_integer(config$M, "M")
  validate_positive_integer(config$K0, "K0")
  if (config$K0 < config$M) stop("K0 must be at least M.")
  if (!length(config$n_values)) stop(label, " requires at least one n value.")
  invisible(lapply(config$n_values, validate_positive_integer, name = "n_values entry"))
  if (config$K0 >= min(config$n_values)) stop("K0 must be smaller than every n value.")
  if (any(config$n_values %% 3L != 0L)) stop(label, " requires every n value divisible by 3 so T=n/3 is integral.")
  if (is.null(names(config$replications)) || !all(as.character(config$n_values) %in% names(config$replications))) {
    stop("replications must be named by every n value.")
  }
  invisible(lapply(config$replications[as.character(config$n_values)], validate_positive_integer, name = "replications entry"))
  if (!is.numeric(config$decay) || length(config$decay) != 1L || config$decay <= 0) stop("decay must be positive.")
  if (config$confidence_level <= 0 || config$confidence_level >= 1) stop("confidence_level must lie strictly between 0 and 1.")
  invisible(config)
}

run_simple_direction_experiment <- function(config = list()) {
  cfg <- merge_config(default_large_domain_config(), config)
  validate_large_domain_config(cfg, "Experiment 3A")
  cfg$domain_rule <- "T-n-div-3"; cfg$basis_rule <- "N-2n"
  experiment <- "experiment_3a_simple_direction"
  run_id <- make_run_id(experiment, cfg, c("n_values", "replications", "domain_rule", "basis_rule", "decay", "M", "K0", "distribution"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list(); all_truth <- list()

  for (scenario_id in seq_along(cfg$n_values)) {
    n <- as.integer(cfg$n_values[scenario_id]); T <- as.integer(n / 3); N <- as.integer(2 * n)
    R <- as.integer(cfg$replications[[as.character(n)]])
    values <- large_domain_spectrum(n, repeated = FALSE, N = N, T = T, decay = cfg$decay)
    truth <- population_spike_truth(values, cfg$M, n, 1L)
    f <- large_domain_target(N)
    all_truth[[scenario_id]] <- transform(truth, experiment = "experiment_3a_simple_direction", n = n, T = T, N = N, theta_true = 0.6)

    worker <- function(replication) {
      set.seed(replication_seed(4L, replication, cfg$master_seed))
      X <- sample_large_domain_coefficients(n, values, cfg$distribution)
      fit <- gram_pca(X, vectors = TRUE, center = FALSE)
      estimate <- dimension_free_inference(fit$values, 1L, cfg$K0, n, cfg$confidence_level, cfg$certification_threshold)
      row <- attach_truth_and_coverage(estimate, truth)
      anchor <- feature_projection(fit, 1L, c(1, rep(0, N - 1L)))
      raw <- orient_projection(feature_projection(fit, 1L, f), anchor)
      corrected <- raw / sqrt(estimate$hat_r2)
      row$experiment <- experiment
      row$scenario <- "simple_spike"; row$score_distribution <- cfg$distribution
      row$n <- n; row$p <- NA_integer_; row$T <- T; row$N <- N; row$m <- NA_integer_
      row$replication <- replication
      row$actual_signal_overlap <- anchor^2
      row$theta_true <- 0.6; row$theta_raw <- raw; row$theta_corrected <- corrected
      row$theta_raw_squared_error <- (raw - 0.6)^2
      row$theta_corrected_squared_error <- (corrected - 0.6)^2
      row$Q_true <- NA_real_; row$Q_raw <- NA_real_; row$Q_corrected <- NA_real_
      row$individual_alignment_repeated <- NA_real_
      row
    }
    all_results[[scenario_id]] <- add_run_metadata(
      bind_rows_base(parallel_map(seq_len(R), worker, cfg$ncores)), cfg, run_id, R
    )
    progress_message(paste("Experiment 3A n=", n), R, R)
  }
  out <- standardize_replicate_output(bind_rows_base(all_results))
  truth_df <- bind_rows_base(all_truth); truth_df$run_id <- run_id; truth_df$M <- cfg$M; truth_df$K0 <- cfg$K0; truth_df$decay <- cfg$decay
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(out, replicate_path); write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(out)
}

run_repeated_eigenspace_experiment <- function(config = list()) {
  cfg <- merge_config(default_large_domain_config(), config)
  validate_large_domain_config(cfg, "Experiment 3B")
  cfg$domain_rule <- "T-n-div-3"; cfg$basis_rule <- "N-2n"
  experiment <- "experiment_3b_repeated_eigenspace"
  run_id <- make_run_id(experiment, cfg, c("n_values", "replications", "domain_rule", "basis_rule", "decay", "M", "K0", "distribution"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list(); all_truth <- list()

  for (scenario_id in seq_along(cfg$n_values)) {
    n <- as.integer(cfg$n_values[scenario_id]); T <- as.integer(n / 3); N <- as.integer(2 * n)
    R <- as.integer(cfg$replications[[as.character(n)]])
    values <- large_domain_spectrum(n, repeated = TRUE, N = N, T = T, decay = cfg$decay)
    truth <- population_spike_truth(values, cfg$M, n, 1L)
    f <- large_domain_target(N)
    all_truth[[scenario_id]] <- transform(truth, experiment = "experiment_3b_repeated_eigenspace", n = n, T = T, N = N, Q_true = 0.45)

    worker <- function(replication) {
      set.seed(replication_seed(5L, replication, cfg$master_seed))
      X <- sample_large_domain_coefficients(n, values, cfg$distribution)
      fit <- gram_pca(X, vectors = TRUE, center = FALSE)
      estimates <- lapply(1:2, function(j) dimension_free_inference(fit$values, j, cfg$K0, n, cfg$confidence_level, cfg$certification_threshold))
      cluster_r2 <- mean(vapply(estimates, function(x) x$hat_r2, numeric(1)))
      projections_f <- vapply(1:2, function(j) feature_projection(fit, j, f), numeric(1))
      raw_Q <- sum(projections_f^2)
      alignment <- feature_projection(fit, 1L, c(1, rep(0, N - 1L)))^2
      data.frame(
        experiment = experiment, scenario = "repeated_spike",
        score_distribution = cfg$distribution, n = n, p = NA_integer_, T = T, N = N, m = NA_integer_,
        replication = replication, spike_index = 1L, K0 = cfg$K0,
        true_alpha = truth$true_alpha, true_psi = truth$true_psi,
        true_Delta = truth$true_Delta, true_r2 = truth$true_r2,
        hat_lambda = mean(vapply(estimates, function(x) x$hat_lambda, numeric(1))),
        hat_alpha = mean(vapply(estimates, function(x) x$hat_alpha, numeric(1))),
        hat_Delta = mean(vapply(estimates, function(x) x$hat_Delta, numeric(1))),
        hat_r2 = cluster_r2, cluster_hat_r2 = cluster_r2,
        se_alpha = NA_real_, se_Delta = NA_real_, se_r2 = NA_real_,
        actual_signal_overlap = NA_real_, theta_true = NA_real_, theta_raw = NA_real_, theta_corrected = NA_real_,
        Q_true = 0.45, Q_raw = raw_Q, Q_corrected = raw_Q / cluster_r2,
        Q_raw_squared_error = (raw_Q - 0.45)^2,
        Q_corrected_squared_error = (raw_Q / cluster_r2 - 0.45)^2,
        individual_alignment_repeated = alignment
      )
    }
    all_results[[scenario_id]] <- add_run_metadata(
      bind_rows_base(parallel_map(seq_len(R), worker, cfg$ncores)), cfg, run_id, R
    )
    progress_message(paste("Experiment 3B n=", n), R, R)
  }
  out <- standardize_replicate_output(bind_rows_base(all_results))
  truth_df <- bind_rows_base(all_truth); truth_df$run_id <- run_id; truth_df$M <- cfg$M; truth_df$K0 <- cfg$K0; truth_df$decay <- cfg$decay
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(out, replicate_path); write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  invisible(out)
}

default_grid_saturation_config <- function() {
  list(
    n = 300L, T = 100L, N = 200L, m_values = c(500L, 1000L, 1500L, 2000L, 3000L),
    replications = 1000L, decay = 1.1, M = 3L, K0 = 5L,
    distribution = "Gaussian", ncores = 1L, master_seed = MASTER_SEED,
    confidence_level = 0.95, certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

run_grid_saturation_experiment <- function(config = list()) {
  cfg <- merge_config(default_grid_saturation_config(), config)
  validate_inference_config(cfg, c("n", "T", "N"))
  if (!length(cfg$m_values)) stop("m_values must not be empty.")
  invisible(lapply(cfg$m_values, validate_positive_integer, name = "m_values entry"))
  if (any(cfg$m_values < 2L)) stop("Every m value must be at least 2.")
  if (!is.numeric(cfg$decay) || length(cfg$decay) != 1L || cfg$decay <= 0) stop("decay must be positive.")
  experiment <- "experiment_3c_grid_saturation"
  run_id <- make_run_id(experiment, cfg, c("n", "T", "N", "m_values", "replications", "decay", "M", "K0", "distribution"))
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  values <- large_domain_spectrum(cfg$n, repeated = FALSE, N = cfg$N, T = cfg$T, decay = cfg$decay)
  truth <- population_spike_truth(values, cfg$M, cfg$n, 1L)
  f <- large_domain_target(cfg$N)

  worker <- function(replication) {
    set.seed(replication_seed(6L, replication, cfg$master_seed))
    coefficients <- sample_large_domain_coefficients(cfg$n, values, cfg$distribution)
    oracle_fit <- gram_pca(coefficients, vectors = TRUE, center = FALSE)
    oracle_est <- dimension_free_inference(oracle_fit$values, 1L, cfg$K0, cfg$n, cfg$confidence_level, cfg$certification_threshold)
    oracle_anchor <- feature_projection(oracle_fit, 1L, c(1, rep(0, cfg$N - 1L)))
    oracle_theta <- orient_projection(feature_projection(oracle_fit, 1L, f), oracle_anchor) / sqrt(oracle_est$hat_r2)

    rows <- lapply(cfg$m_values, function(m) {
      X_grid <- coefficient_to_grid(coefficients, cfg$T, m)
      f_grid <- direction_to_grid(f, cfg$T, m)
      anchor_grid <- direction_to_grid(c(1, rep(0, cfg$N - 1L)), cfg$T, m)
      grid_fit <- gram_pca(X_grid, vectors = TRUE, center = FALSE)
      grid_est <- dimension_free_inference(grid_fit$values, 1L, cfg$K0, cfg$n, cfg$confidence_level, cfg$certification_threshold)
      grid_anchor <- feature_projection(grid_fit, 1L, anchor_grid)
      grid_theta <- orient_projection(feature_projection(grid_fit, 1L, f_grid), grid_anchor) / sqrt(grid_est$hat_r2)
      data.frame(
        experiment = experiment, scenario = "paired_oracle_discrete",
        score_distribution = cfg$distribution, n = cfg$n, p = NA_integer_, T = cfg$T, N = cfg$N, m = m,
        replication = replication, spike_index = 1L, K0 = cfg$K0,
        true_alpha = truth$true_alpha, true_psi = truth$true_psi, true_Delta = truth$true_Delta, true_r2 = truth$true_r2,
        oracle_hat_alpha = oracle_est$hat_alpha, oracle_hat_Delta = oracle_est$hat_Delta,
        oracle_hat_r2 = oracle_est$hat_r2, oracle_hat_theta = oracle_theta,
        grid_hat_alpha = grid_est$hat_alpha, grid_hat_Delta = grid_est$hat_Delta,
        grid_hat_r2 = grid_est$hat_r2, grid_hat_theta = grid_theta,
        hat_lambda = grid_est$hat_lambda, hat_alpha = grid_est$hat_alpha,
        hat_Delta = grid_est$hat_Delta, hat_r2 = grid_est$hat_r2,
        se_alpha = grid_est$se_alpha, se_Delta = grid_est$se_Delta, se_r2 = grid_est$se_r2,
        theta_true = 0.6, theta_raw = NA_real_, theta_corrected = grid_theta,
        discrepancy_alpha = sqrt(cfg$n) * abs(grid_est$hat_alpha - oracle_est$hat_alpha),
        discrepancy_Delta = sqrt(cfg$n) * abs(grid_est$hat_Delta - oracle_est$hat_Delta),
        discrepancy_r2 = sqrt(cfg$n) * abs(grid_est$hat_r2 - oracle_est$hat_r2),
        discrepancy_theta = sqrt(cfg$n) * abs(grid_theta - oracle_theta)
      )
    })
    bind_rows_base(rows)
  }
  out <- add_run_metadata(bind_rows_base(parallel_map(seq_len(cfg$replications), worker, cfg$ncores)), cfg, run_id, cfg$replications)
  out <- standardize_replicate_output(out)
  truth_df <- transform(truth, experiment = experiment, run_id = run_id, n = cfg$n, T = cfg$T, N = cfg$N, M = cfg$M, K0 = cfg$K0, decay = cfg$decay)
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(out, replicate_path); write_csv_atomic(truth_df, truth_path)
  register_experiment_run(experiment, run_id, list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root)
  progress_message("Experiment 3C", cfg$replications, cfg$replications)
  invisible(out)
}
