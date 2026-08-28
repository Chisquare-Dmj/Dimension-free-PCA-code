# Experiment 2B: triangular-array functional inference asymptotics.

FUNCTIONAL_ASYMPTOTIC_P_RULE <- "p=nearest-multiple-10-to-n/3"

functional_asymptotic_p <- function(n) {
  # F2 uses blocks of size 10. Rounding n/3 to the nearest block-compatible
  # dimension keeps p/n asymptotically equal to 1/3 for arbitrary sample sizes.
  as.integer(pmax(10L, 10L * round(as.numeric(n) / 30)))
}

default_functional_asymptotic_config <- function() {
  list(
    n_values = c(150L, 300L, 600L),
    replications = c(`150` = 2000L, `300` = 1500L, `600` = 1000L),
    p_rule = FUNCTIONAL_ASYMPTOTIC_P_RULE,
    N = 200L, decay = 1.1, M = 3L, K0 = 3L,
    scenarios = "F3_Haar",
    distributions = c(F3_Haar = "Uniform"),
    spike_indices = 2L, population_seed = 201L,
    ncores = 1L, master_seed = MASTER_SEED, confidence_level = 0.95,
    certification_threshold = 0.20, root = PROJECT_ROOT
  )
}

prepare_functional_asymptotic_config <- function(config) {
  config$p_rule <- FUNCTIONAL_ASYMPTOTIC_P_RULE
  if (length(config$replications) == length(config$n_values) &&
      (is.null(names(config$replications)) ||
       !all(as.character(config$n_values) %in% names(config$replications)))) {
    names(config$replications) <- as.character(config$n_values)
  }
  config$design_grid <- paste(
    paste0(config$scenarios, "-", unname(config$distributions[config$scenarios])),
    collapse = "-"
  )
  config
}

validate_functional_asymptotic_config <- function(config) {
  validate_positive_integer(config$N, "N")
  validate_positive_integer(config$M, "M")
  validate_positive_integer(config$K0, "K0")
  if (config$K0 != config$M) stop("Experiment 2B baseline requires K0=M; over-deflation is studied separately.")
  if (!length(config$n_values)) stop("n_values must not be empty.")
  invisible(lapply(config$n_values, validate_positive_integer, name = "n_values entry"))
  p_values <- functional_asymptotic_p(config$n_values)
  if ("F2_block" %in% config$scenarios && any(p_values %% 10L != 0L)) stop("Experiment 2B F2 requires p divisible by 10.")
  if (config$K0 >= min(config$n_values)) stop("K0 must be smaller than every n value.")
  if (is.null(names(config$replications)) || !all(as.character(config$n_values) %in% names(config$replications))) {
    stop(
      "replications must be named by every n value; n_values=",
      paste(config$n_values, collapse = ","), "; replication_names=",
      paste(names(config$replications), collapse = ","),
      "; replication_length=", length(config$replications), "."
    )
  }
  invisible(lapply(config$replications[as.character(config$n_values)], validate_positive_integer, name = "replications entry"))
  if (!all(config$scenarios %in% names(config$distributions))) stop("distributions must name every requested scenario.")
  if (!length(config$spike_indices) || any(config$spike_indices < 1L | config$spike_indices > config$M)) stop("spike_indices must lie between 1 and M.")
  if (!is.numeric(config$decay) || length(config$decay) != 1L || config$decay <= 0) stop("decay must be positive.")
  if (config$confidence_level <= 0 || config$confidence_level >= 1) stop("confidence_level must lie strictly between 0 and 1.")
  invisible(config)
}

run_functional_asymptotic_experiment <- function(config = list()) {
  cfg <- prepare_functional_asymptotic_config(
    merge_config(default_functional_asymptotic_config(), config)
  )
  validate_functional_asymptotic_config(cfg)
  experiment <- "experiment_2b_functional_asymptotics"
  run_id <- make_run_id(
    experiment, cfg,
    c("n_values", "replications", "p_rule", "N", "decay", "M", "K0", "spike_indices", "design_grid")
  )
  config_path <- save_config_table(cfg, experiment, run_id, cfg$root)
  all_results <- list(); all_truth <- list(); result_index <- 1L

  for (n_index in seq_along(cfg$n_values)) {
    n <- as.integer(cfg$n_values[n_index])
    p <- functional_asymptotic_p(n)
    replications <- as.integer(cfg$replications[[as.character(n)]])
    for (scenario in cfg$scenarios) {
      distribution <- unname(cfg$distributions[[scenario]])
      population <- build_functional_population(
        scenario, p, cfg$N, cfg$decay, cfg$population_seed
      )
      truth <- population_spike_truth(population$values, cfg$M, n, cfg$spike_indices)
      derivatives <- lapply(cfg$spike_indices, function(j) {
        values <- population_derivatives(population$values, cfg$M, n, j)
        names(values) <- paste0("true_", names(values))
        as.data.frame(as.list(values))
      })
      derivatives <- bind_rows_base(derivatives); derivatives$spike_index <- cfg$spike_indices
      truth <- merge(truth, derivatives, by = "spike_index", sort = FALSE)
      scenario_truth <- transform(
        truth, experiment = experiment, run_id = run_id,
        scenario = scenario, score_distribution = distribution,
        n = n, p = p, N = cfg$N, M = cfg$M, K0 = cfg$K0,
        decay = cfg$decay, population_seed = cfg$population_seed
      )
      all_truth[[result_index]] <- scenario_truth
      target_directions <- lapply(
        cfg$spike_indices, function(j) functional_population_direction(population, j)
      )

      worker <- function(replication) {
        set.seed(replication_seed(9L, replication, cfg$master_seed))
        X <- sample_functional_coefficients(n, population, distribution)
        fit <- gram_pca(X, vectors = TRUE, center = FALSE)
        rows <- lapply(seq_along(cfg$spike_indices), function(index) {
          j <- cfg$spike_indices[index]
          estimate <- dimension_free_inference(
            fit$values, j, cfg$K0, n,
            cfg$confidence_level, cfg$certification_threshold
          )
          row <- attach_truth_and_coverage(
            estimate, scenario_truth[scenario_truth$spike_index == j, ]
          )
          fpca_interval <- fpca_eigenvalue_interval(fit, j, cfg$confidence_level)
          row <- cbind(row, fpca_interval)
          row$fpca_alpha_covered <- as.integer(
            row$fpca_alpha_ci_lower <= row$true_alpha &
              row$true_alpha <= row$fpca_alpha_ci_upper
          )
          overlap <- feature_projection(fit, j, target_directions[[index]])
          row$experiment <- experiment; row$scenario <- scenario
          row$score_distribution <- distribution
          row$n <- n; row$p <- p; row$T <- NA_integer_; row$N <- cfg$N; row$m <- NA_integer_
          row$replication <- replication; row$actual_signal_overlap <- overlap^2
          row
        })
        bind_rows_base(rows)
      }
      all_results[[result_index]] <- add_run_metadata(
        bind_rows_base(parallel_map(seq_len(replications), worker, cfg$ncores)),
        cfg, run_id, replications, cfg$population_seed
      )
      progress_message(
        paste("Experiment 2B", scenario, distribution, "n=", n),
        replications, replications
      )
      result_index <- result_index + 1L
    }
  }

  output <- standardize_replicate_output(bind_rows_base(all_results))
  truth_df <- bind_rows_base(all_truth)
  replicate_path <- run_artifact_path("replicate", run_id, "replicates", "csv", cfg$root)
  truth_path <- run_artifact_path("truth", run_id, "truth", "csv", cfg$root)
  write_csv_atomic(output, replicate_path); write_csv_atomic(truth_df, truth_path)
  register_experiment_run(
    experiment, run_id,
    list(replicate = replicate_path, truth = truth_path, config = config_path), cfg$root
  )
  invisible(output)
}
