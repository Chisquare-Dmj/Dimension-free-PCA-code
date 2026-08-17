# Saved-result summaries and LaTeX table writers.

group_apply <- function(data, group_columns, FUN) {
  key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  pieces <- split(data, key, drop = TRUE)
  rows <- lapply(pieces, function(piece) cbind(piece[1, group_columns, drop = FALSE], FUN(piece)))
  bind_rows_base(rows)
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)

summarize_panel_phase <- function(data) {
  group_apply(data, c("scenario", "spike_index"), function(x) {
    data.frame(
      replications = length(unique(x$replication)),
      true_alpha = x$true_alpha[1],
      true_psi = x$true_psi[1],
      true_Delta = x$true_Delta[1],
      true_r2 = x$true_r2[1],
      mean_hat_lambda = mean_or_na(x$hat_lambda),
      sd_hat_lambda = sd_or_na(x$hat_lambda),
      mean_sample_gap_to_lambda4 = mean_or_na(x$sample_gap_to_lambda4)
    )
  })
}

summarize_inference_performance <- function(data) {
  group_apply(data, c("scenario", "score_distribution", "spike_index"), function(x) {
    data.frame(
      replications = nrow(x),
      true_alpha = x$true_alpha[1],
      alpha_relative_bias_percent = 100 * mean((x$hat_alpha - x$true_alpha) / x$true_alpha, na.rm = TRUE),
      alpha_relative_rmse_percent = 100 * sqrt(mean(((x$hat_alpha - x$true_alpha) / x$true_alpha)^2, na.rm = TRUE)),
      empirical_sd_alpha = sd_or_na(x$hat_alpha),
      mean_se_alpha = mean_or_na(x$se_alpha),
      alpha_sd_over_mean_se = sd_or_na(x$hat_alpha) / mean_or_na(x$se_alpha),
      alpha_coverage_percent = 100 * mean_or_na(x$cover_alpha),
      true_Delta = x$true_Delta[1],
      Delta_bias = mean(x$hat_Delta - x$true_Delta, na.rm = TRUE),
      Delta_rmse = sqrt(mean((x$hat_Delta - x$true_Delta)^2, na.rm = TRUE)),
      empirical_sd_Delta = sd_or_na(x$hat_Delta),
      mean_se_Delta = mean_or_na(x$se_Delta),
      Delta_sd_over_mean_se = sd_or_na(x$hat_Delta) / mean_or_na(x$se_Delta),
      Delta_coverage_percent = 100 * mean_or_na(x$cover_Delta),
      true_r2 = x$true_r2[1],
      r2_bias = mean(x$hat_r2 - x$true_r2, na.rm = TRUE),
      r2_rmse = sqrt(mean((x$hat_r2 - x$true_r2)^2, na.rm = TRUE)),
      empirical_sd_r2 = sd_or_na(x$hat_r2),
      mean_se_r2 = mean_or_na(x$se_r2),
      r2_sd_over_mean_se = sd_or_na(x$hat_r2) / mean_or_na(x$se_r2),
      r2_coverage_percent = 100 * mean_or_na(x$cover_r2),
      mean_actual_signal_overlap = mean_or_na(x$actual_signal_overlap)
    )
  })
}

summarize_asymptotic_inference <- function(data) {
  group_apply(data, c("scenario", "score_distribution", "n", "p", "spike_index"), function(x) {
    alpha_relative_error <- (x$hat_alpha - x$true_alpha) / x$true_alpha
    alpha_relative_rmse <- 100 * sqrt(mean(alpha_relative_error^2, na.rm = TRUE))
    Delta_rmse <- sqrt(mean((x$hat_Delta - x$true_Delta)^2, na.rm = TRUE))
    r2_rmse <- sqrt(mean((x$hat_r2 - x$true_r2)^2, na.rm = TRUE))
    data.frame(
      replications = nrow(x),
      true_alpha = x$true_alpha[1], true_Delta = x$true_Delta[1], true_r2 = x$true_r2[1],
      alpha_relative_bias_percent = 100 * mean(alpha_relative_error, na.rm = TRUE),
      alpha_relative_rmse_percent = alpha_relative_rmse,
      alpha_scaled_relative_rmse_percent = sqrt(x$n[1]) * alpha_relative_rmse,
      empirical_sd_alpha = sd_or_na(x$hat_alpha),
      alpha_sd_over_mean_se = sd_or_na(x$hat_alpha) / mean_or_na(x$se_alpha),
      alpha_coverage_percent = 100 * mean_or_na(x$cover_alpha),
      mean_se_alpha = mean_or_na(x$se_alpha),
      sqrt_n_mean_se_alpha = sqrt(x$n[1]) * mean_or_na(x$se_alpha),
      Delta_bias = mean(x$hat_Delta - x$true_Delta, na.rm = TRUE),
      Delta_rmse = Delta_rmse,
      Delta_scaled_rmse = sqrt(x$n[1]) * Delta_rmse,
      empirical_sd_Delta = sd_or_na(x$hat_Delta),
      Delta_sd_over_mean_se = sd_or_na(x$hat_Delta) / mean_or_na(x$se_Delta),
      Delta_coverage_percent = 100 * mean_or_na(x$cover_Delta),
      mean_se_Delta = mean_or_na(x$se_Delta),
      sqrt_n_mean_se_Delta = sqrt(x$n[1]) * mean_or_na(x$se_Delta),
      r2_bias = mean(x$hat_r2 - x$true_r2, na.rm = TRUE),
      r2_rmse = r2_rmse,
      r2_scaled_rmse = sqrt(x$n[1]) * r2_rmse,
      empirical_sd_r2 = sd_or_na(x$hat_r2),
      r2_sd_over_mean_se = sd_or_na(x$hat_r2) / mean_or_na(x$se_r2),
      r2_coverage_percent = 100 * mean_or_na(x$cover_r2),
      mean_se_r2 = mean_or_na(x$se_r2),
      sqrt_n_mean_se_r2 = sqrt(x$n[1]) * mean_or_na(x$se_r2)
    )
  })
}

summarize_phase_margin <- function(data) {
  group_apply(data, c("phase_region", "target_Delta"), function(x) {
    data.frame(
      replications = nrow(x), true_Delta = x$true_Delta[1],
      mean_hat_Delta = mean_or_na(x$hat_Delta),
      Delta_bias = mean(x$hat_Delta - x$true_Delta, na.rm = TRUE),
      Delta_coverage_percent = 100 * mean_or_na(x$cover_Delta),
      certification_probability = mean_or_na(x$phase_certified),
      naive_relative_bias = mean_or_na(x$relative_bias_naive),
      corrected_relative_bias = mean_or_na(x$relative_bias_corrected),
      scaled_alpha_variance = x$n[1] * var(x$hat_alpha, na.rm = TRUE) / x$true_alpha[1]^2,
      theoretical_scaled_variance = 2 / x$true_Delta[1],
      mean_alpha_ci_length = mean_or_na(x$ci_length_alpha)
    )
  })
}

summarize_direction <- function(data) {
  group_apply(data, "n", function(x) {
    data.frame(
      theta_true = x$theta_true[1],
      raw_mean = mean_or_na(x$theta_raw), raw_bias = mean(x$theta_raw - x$theta_true, na.rm = TRUE),
      raw_empirical_sd = sd_or_na(x$theta_raw),
      raw_rmse = sqrt(mean((x$theta_raw - x$theta_true)^2, na.rm = TRUE)),
      corrected_mean = mean_or_na(x$theta_corrected), corrected_bias = mean(x$theta_corrected - x$theta_true, na.rm = TRUE),
      corrected_empirical_sd = sd_or_na(x$theta_corrected),
      corrected_rmse = sqrt(mean((x$theta_corrected - x$theta_true)^2, na.rm = TRUE))
    )
  })
}

summarize_repeated <- function(data) {
  group_apply(data, "n", function(x) {
    data.frame(
      Q_true = x$Q_true[1],
      raw_mean = mean_or_na(x$Q_raw), raw_bias = mean(x$Q_raw - x$Q_true, na.rm = TRUE),
      raw_empirical_sd = sd_or_na(x$Q_raw),
      raw_rmse = sqrt(mean((x$Q_raw - x$Q_true)^2, na.rm = TRUE)),
      corrected_mean = mean_or_na(x$Q_corrected),
      corrected_bias = mean(x$Q_corrected - x$Q_true, na.rm = TRUE),
      corrected_empirical_sd = sd_or_na(x$Q_corrected),
      corrected_rmse = sqrt(mean((x$Q_corrected - x$Q_true)^2, na.rm = TRUE)),
      alignment_mean = mean_or_na(x$individual_alignment_repeated),
      alignment_sd = sd_or_na(x$individual_alignment_repeated)
    )
  })
}

# Compare two estimators of the same scalar population target. These summaries
# deliberately keep finite-sample Monte Carlo evidence separate from any
# asymptotic statement, which must be justified by the accompanying theory.
summarize_two_method_target <- function(data, group_columns, truth_column,
                                        fpca_column, proposed_column,
                                        include_relative_bias = TRUE) {
  rows <- list()
  for (method in c("FPCA", "Proposed")) {
    estimate_column <- if (method == "FPCA") fpca_column else proposed_column
    method_summary <- group_apply(data, group_columns, function(x) {
      truth <- x[[truth_column]][1]
      estimate <- x[[estimate_column]]
      row <- data.frame(
        replications = sum(is.finite(estimate)),
        truth = truth,
        mean = mean_or_na(estimate),
        bias = mean_or_na(estimate - truth),
        empirical_sd = sd_or_na(estimate),
        rmse = sqrt(mean((estimate - truth)^2, na.rm = TRUE))
      )
      if (include_relative_bias) {
        row$relative_bias <- row$bias / truth
        row$relative_bias_percent <- 100 * row$relative_bias
      }
      row
    })
    method_summary$method <- method
    rows[[method]] <- method_summary
  }
  result <- bind_rows_base(rows)
  result[, c(group_columns, "method", setdiff(names(result), c(group_columns, "method"))), drop = FALSE]
}

summarize_large_domain_eigenvalue_comparison <- function(data) {
  summarize_two_method_target(
    data, "n", "true_alpha", "hat_lambda", "hat_alpha",
    include_relative_bias = TRUE
  )
}

summarize_large_domain_direction_comparison <- function(data) {
  summarize_two_method_target(
    data, "n", "theta_true", "theta_raw", "theta_corrected",
    include_relative_bias = FALSE
  )
}

summarize_large_domain_eigenspace_comparison <- function(data) {
  summarize_two_method_target(
    data, "n", "Q_true", "Q_raw", "Q_corrected",
    include_relative_bias = FALSE
  )
}

summarize_functional_eigenvalue_comparison <- function(data) {
  summarize_two_method_target(
    data, c("scenario", "score_distribution", "n", "p", "spike_index"),
    "true_alpha", "hat_lambda", "hat_alpha", include_relative_bias = TRUE
  )
}

summarize_functional_reliability_check <- function(data) {
  group_apply(
    data,
    c("scenario", "score_distribution", "n", "p", "spike_index"),
    function(x) data.frame(
      replications = nrow(x),
      true_r2 = x$true_r2[1],
      mean_empirical_squared_overlap = mean_or_na(x$actual_signal_overlap),
      mean_proposed_r2 = mean_or_na(x$hat_r2),
      mean_difference = mean_or_na(x$hat_r2) - mean_or_na(x$actual_signal_overlap)
    )
  )
}

summarize_large_domain_individual_alignment <- function(data) {
  group_apply(data, "n", function(x) data.frame(
    replications = nrow(x),
    mean_squared_alignment = mean_or_na(x$individual_alignment_repeated),
    sd_squared_alignment = sd_or_na(x$individual_alignment_repeated)
  ))
}

write_named_simulation_summary <- function(data, filename, root = PROJECT_ROOT) {
  path <- result_path("summary", filename, root)
  write_csv_atomic(data, path)
  path
}

generate_fpca_comparison_summaries <- function(root = PROJECT_ROOT, run_ids = list()) {
  experiments <- c(
    direction = "experiment_3a_simple_direction",
    eigenspace = "experiment_3b_repeated_eigenspace",
    functional = "experiment_2b_functional_asymptotics"
  )
  source_ids <- vapply(experiments, function(experiment) {
    resolved_run_id(experiment, root, run_ids[[experiment]])
  }, character(1))
  direction <- read_required_csv(resolve_run_artifact(
    experiments[["direction"]], "replicate", root,
    run_ids[[experiments[["direction"]]]]
  ))
  eigenspace <- read_required_csv(resolve_run_artifact(
    experiments[["eigenspace"]], "replicate", root,
    run_ids[[experiments[["eigenspace"]]]]
  ))
  functional <- read_required_csv(resolve_run_artifact(
    experiments[["functional"]], "replicate", root,
    run_ids[[experiments[["functional"]]]]
  ))

  outputs <- list(
    large_domain_fpca_eigenvalue = summarize_large_domain_eigenvalue_comparison(direction),
    large_domain_direction_comparison = summarize_large_domain_direction_comparison(direction),
    large_domain_eigenspace_comparison = summarize_large_domain_eigenspace_comparison(eigenspace),
    large_domain_individual_alignment = summarize_large_domain_individual_alignment(eigenspace),
    functional_eigenvalue_fpca_vs_proposed = summarize_functional_eigenvalue_comparison(functional),
    functional_reliability_check = summarize_functional_reliability_check(functional)
  )
  outputs$large_domain_fpca_eigenvalue$source_run_id <- source_ids[["direction"]]
  outputs$large_domain_direction_comparison$source_run_id <- source_ids[["direction"]]
  outputs$large_domain_eigenspace_comparison$source_run_id <- source_ids[["eigenspace"]]
  outputs$large_domain_individual_alignment$source_run_id <- source_ids[["eigenspace"]]
  outputs$functional_eigenvalue_fpca_vs_proposed$source_run_id <- source_ids[["functional"]]
  outputs$functional_reliability_check$source_run_id <- source_ids[["functional"]]

  paths <- mapply(
    function(data, name) write_named_simulation_summary(data, paste0(name, ".csv"), root),
    outputs, names(outputs), SIMPLIFY = FALSE
  )
  functional_replicates <- functional[, c(
    "n", "replication", "true_alpha", "hat_lambda", "hat_alpha",
    "hat_Delta", "hat_r2", "actual_signal_overlap"
  ), drop = FALSE]
  names(functional_replicates) <- c(
    "n", "replication", "alpha_true", "lambda_fpca", "alpha_proposed",
    "Delta_hat", "r2_hat", "true_overlap"
  )
  functional_replicates$source_run_id <- source_ids[["functional"]]
  paths$functional_eigenvalue_fpca_vs_proposed_replicates <- result_path(
    "replicate", "functional_eigenvalue_fpca_vs_proposed_replicates.csv", root
  )
  write_csv_atomic(functional_replicates,
                   paths$functional_eigenvalue_fpca_vs_proposed_replicates)
  combined_parts <- list()
  for (name in c(
    "large_domain_fpca_eigenvalue", "large_domain_direction_comparison",
    "large_domain_eigenspace_comparison", "functional_eigenvalue_fpca_vs_proposed"
  )) {
    piece <- outputs[[name]]
    piece$comparison <- name
    combined_parts[[name]] <- piece
  }
  combined <- bind_rows_base(combined_parts)
  paths$simulation_fpca_comparison_summary <- write_named_simulation_summary(
    combined, "simulation_fpca_comparison_summary.csv", root
  )
  invisible(list(data = outputs, paths = paths, source_run_ids = source_ids))
}

summarize_grid <- function(data) {
  quantities <- c("alpha", "Delta", "r2", "theta")
  rows <- list(); index <- 1L
  for (m in sort(unique(data$m))) {
    piece <- data[data$m == m, ]
    for (quantity in quantities) {
      values <- piece[[paste0("discrepancy_", quantity)]]
      rows[[index]] <- data.frame(
        m = m, m_over_n = m / piece$n[1], quantity = quantity,
        median = median(values, na.rm = TRUE),
        q90 = unname(quantile(values, 0.90, na.rm = TRUE))
      )
      index <- index + 1L
    }
  }
  bind_rows_base(rows)
}

summarize_k0 <- function(data) {
  group_apply(data, c("score_distribution", "n", "p", "K0", "spike_index"), function(x) {
    alpha_error <- x$overdeflation_error_alpha
    Delta_error <- x$overdeflation_error_Delta
    r2_error <- x$overdeflation_error_r2
    data.frame(
      replications = length(unique(x$replication)),
      overdeflation_excess = x$overdeflation_excess[1],
      true_alpha = x$true_alpha[1],
      mean_hat_alpha = mean_or_na(x$hat_alpha),
      coverage_alpha = mean_or_na(x$cover_alpha),
      alpha_overdeflation_bias = mean_or_na(alpha_error),
      alpha_overdeflation_rmse = sqrt(mean(alpha_error^2, na.rm = TRUE)),
      alpha_sqrt_n_overdeflation_rmse = sqrt(x$n[1]) * sqrt(mean(alpha_error^2, na.rm = TRUE)),
      alpha_n_overdeflation_rmse = x$n[1] * sqrt(mean(alpha_error^2, na.rm = TRUE)),
      true_Delta = x$true_Delta[1],
      mean_hat_Delta = mean_or_na(x$hat_Delta),
      coverage_Delta = mean_or_na(x$cover_Delta),
      Delta_overdeflation_bias = mean_or_na(Delta_error),
      Delta_overdeflation_rmse = sqrt(mean(Delta_error^2, na.rm = TRUE)),
      Delta_sqrt_n_overdeflation_rmse = sqrt(x$n[1]) * sqrt(mean(Delta_error^2, na.rm = TRUE)),
      Delta_n_overdeflation_rmse = x$n[1] * sqrt(mean(Delta_error^2, na.rm = TRUE)),
      true_r2 = x$true_r2[1],
      mean_hat_r2 = mean_or_na(x$hat_r2),
      coverage_r2 = mean_or_na(x$cover_r2),
      r2_overdeflation_bias = mean_or_na(r2_error),
      r2_overdeflation_rmse = sqrt(mean(r2_error^2, na.rm = TRUE)),
      r2_sqrt_n_overdeflation_rmse = sqrt(x$n[1]) * sqrt(mean(r2_error^2, na.rm = TRUE)),
      r2_n_overdeflation_rmse = x$n[1] * sqrt(mean(r2_error^2, na.rm = TRUE))
    )
  })
}

full_cumulant_sd_multiplier <- function(fourth_cumulant, true_Delta) {
  variance_ratio <- 1 + fourth_cumulant * true_Delta / 2
  if (any(variance_ratio < 0, na.rm = TRUE)) {
    stop("The full-cumulant variance multiplier must be nonnegative.")
  }
  sqrt(variance_ratio)
}

summarize_universality <- function(data) {
  group_apply(data, c("score_distribution", "spike_index"), function(x) {
    # Theorem 4.1 gives one common first-order fluctuation. The fourth-cumulant
    # term therefore rescales all three observable universal standard errors by
    # the same population factor.
    multiplier <- full_cumulant_sd_multiplier(x$fourth_cumulant, x$true_Delta)
    full_se_alpha <- x$se_alpha * multiplier
    full_se_Delta <- x$se_Delta * multiplier
    full_se_r2 <- x$se_r2 * multiplier
    z <- qnorm(0.975)
    data.frame(
      empirical_sd_alpha = sd_or_na(x$hat_alpha),
      mean_universal_se_alpha = mean_or_na(x$se_alpha),
      mean_full_se_alpha = mean_or_na(full_se_alpha),
      universal_sd_ratio_alpha = mean_or_na(x$se_alpha) / sd_or_na(x$hat_alpha),
      full_sd_ratio_alpha = mean_or_na(full_se_alpha) / sd_or_na(x$hat_alpha),
      empirical_sd_Delta = sd_or_na(x$hat_Delta),
      mean_universal_se_Delta = mean_or_na(x$se_Delta),
      mean_full_se_Delta = mean_or_na(full_se_Delta),
      universal_sd_ratio_Delta = mean_or_na(x$se_Delta) / sd_or_na(x$hat_Delta),
      full_sd_ratio_Delta = mean_or_na(full_se_Delta) / sd_or_na(x$hat_Delta),
      empirical_sd_r2 = sd_or_na(x$hat_r2),
      mean_universal_se_r2 = mean_or_na(x$se_r2),
      mean_full_se_r2 = mean_or_na(full_se_r2),
      universal_sd_ratio_r2 = mean_or_na(x$se_r2) / sd_or_na(x$hat_r2),
      full_sd_ratio_r2 = mean_or_na(full_se_r2) / sd_or_na(x$hat_r2),
      universal_alpha_coverage = mean_or_na(x$cover_alpha),
      full_alpha_coverage = mean_or_na(abs(x$hat_alpha - x$true_alpha) <= z * full_se_alpha),
      universal_Delta_coverage = mean_or_na(x$cover_Delta),
      full_Delta_coverage = mean_or_na(abs(x$hat_Delta - x$true_Delta) <= z * full_se_Delta),
      universal_r2_coverage = mean_or_na(x$cover_r2),
      full_r2_coverage = mean_or_na(abs(x$hat_r2 - x$true_r2) <= z * full_se_r2)
    )
  })
}

escape_latex <- function(x) {
  x <- gsub("_", "\\_", as.character(x), fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x
}

write_latex_table <- function(data, path, digits = 3L, caption = NULL, label = NULL) {
  formatted <- data
  numeric_columns <- vapply(formatted, is.numeric, logical(1))
  formatted[numeric_columns] <- lapply(formatted[numeric_columns], function(x) ifelse(is.na(x), "", formatC(x, digits = digits, format = "f")))
  formatted[!numeric_columns] <- lapply(formatted[!numeric_columns], escape_latex)
  alignment <- paste0("l", paste(rep("r", ncol(formatted) - 1L), collapse = ""))
  lines <- c("\\begin{table}[htbp]", "\\centering")
  if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", caption, "}"))
  if (!is.null(label)) lines <- c(lines, paste0("\\label{", label, "}"))
  lines <- c(lines, paste0("\\begin{tabular}{", alignment, "}"), "\\toprule")
  lines <- c(lines, paste(escape_latex(names(formatted)), collapse = " & "), "\\\\", "\\midrule")
  body <- apply(formatted, 1L, function(row) paste0(paste(row, collapse = " & "), " \\\\"))
  lines <- c(lines, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

write_asymptotic_inference_table <- function(data, path) {
  panel_specifications <- list(
    "Panel A: Population spike alpha (relative bias and RMSE in percent)" = c("alpha_relative_bias_percent", "alpha_relative_rmse_percent", "alpha_scaled_relative_rmse_percent", "empirical_sd_alpha", "mean_se_alpha", "alpha_sd_over_mean_se", "alpha_coverage_percent"),
    "Panel B: Phase margin Delta" = c("Delta_bias", "Delta_rmse", "Delta_scaled_rmse", "empirical_sd_Delta", "mean_se_Delta", "Delta_sd_over_mean_se", "Delta_coverage_percent"),
    "Panel C: Reliability r2" = c("r2_bias", "r2_rmse", "r2_scaled_rmse", "empirical_sd_r2", "mean_se_r2", "r2_sd_over_mean_se", "r2_coverage_percent")
  )
  lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Functional inference along the high-complexity asymptotic sequence}",
    "\\label{tab:functional-asymptotics}",
    "\\begin{tabular}{llrrrrrrrrrr}", "\\toprule",
    "Scenario & Score & $n$ & $p$ & Spike & Bias & RMSE & $\\sqrt{n}$ RMSE & Emp. SD & Mean SE & SD/SE & Coverage (\\%) \\\\",
    "\\midrule"
  )
  for (panel_name in names(panel_specifications)) {
    columns <- panel_specifications[[panel_name]]
    lines <- c(lines, paste0("\\multicolumn{12}{l}{\\textit{", panel_name, "}} \\\\"))
    for (i in seq_len(nrow(data))) {
      values <- c(
        escape_latex(data$scenario[i]), escape_latex(data$score_distribution[i]),
        data$n[i], data$p[i], data$spike_index[i],
        vapply(data[i, columns, drop = FALSE], function(x) formatC(x, digits = 3L, format = "f"), character(1))
      )
      lines <- c(lines, paste0(paste(values, collapse = " & "), " \\\\"))
    }
    lines <- c(lines, "\\addlinespace")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

write_k0_robustness_table <- function(data, path) {
  panel_specifications <- list(
    "Panel A: Population spike alpha" = c(
      "alpha_overdeflation_bias", "alpha_overdeflation_rmse",
      "alpha_sqrt_n_overdeflation_rmse", "alpha_n_overdeflation_rmse",
      "coverage_alpha"
    ),
    "Panel B: Phase margin Delta" = c(
      "Delta_overdeflation_bias", "Delta_overdeflation_rmse",
      "Delta_sqrt_n_overdeflation_rmse", "Delta_n_overdeflation_rmse",
      "coverage_Delta"
    ),
    "Panel C: Reliability r2" = c(
      "r2_overdeflation_bias", "r2_overdeflation_rmse",
      "r2_sqrt_n_overdeflation_rmse", "r2_n_overdeflation_rmse",
      "coverage_r2"
    )
  )
  lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Sensitivity to the number of leading components removed from the bulk estimate}",
    "\\label{tab:k0-robustness}",
    "\\begin{tabular}{lrrrrrrrrr}", "\\toprule",
    "Score & $n$ & $p$ & $K_0$ & $K_0-M$ & Bias & RMSE & $\\sqrt{n}$ RMSE & $n$ RMSE & Coverage (\\%) \\\\",
    "\\midrule"
  )
  for (panel_name in names(panel_specifications)) {
    columns <- panel_specifications[[panel_name]]
    lines <- c(lines, paste0("\\multicolumn{10}{l}{\\textit{", panel_name, "}} \\\\"))
    for (i in seq_len(nrow(data))) {
      metrics <- as.numeric(unlist(data[i, columns, drop = FALSE], use.names = FALSE))
      metrics[5] <- 100 * metrics[5]
      values <- c(
        escape_latex(data$score_distribution[i]),
        data$n[i], data$p[i], data$K0[i], data$overdeflation_excess[i],
        formatC(metrics, digits = 3L, format = "f")
      )
      lines <- c(lines, paste0(paste(values, collapse = " & "), " \\\\"))
    }
    lines <- c(lines, "\\addlinespace")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

write_universality_table <- function(data, path) {
  panel_specifications <- list(
    "Panel A: Population spike alpha" = c(
      "empirical_sd_alpha", "mean_universal_se_alpha",
      "universal_sd_ratio_alpha", "mean_full_se_alpha", "full_sd_ratio_alpha"
    ),
    "Panel B: Phase margin Delta" = c(
      "empirical_sd_Delta", "mean_universal_se_Delta",
      "universal_sd_ratio_Delta", "mean_full_se_Delta", "full_sd_ratio_Delta"
    ),
    "Panel C: Reliability r2" = c(
      "empirical_sd_r2", "mean_universal_se_r2",
      "universal_sd_ratio_r2", "mean_full_se_r2", "full_sd_ratio_r2"
    )
  )
  lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{Universality and full fourth-cumulant standard-deviation calibration}",
    "\\label{tab:universality-sd}",
    "\\begin{tabular}{lrrrrrr}", "\\toprule",
    "Score & Spike & Emp. SD & Universal SD & $R_U$ & Full SD & $R_F$ \\\\",
    "\\midrule"
  )
  for (panel_name in names(panel_specifications)) {
    columns <- panel_specifications[[panel_name]]
    lines <- c(lines, paste0("\\multicolumn{7}{l}{\\textit{", panel_name, "}} \\\\"))
    for (i in seq_len(nrow(data))) {
      metrics <- as.numeric(unlist(data[i, columns, drop = FALSE], use.names = FALSE))
      values <- c(
        escape_latex(data$score_distribution[i]), data$spike_index[i],
        formatC(metrics, digits = 3L, format = "f")
      )
      lines <- c(lines, paste0(paste(values, collapse = " & "), " \\\\"))
    }
    lines <- c(lines, "\\addlinespace")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

summary_specifications <- function() {
  list(
    experiment_1a_panel_phase = list(FUN = summarize_panel_phase, output = "panel_phase_transition"),
    experiment_1b_phase_margin = list(FUN = summarize_phase_margin, output = "phase_margin_certification"),
    experiment_2_functional_inference = list(FUN = summarize_inference_performance, output = "inference_performance"),
    experiment_2b_functional_asymptotics = list(FUN = summarize_asymptotic_inference, output = "functional_asymptotics", writer = write_asymptotic_inference_table),
    experiment_3a_simple_direction = list(FUN = summarize_direction, output = "simple_direction"),
    experiment_3b_repeated_eigenspace = list(FUN = summarize_repeated, output = "repeated_eigenspace"),
    experiment_3c_grid_saturation = list(FUN = summarize_grid, output = "grid_saturation"),
    robustness_k0 = list(FUN = summarize_k0, output = "k0_robustness", writer = write_k0_robustness_table),
    robustness_universality = list(FUN = summarize_universality, output = "universality_stress", writer = write_universality_table)
  )
}

generate_experiment_summary <- function(experiment, root = PROJECT_ROOT, run_id = NULL) {
  specification <- summary_specifications()[[experiment]]
  if (is.null(specification)) stop("No summary specification for experiment: ", experiment)
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source_run_id <- resolved_run_id(experiment, root, run_id)
  summary <- specification$FUN(data)
  summary$source_run_id <- source_run_id
  stem <- paste0(source_run_id, "__", specification$output)
  csv_path <- result_path("summary", paste0(stem, ".csv"), root)
  table_path <- result_path("table", paste0(stem, ".tex"), root)
  write_csv_atomic(summary, csv_path)
  if (!is.null(specification$writer)) {
    specification$writer(summary, table_path)
  } else {
    write_latex_table(summary, table_path, caption = gsub("_", " ", stem, fixed = TRUE))
  }
  invisible(list(data = summary, csv = csv_path, table = table_path, run_id = source_run_id))
}

generate_all_summaries <- function(root = PROJECT_ROOT, run_ids = list(), experiments = names(summary_specifications())) {
  outputs <- lapply(experiments, function(experiment) {
    generate_experiment_summary(experiment, root, run_ids[[experiment]])
  })
  names(outputs) <- experiments
  invisible(outputs)
}
