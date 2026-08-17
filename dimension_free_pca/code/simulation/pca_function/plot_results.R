# Publication figures generated from saved CSV results. Each exported PDF is a
# single statistical panel; LaTeX controls any later multi-panel arrangement.

FIGURE3_SCENARIO <- "F3_Haar"
FIGURE3_DISTRIBUTION <- "Uniform"
FIGURE3_SPIKE_INDEX <- 2L
FIGURE3_N <- 300L

publication_par <- function() {
  par(
    mar = c(5.1, 5.3, 1.0, 1.0),
    cex = 1.08, cex.axis = 1.14, cex.lab = 1.25,
    mgp = c(3.2, 0.9, 0), tcl = -0.28, family = "sans"
  )
}

append_figure_manifest <- function(row, root = PROJECT_ROOT) {
  path <- result_path("summary", "simulation_figure_output_manifest.csv", root)
  manifest <- if (file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else data.frame()
  if (nrow(manifest)) manifest <- manifest[manifest$filename != row$filename, , drop = FALSE]
  manifest <- if (nrow(manifest)) bind_rows_base(list(manifest, row)) else row
  write_csv_atomic(manifest, path)
  refresh_figure_output_summary(root)
  invisible(path)
}

save_publication_figure <- function(name, plot_data, draw, root = PROJECT_ROOT,
                                    width = 6.4, height = 5.3,
                                    source_runs = character(), setting = "",
                                    quantity = "", comparison = "") {
  data_path <- result_path("summary", paste0(name, ".csv"), root)
  write_csv_atomic(plot_data, data_path)
  pdf_path <- result_path("figure", paste0(name, ".pdf"), root)
  png_path <- result_path("figure", paste0(name, ".png"), root)
  render <- function() {
    old <- publication_par()
    on.exit(par(old), add = TRUE)
    draw(read_required_csv(data_path))
  }
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  render(); grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = 200L)
  render(); grDevices::dev.off()
  source_text <- paste(paste(names(source_runs), unname(source_runs), sep = "="), collapse = ";")
  append_figure_manifest(data.frame(
    filename = relative_to_root(pdf_path, root),
    plot_data_file = relative_to_root(data_path, root),
    setting = setting, quantity = quantity, comparison = comparison,
    method_scope = figure_method_scope(quantity, comparison),
    source_run_ids = source_text, stringsAsFactors = FALSE
  ), root)
  invisible(c(pdf = pdf_path, png = png_path, data = data_path))
}

safe_density <- function(x) density(x[is.finite(x)], n = 512)

density_rows <- function(values, curve) {
  d <- safe_density(values)
  data.frame(curve = curve, x = d$x, density = d$y)
}

plot_phase_spectrum <- function(spikes, bulk, scenario, filename, source, root) {
  rows <- list(density_rows(bulk$bulk_eigenvalue[bulk$scenario == scenario], "Bulk spectrum"))
  spike3 <- spikes[spikes$scenario == scenario & spikes$spike_index == 3L, ]
  rows[[2L]] <- density_rows(spike3$hat_lambda - spike3$sample_gap_to_lambda4,
                            "Bulk edge: lambda_4")
  for (j in 1:3) {
    rows[[j + 2L]] <- density_rows(
      spikes$hat_lambda[spikes$scenario == scenario & spikes$spike_index == j],
      paste("Spike", j)
    )
  }
  plot_data <- bind_rows_base(rows)
  palette <- c("grey45", "black", "#0072B2", "#D55E00", "#009E73")
  lty <- c(1, 2, 1, 1, 1)
  save_publication_figure(
    filename, plot_data,
    function(x) {
      curves <- unique(x$curve)
      plot(NA, xlim = range(x$x), ylim = c(0, 1.06 * max(x$density)),
           xlab = "Sample eigenvalue", ylab = "Density")
      for (i in seq_along(curves)) {
        piece <- x[x$curve == curves[i], ]
        lines(piece$x, piece$density, col = palette[i], lty = lty[i], lwd = 2)
      }
      legend("topright", curves, col = palette, lty = lty, lwd = 2, bty = "n")
    }, root, source_runs = setNames(source, "experiment_1a_panel_phase"),
    setting = "n=300; two distant spikes and one close spike",
    quantity = "sample spectral density", comparison = "bulk edge and three sample spikes"
  )
}

plot_panel_phase_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "experiment_1a_panel_phase"
  spikes <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  bulk <- read_required_csv(resolve_run_artifact(experiment, "bulk", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  scenarios <- unique(spikes$scenario)
  if (length(scenarios) < 2L) stop("Experiment 1A requires two panel scenarios.")
  list(
    independent = plot_phase_spectrum(spikes, bulk, scenarios[1],
                                      "phase_panel_independent_spectrum", source, root),
    correlated = plot_phase_spectrum(spikes, bulk, scenarios[2],
                                     "phase_panel_correlated_spectrum", source, root)
  )
}

plot_phase_margin_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "experiment_1b_phase_margin"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  summary <- summarize_phase_margin(data)
  summary <- summary[summary$phase_region == "main", ]
  summary <- summary[order(summary$true_Delta), ]
  bias <- data.frame(
    phase_margin = rep(summary$true_Delta, 2L),
    method = rep(c("FPCA", "Proposed"), each = nrow(summary)),
    relative_bias_percent = 100 * c(summary$naive_relative_bias,
                                    summary$corrected_relative_bias)
  )
  bias_plot <- save_publication_figure(
    "phase_spike_bias", bias,
    function(x) {
      methods <- c("FPCA", "Proposed"); colors <- c("#D55E00", "#0072B2")
      ylim <- range(c(0, x$relative_bias_percent), finite = TRUE)
      plot(NA, xlim = range(x$phase_margin), ylim = ylim,
           xlab = "Phase margin", ylab = "Relative bias (%)")
      for (i in seq_along(methods)) {
        z <- x[x$method == methods[i], ]
        lines(z$phase_margin, z$relative_bias_percent, type = "b", pch = 15 + i,
              col = colors[i], lwd = 2)
      }
      abline(h = 0, lty = 3, col = "grey45")
      legend("topright", methods, col = colors, pch = 16:17, lty = 1, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "n=p=300; M=K0=1; main phase margins 0.4--0.9",
    quantity = "relative population-spike bias", comparison = "FPCA and Proposed"
  )
  variance <- data.frame(
    phase_margin = rep(summary$true_Delta, 2L),
    curve = rep(c("Empirical", "Theory"), each = nrow(summary)),
    scaled_variance = c(summary$scaled_alpha_variance,
                        summary$theoretical_scaled_variance)
  )
  variance_plot <- save_publication_figure(
    "phase_spike_variance", variance,
    function(x) {
      colors <- c("#0072B2", "#D55E00"); curves <- c("Empirical", "Theory")
      plot(NA, xlim = range(x$phase_margin), ylim = range(c(0, x$scaled_variance)),
           xlab = "Phase margin", ylab = expression(n*Var(hat(alpha))/alpha^2))
      for (i in seq_along(curves)) {
        z <- x[x$curve == curves[i], ]
        lines(z$phase_margin, z$scaled_variance, type = if (i == 1L) "b" else "l",
              pch = if (i == 1L) 16 else NA, lty = i, col = colors[i], lwd = 2)
      }
      legend("topright", c("Empirical", expression(2/Delta)), col = colors,
             pch = c(16, NA), lty = 1:2, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "n=p=300; M=K0=1; main phase margins 0.4--0.9",
    quantity = "scaled population-spike variance", comparison = "empirical and theoretical"
  )
  list(bias = bias_plot, variance = variance_plot)
}

plot_functional_inference_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "experiment_2_functional_inference"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  summary <- summarize_inference_performance(data)
  summary <- summary[summary$spike_index == 2L, ]
  plot_data <- data.frame(
    setting = paste(summary$scenario, summary$score_distribution, sep = " / "),
    alpha = summary$alpha_coverage_percent,
    Delta = summary$Delta_coverage_percent,
    r2 = summary$r2_coverage_percent
  )
  long <- reshape(plot_data, varying = c("alpha", "Delta", "r2"),
                  v.names = "coverage_percent", timevar = "quantity",
                  times = c("alpha", "Delta", "r2"), direction = "long")
  long$id <- NULL; rownames(long) <- NULL
  save_publication_figure(
    "functional_inference_coverage", long,
    function(x) {
      quantities <- c("alpha", "Delta", "r2"); colors <- c("#0072B2", "#D55E00", "#009E73")
      positions <- seq_len(length(unique(x$setting)))
      settings <- unique(x$setting)
      plot(NA, xlim = c(0.5, length(settings) + 0.5), ylim = range(c(90, x$coverage_percent)),
           xaxt = "n", xlab = "Covariance / score setting", ylab = "Coverage (%)")
      axis(1, positions, settings, las = 2, cex.axis = 0.75)
      for (i in seq_along(quantities)) {
        z <- x[x$quantity == quantities[i], ]; z <- z[match(settings, z$setting), ]
        lines(positions, z$coverage_percent, type = "b", pch = 14 + i, col = colors[i], lwd = 2)
      }
      abline(h = 95, lty = 3)
      legend("bottomleft", quantities, col = colors, pch = 15:17, lty = 1, bty = "n")
    }, root, width = 7.5, height = 5.8, source_runs = setNames(source, experiment),
    setting = "n=300; spike 2; covariance and score robustness grid",
    quantity = "Wald coverage", comparison = "alpha, Delta, and r2"
  )
}

qq_plot_data <- function(values) {
  values <- sort(values[is.finite(values)])
  data.frame(normal_quantile = qnorm(ppoints(length(values))), empirical_quantile = values)
}

plot_figure_3 <- function(root = PROJECT_ROOT, run_ids = list()) {
  experiment <- "experiment_2b_functional_asymptotics"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_ids[[experiment]]))
  source <- resolved_run_id(experiment, root, run_ids[[experiment]])
  x <- data[data$scenario == FIGURE3_SCENARIO &
              tolower(data$score_distribution) == tolower(FIGURE3_DISTRIBUTION) &
              data$spike_index == FIGURE3_SPIKE_INDEX & data$n == FIGURE3_N, ]
  z <- list(
    alpha = (x$hat_alpha - x$true_alpha) / x$se_alpha,
    delta = (x$hat_Delta - x$true_Delta) / x$se_Delta,
    r2 = (x$hat_r2 - x$true_r2) / x$se_r2
  )
  outputs <- list()
  for (quantity in names(z)) {
    outputs[[quantity]] <- save_publication_figure(
      paste0("qq_", quantity), qq_plot_data(z[[quantity]]),
      function(y) {
        plot(y$normal_quantile, y$empirical_quantile, pch = 16, cex = 0.48,
             col = grDevices::adjustcolor("#0072B2", 0.55),
             xlab = "Normal quantiles", ylab = "Empirical quantiles")
        abline(0, 1, col = "#D55E00", lwd = 2)
      }, root, source_runs = setNames(source, experiment),
      setting = "Haar covariance; standardized Uniform scores; n=300; spike 2",
      quantity = paste("studentized error of", quantity), comparison = "standard normal reference"
    )
  }
  joint <- data.frame(
    standardized_alpha = rep(z$alpha, 2L),
    quantity = rep(c("Delta", "r2"), each = length(z$alpha)),
    standardized_error = c(z$delta, z$r2)
  )
  outputs$joint <- save_publication_figure(
    "joint_standardized_errors", joint,
    function(y) {
      colors <- c("#0072B2", "#009E73"); quantities <- c("Delta", "r2")
      limits <- range(c(y$standardized_alpha, y$standardized_error), finite = TRUE)
      plot(NA, xlim = limits, ylim = limits,
           xlab = "Standardized error of alpha", ylab = "Standardized error")
      for (i in seq_along(quantities)) {
        a <- y[y$quantity == quantities[i], ]
        points(a$standardized_alpha, a$standardized_error, pch = 14 + i, cex = 0.45,
               col = grDevices::adjustcolor(colors[i], 0.35))
      }
      abline(0, 1, col = "#D55E00", lwd = 2)
      legend("topleft", quantities, col = colors, pch = 15:16, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "Haar covariance; standardized Uniform scores; n=300; spike 2",
    quantity = "joint studentized errors", comparison = "Delta and r2 against alpha"
  )
  invisible(outputs)
}

draw_method_means <- function(x, ylab) {
  methods <- c("FPCA", "Proposed"); colors <- c("#D55E00", "#0072B2")
  truth <- unique(x$truth)
  ylim <- range(c(x$mean, truth), finite = TRUE)
  span <- diff(ylim)
  lower_padding <- max(0.01, 0.08 * span)
  upper_padding <- max(0.01, 0.35 * span)
  plot(NA, xlim = range(x$n), ylim = ylim + c(-lower_padding, upper_padding),
       xlab = "n", ylab = ylab)
  for (i in seq_along(methods)) {
    z <- x[x$method == methods[i], ]; z <- z[order(z$n), ]
    lines(z$n, z$mean, type = "b", pch = 15 + i, col = colors[i], lwd = 2)
  }
  abline(h = truth[1], lty = 2, col = "grey30", lwd = 2)
  legend("topright", c(methods, "Population truth"),
         col = c(colors, "grey30"), pch = c(16, 17, NA), lty = c(1, 1, 2), bty = "n")
}

plot_fpca_comparison_results <- function(root = PROJECT_ROOT, run_ids = list()) {
  bundle <- generate_fpca_comparison_summaries(root, run_ids)
  d <- bundle$data; s <- bundle$source_run_ids
  list(
    large_eigenvalue = save_publication_figure(
      "large_domain_eigenvalue_comparison", d$large_domain_fpca_eigenvalue,
      function(x) draw_method_means(x, "Population-spike estimate"), root,
      source_runs = c(experiment_3a_simple_direction = s[["direction"]]),
      setting = "large-domain simple-spike experiment; n=150,300,600",
      quantity = "Monte Carlo mean population-spike estimate", comparison = "FPCA and Proposed"
    ),
    direction = save_publication_figure(
      "large_domain_direction_comparison", d$large_domain_direction_comparison,
      function(x) draw_method_means(x, "Direction functional estimate"), root,
      source_runs = c(experiment_3a_simple_direction = s[["direction"]]),
      setting = "large-domain simple-direction experiment; truth=0.6",
      quantity = "Monte Carlo mean direction functional", comparison = "FPCA and Proposed"
    ),
    eigenspace = save_publication_figure(
      "large_domain_eigenspace_comparison", d$large_domain_eigenspace_comparison,
      function(x) draw_method_means(x, "Eigenspace functional estimate"), root,
      source_runs = c(experiment_3b_repeated_eigenspace = s[["eigenspace"]]),
      setting = "large-domain repeated-eigenspace experiment; truth=0.45",
      quantity = "Monte Carlo mean projector functional", comparison = "FPCA and Proposed"
    ),
    alignment = save_publication_figure(
      "large_domain_individual_alignment", d$large_domain_individual_alignment,
      function(x) {
        plot(x$n, x$sd_squared_alignment, type = "b", pch = 16, col = "#0072B2", lwd = 2,
             xlab = "n", ylab = "SD of squared alignment")
      }, root, source_runs = c(experiment_3b_repeated_eigenspace = s[["eigenspace"]]),
      setting = "large-domain repeated-eigenspace experiment",
      quantity = "Monte Carlo SD of individual squared alignment", comparison = "none"
    ),
    functional_eigenvalue = save_publication_figure(
      "functional_eigenvalue_fpca_vs_proposed", d$functional_eigenvalue_fpca_vs_proposed,
      function(x) draw_method_means(x, "Population-spike estimate"), root,
      source_runs = c(experiment_2b_functional_asymptotics = s[["functional"]]),
      setting = "Haar covariance; Uniform scores; spike 2; n=150,300,600",
      quantity = "Monte Carlo mean population-spike estimate", comparison = "FPCA and Proposed"
    ),
    reliability = save_publication_figure(
      "functional_reliability_check", d$functional_reliability_check,
      function(x) {
        ylim <- range(c(x$mean_empirical_squared_overlap, x$mean_proposed_r2, x$true_r2))
        plot(x$n, x$mean_empirical_squared_overlap, type = "b", pch = 16,
             col = "#D55E00", lwd = 2, ylim = ylim,
             xlab = "n", ylab = "Squared alignment")
        lines(x$n, x$mean_proposed_r2, type = "b", pch = 17, col = "#0072B2", lwd = 2)
        lines(x$n, x$true_r2, lty = 2, col = "grey30", lwd = 2)
        legend("bottomright", c("Empirical overlap", "Proposed r2 estimate", "Population r2"),
               col = c("#D55E00", "#0072B2", "grey30"), pch = c(16, 17, NA),
               lty = c(1, 1, 2), bty = "n")
      }, root, source_runs = c(experiment_2b_functional_asymptotics = s[["functional"]]),
      setting = "Haar covariance; Uniform scores; spike 2; n=150,300,600",
      quantity = "squared population overlap and reliability", comparison = "empirical overlap and Proposed r2"
    )
  )
}

plot_simple_direction_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  ids <- list(experiment_3a_simple_direction = run_id)
  plot_fpca_comparison_results(root, ids)$direction
}

plot_repeated_eigenspace_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  ids <- list(experiment_3b_repeated_eigenspace = run_id)
  plot_fpca_comparison_results(root, ids)$eigenspace
}

plot_grid_saturation_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "experiment_3c_grid_saturation"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  summary <- summarize_grid(data)
  save_publication_figure(
    "large_domain_discretization_error", summary,
    function(x) {
      quantities <- c("alpha", "Delta", "r2", "theta")
      colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
      plot(NA, xlim = range(x$m_over_n), ylim = range(c(0, x$q90)),
           xlab = "m/n", ylab = expression(sqrt(n)*abs(hat(T)[grid]-hat(T)[oracle])))
      for (i in seq_along(quantities)) {
        z <- x[x$quantity == quantities[i], ]; z <- z[order(z$m_over_n), ]
        lines(z$m_over_n, z$median, type = "b", pch = 14 + i, col = colors[i], lwd = 2)
      }
      legend("topright", quantities, col = colors, pch = 15:18, lty = 1, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "n=300; dense-grid functional experiment",
    quantity = "median root-n discretization error", comparison = "grid and oracle implementations"
  )
}

plot_k0_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "robustness_k0"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  s <- summarize_k0(data)
  s <- s[s$score_distribution == "Uniform" & s$K0 == 5L & s$spike_index == 2L, ]
  plot_data <- data.frame(
    n = rep(s$n, 3L), quantity = rep(c("alpha", "Delta", "r2"), each = nrow(s)),
    n_paired_rmse = c(s$alpha_n_overdeflation_rmse, s$Delta_n_overdeflation_rmse,
                      s$r2_n_overdeflation_rmse)
  )
  save_publication_figure(
    "overdeflation_rmse", plot_data,
    function(x) {
      quantities <- c("alpha", "Delta", "r2"); colors <- c("#0072B2", "#D55E00", "#009E73")
      plot(NA, xlim = range(x$n), ylim = range(c(0, x$n_paired_rmse)),
           xlab = "n", ylab = expression(n %*% " paired RMSE"))
      for (i in seq_along(quantities)) {
        z <- x[x$quantity == quantities[i], ]
        lines(z$n, z$n_paired_rmse, type = "b", pch = 14 + i, col = colors[i], lwd = 2)
      }
      legend("topright", quantities, col = colors, pch = 15:17, lty = 1, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "Haar covariance; Uniform scores; spike 2; K0=5 versus baseline K0=M=3",
    quantity = "n times paired K0-sensitivity RMSE", comparison = "alpha, Delta, and r2"
  )
}

plot_universality_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  experiment <- "robustness_universality"
  data <- read_required_csv(resolve_run_artifact(experiment, "replicate", root, run_id))
  source <- resolved_run_id(experiment, root, run_id)
  s <- summarize_universality(data); s <- s[s$spike_index == 2L, ]
  rows <- list(); index <- 1L
  for (quantity in c("alpha", "Delta", "r2")) for (variance in c("Universal", "Full cumulant")) {
    column <- paste0(if (variance == "Universal") "universal" else "full", "_sd_ratio_", quantity)
    rows[[index]] <- data.frame(score_distribution = s$score_distribution,
                               quantity = quantity, variance = variance,
                               theoretical_sd_over_empirical_sd = s[[column]])
    index <- index + 1L
  }
  plot_data <- bind_rows_base(rows)
  save_publication_figure(
    "fourth_cumulant_variance", plot_data,
    function(x) {
      distributions <- unique(x$score_distribution); positions <- seq_along(distributions)
      quantities <- c("alpha", "Delta", "r2"); colors <- c("#0072B2", "#D55E00", "#009E73")
      plot(NA, xlim = c(0.5, length(distributions) + 0.5),
           ylim = range(c(1, x$theoretical_sd_over_empirical_sd)), xaxt = "n",
           xlab = "Score distribution", ylab = "Theoretical SD / empirical SD")
      axis(1, positions, distributions)
      for (i in seq_along(quantities)) for (variance in c("Universal", "Full cumulant")) {
        z <- x[x$quantity == quantities[i] & x$variance == variance, ]
        z <- z[match(distributions, z$score_distribution), ]
        lines(positions, z$theoretical_sd_over_empirical_sd, type = "b",
              pch = if (variance == "Universal") 16 else 17,
              lty = if (variance == "Universal") 1 else 2, col = colors[i], lwd = 2)
      }
      abline(h = 1, lty = 3)
      legend("topleft", quantities, col = colors, lty = 1, bty = "n")
      legend("bottomright", c("Universal", "Full cumulant"), col = "grey25",
             pch = c(16, 17), lty = 1:2, bty = "n")
    }, root, source_runs = setNames(source, experiment),
    setting = "diagonal functional covariance; spike 2; Gaussian, t12, and Uniform scores",
    quantity = "theoretical-to-empirical SD ratio", comparison = "universal and full fourth-cumulant variance"
  )
}

plot_functional_asymptotic_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  ids <- list(experiment_2b_functional_asymptotics = run_id)
  list(comparison = plot_fpca_comparison_results(root, ids)$functional_eigenvalue,
       reliability = plot_fpca_comparison_results(root, ids)$reliability)
}

plot_functional_asymptotic_bundle <- function(root = PROJECT_ROOT, run_id = NULL) {
  list(asymptotic = plot_functional_asymptotic_results(root, run_id),
       joint = plot_figure_3(root, list(experiment_2b_functional_asymptotics = run_id)))
}

plot_figure_2 <- function(root = PROJECT_ROOT, run_ids = list()) {
  list(panel = plot_panel_phase_results(root, run_ids[["experiment_1a_panel_phase"]]),
       margin = plot_phase_margin_results(root, run_ids[["experiment_1b_phase_margin"]]))
}

plot_figure_4 <- function(root = PROJECT_ROOT, run_ids = list()) {
  list(comparisons = plot_fpca_comparison_results(root, run_ids),
       discretization = plot_grid_saturation_results(root, run_ids[["experiment_3c_grid_saturation"]]))
}

plot_robustness_results <- function(root = PROJECT_ROOT, run_ids = list()) {
  list(overdeflation = plot_k0_results(root, run_ids[["robustness_k0"]]),
       cumulant = plot_universality_results(root, run_ids[["robustness_universality"]]))
}

plot_experiment_results <- function(experiment, root = PROJECT_ROOT, run_id = NULL) {
  switch(
    experiment,
    experiment_1a_panel_phase = plot_panel_phase_results(root, run_id),
    experiment_1b_phase_margin = plot_phase_margin_results(root, run_id),
    experiment_2_functional_inference = plot_functional_inference_results(root, run_id),
    experiment_2b_functional_asymptotics = plot_functional_asymptotic_bundle(root, run_id),
    experiment_3a_simple_direction = plot_simple_direction_results(root, run_id),
    experiment_3b_repeated_eigenspace = plot_repeated_eigenspace_results(root, run_id),
    experiment_3c_grid_saturation = plot_grid_saturation_results(root, run_id),
    robustness_k0 = plot_k0_results(root, run_id),
    robustness_universality = plot_universality_results(root, run_id),
    stop("No plot specification for experiment: ", experiment)
  )
}

plot_selected_results <- function(experiments, root = PROJECT_ROOT, run_ids = list(),
                                  include_comparison_exports = TRUE) {
  for (experiment in experiments) plot_experiment_results(experiment, root, run_ids[[experiment]])
  # This call creates only additional single-panel FPCA comparison exports.
  if (include_comparison_exports &&
      all(c("experiment_2b_functional_asymptotics", "experiment_3a_simple_direction",
            "experiment_3b_repeated_eigenspace") %in% experiments)) {
    plot_fpca_comparison_results(root, run_ids)
  }
  invisible(NULL)
}

plot_all_saved_results <- function(root = PROJECT_ROOT, run_ids = list()) {
  experiments <- names(summary_specifications())
  generate_all_summaries(root, run_ids, experiments)
  generate_fpca_comparison_summaries(root, run_ids)
  plot_selected_results(experiments, root, run_ids, include_comparison_exports = FALSE)
  invisible(NULL)
}
