# Shared infrastructure for all real-data applications.
#
# Rows of every input matrix are independent observational units. Columns are
# weighted discretized Hilbert coordinates. All real-data
# fits are centered across observations and retain the covariance divisor n.

required_real_data_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      ". Install them with Rscript code/install_packages.R."
    )
  }
  invisible(TRUE)
}

sha256_file <- function(path) {
  if (!nzchar(Sys.which("sha256sum"))) {
    stop("The sha256sum command is required to verify downloaded files.")
  }
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(output)) stop("Could not compute SHA256 for ", path)
  strsplit(output[1], "[[:space:]]+")[[1]][1]
}

real_data_paths <- function(root = PROJECT_ROOT, dataset = NULL) {
  base <- list(
    data = file.path(root, "data"),
    output = file.path(root, "output", "real_data"),
    manuscript = file.path(root, "manuscript", "real_data"),
    manifest = file.path(root, "output", "real_data", "manifest")
  )
  if (!is.null(dataset)) {
    base$data <- file.path(base$data, dataset)
    base$output <- file.path(base$output, dataset)
    base$manuscript <- file.path(base$manuscript, dataset)
  }
  invisible(lapply(base, dir.create, recursive = TRUE, showWarnings = FALSE))
  base
}

write_json_atomic <- function(x, path) {
  required_real_data_packages("jsonlite")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
  invisible(path)
}

register_real_data_run <- function(dataset, run_id, artifacts, root = PROJECT_ROOT) {
  paths <- real_data_paths(root)
  row <- data.frame(
    dataset = dataset,
    run_id = run_id,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  for (name in names(artifacts)) row[[name]] <- relative_to_root(artifacts[[name]], root)
  registry_path <- file.path(paths$manifest, "run_registry.csv")
  registry <- if (file.exists(registry_path)) {
    read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }
  if (nrow(registry)) {
    registry <- registry[!(registry$dataset == dataset & registry$run_id == run_id), , drop = FALSE]
    registry <- bind_rows_base(list(registry, row))
  } else {
    registry <- row
  }
  write_csv_atomic(registry, registry_path)
  write_csv_atomic(row, file.path(paths$manifest, paste0(dataset, "__latest.csv")))
  invisible(row)
}

resolve_real_data_run <- function(dataset, artifact, root = PROJECT_ROOT, run_id = NULL) {
  paths <- real_data_paths(root)
  manifest_path <- if (is.null(run_id)) {
    file.path(paths$manifest, paste0(dataset, "__latest.csv"))
  } else {
    file.path(paths$manifest, "run_registry.csv")
  }
  manifest <- read_required_csv(manifest_path)
  if (!is.null(run_id)) {
    manifest <- manifest[manifest$dataset == dataset & manifest$run_id == run_id, , drop = FALSE]
  }
  if (nrow(manifest) != 1L) stop("Could not uniquely resolve real-data run for ", dataset, ".")
  if (!artifact %in% names(manifest) || is.na(manifest[[artifact]]) || !nzchar(manifest[[artifact]])) {
    stop("Artifact '", artifact, "' is not registered for real-data run ", manifest$run_id)
  }
  path <- manifest[[artifact]]
  if (!grepl("^/", path)) path <- file.path(root, path)
  if (!file.exists(path)) stop("Registered real-data artifact is missing: ", path)
  path
}

resolved_real_data_run_id <- function(dataset, root = PROJECT_ROOT, run_id = NULL) {
  paths <- real_data_paths(root)
  manifest_path <- if (is.null(run_id)) {
    file.path(paths$manifest, paste0(dataset, "__latest.csv"))
  } else {
    file.path(paths$manifest, "run_registry.csv")
  }
  manifest <- read_required_csv(manifest_path)
  if (!is.null(run_id)) manifest <- manifest[manifest$dataset == dataset & manifest$run_id == run_id, , drop = FALSE]
  if (nrow(manifest) != 1L) stop("Could not uniquely resolve real-data run ID for ", dataset, ".")
  manifest$run_id
}

normalize_direction <- function(direction) {
  direction <- as.numeric(direction)
  norm_value <- sqrt(sum(direction^2))
  if (!is.finite(norm_value) || norm_value <= 0) stop("A scientific direction must have positive finite norm.")
  direction / norm_value
}

spectral_effective_rank <- function(eigenvalues) {
  eigenvalues <- as.numeric(eigenvalues)
  if (!length(eigenvalues) || !is.finite(eigenvalues[1]) || eigenvalues[1] <= 0) return(NA_real_)
  sum(eigenvalues) / eigenvalues[1]
}

component_separation_diagnostics <- function(eigenvalues, indices) {
  values <- sort(as.numeric(eigenvalues), decreasing = TRUE)
  rows <- lapply(indices, function(j) {
    previous_gap <- if (j == 1L) Inf else values[j - 1L] - values[j]
    next_gap <- if (j == length(values)) Inf else values[j] - values[j + 1L]
    # A shared adjacent pair must receive the same relative-gap value from
    # both sides. Each gap is therefore scaled by the larger eigenvalue in
    # that pair, rather than by the eigenvalue of the row being evaluated.
    previous_relative_gap <- if (j == 1L) Inf else previous_gap / values[j - 1L]
    next_relative_gap <- if (j == length(values)) Inf else next_gap / values[j]
    data.frame(
      spike_index = j,
      sample_gap_previous = previous_gap,
      sample_gap_next = next_gap,
      minimum_neighbor_gap = min(previous_gap, next_gap),
      previous_relative_gap = previous_relative_gap,
      next_relative_gap = next_relative_gap,
      relative_neighbor_gap = min(previous_relative_gap, next_relative_gap)
    )
  })
  bind_rows_base(rows)
}

# Legacy descriptive gap diagnostic retained for backward compatibility only.
# These 0.05-threshold fields are not inferential validity gates and are omitted
# from the compact manuscript-facing PC1--PC6 output.
refresh_component_separation_flags <- function(table, minimum_relative_gap = 0.05) {
  required <- c("hat_lambda", "sample_gap_previous", "sample_gap_next", "phase_lower_bound")
  missing <- setdiff(required, names(table))
  if (length(missing)) stop("Cannot refresh separation flags; missing columns: ", paste(missing, collapse = ", "))
  table$previous_relative_gap <- ifelse(
    is.infinite(table$sample_gap_previous), Inf,
    table$sample_gap_previous / (table$hat_lambda + table$sample_gap_previous)
  )
  table$next_relative_gap <- ifelse(
    is.infinite(table$sample_gap_next), Inf,
    table$sample_gap_next / table$hat_lambda
  )
  table$minimum_neighbor_gap <- pmin(table$sample_gap_previous, table$sample_gap_next)
  table$relative_neighbor_gap <- pmin(
    table$previous_relative_gap, table$next_relative_gap
  )
  table$phase_interior <- table$phase_lower_bound > 0
  table$separation_pass <- table$relative_neighbor_gap >= minimum_relative_gap
  table$regular_inference_valid <- table$phase_interior & table$separation_pass
  table$inference_status <- ifelse(
    table$regular_inference_valid, "regular_separated",
    ifelse(!table$phase_interior, "phase_not_certified", "weak_sample_separation")
  )
  table
}

infer_real_components <- function(fit, K0, J, confidence_level = 0.95,
                                  certification_threshold = 0.20,
                                  minimum_relative_gap = 0.05,
                                  directions = list(), orientation_direction = NULL,
                                  equality_draws = 100000L,
                                  equality_seed = 82471L) {
  n <- fit$n
  J <- min(as.integer(J), as.integer(K0), length(fit$values) - 1L)
  if (J < 1L) stop("J must leave at least one empirical bulk eigenvalue.")
  rows <- lapply(seq_len(J), function(j) {
    proposed <- dimension_free_inference(
      fit$values, j, K0, n,
      confidence_level = confidence_level,
      certification_threshold = certification_threshold
    )
    cbind(proposed, fpca_eigenvalue_interval(fit, j, confidence_level))
  })
  inference <- bind_rows_base(rows)
  diagnostics <- component_separation_diagnostics(fit$values, seq_len(J))
  inference <- merge(inference, diagnostics, by = "spike_index", sort = FALSE)
  inference <- inference[order(inference$spike_index), , drop = FALSE]
  inference <- refresh_component_separation_flags(inference, minimum_relative_gap)

  gap_table <- adjacent_gap_inference(
    fit, K0, J, confidence_level = confidence_level,
    equality_draws = equality_draws, equality_seed = equality_seed
  )
  inference$phase_valid <- inference$phase_interior
  inference$fpca_simple_gap_wald_valid <- TRUE
  inference$fpca_equality_rejected <- TRUE
  inference$proposed_equality_rejected <- TRUE
  inference$proposed_gap_local_inference_valid <- TRUE
  inference$prespecified_cluster_member <- FALSE
  for (j in seq_len(J)) {
    adjacent <- gap_table[gap_table$upper_rank %in% c(j - 1L, j), , drop = FALSE]
    if (nrow(adjacent)) {
      inference$fpca_simple_gap_wald_valid[j] <- all(adjacent$fpca_gap_wald_regular)
      inference$fpca_equality_rejected[j] <- all(adjacent$fpca_equality_rejected)
      inference$proposed_equality_rejected[j] <- all(adjacent$proposed_equality_rejected)
      inference$proposed_gap_local_inference_valid[j] <-
        all(adjacent$proposed_gap_local_inference_valid)
    }
  }
  feature_vectors <- lapply(seq_len(J), function(j) feature_eigenvector(fit, j))
  if (!is.null(orientation_direction)) {
    orientation_direction <- normalize_direction(orientation_direction)
    for (j in seq_len(J)) {
      if (sum(feature_vectors[[j]] * orientation_direction) < 0) {
        feature_vectors[[j]] <- -feature_vectors[[j]]
      }
    }
  }
  if (length(directions)) {
    directions <- lapply(directions, normalize_direction)
    for (name in names(directions)) {
      raw <- vapply(feature_vectors, function(v) sum(v * directions[[name]]), numeric(1))
      inference[[paste0("raw_", name)]] <- raw
      inference[[paste0("corrected_", name)]] <- ifelse(
        inference$hat_r2 > 0, raw / sqrt(inference$hat_r2), NA_real_
      )
    }
  }
  list(table = inference, vectors = feature_vectors, gaps = gap_table)
}

adjacent_gap_inference <- function(fit, K0, J, confidence_level = 0.95,
                                   equality_draws = 100000L,
                                   equality_seed = 82471L,
                                   known_distinct = NULL) {
  J <- min(as.integer(J), as.integer(K0), length(fit$values) - 1L)
  alpha <- 1 - confidence_level
  bind_rows_base(lapply(seq_len(J - 1L), function(j) {
    gap <- multiplicity2_gap_inference(fit$values, j, fit$n, K0, confidence_level)
    classical <- fpca_pair_inference(
      fit, j, confidence_level, equality_draws, equality_seed + j,
      known_distinct = known_distinct
    )
    pooled <- pooled_cluster_inference(fit$values, j:(j + 1L), K0, fit$n, confidence_level)
    gap <- cbind(classical, gap)
    gap$pair <- paste0("PC", j, "-PC", j + 1L)
    gap$upper_rank <- j
    gap$fpca_p_equal <- gap$fpca_p_equal_general
    gap$proposed_p_equal <- gap$p_proposed
    gap$sample_HC_p_equal <- gap$p_sample_HC
    gap$fpca_reject_equal <- !is.na(gap$fpca_p_equal_general) && gap$fpca_p_equal_general < alpha
    gap$proposed_reject_equal <- !is.na(gap$p_proposed) && gap$p_proposed < alpha
    gap$phase_separation_certified <- pooled$pooled_Delta_one_sided_lower > 0
    gap$gap_inference_valid <- gap$inference_valid &
      pooled$pooled_inference_valid & gap$phase_separation_certified
    gap$phase_valid <- gap$phase_separation_certified
    gap$fpca_simple_gap_wald_valid <- gap$fpca_gap_wald_regular
    gap$fpca_equality_rejected <- gap$fpca_reject_equal
    gap$proposed_equality_rejected <- gap$proposed_reject_equal
    gap$prespecified_cluster_member <- FALSE
    gap$equality_evidence <- ifelse(
      gap$fpca_reject_equal & gap$proposed_reject_equal, "both_reject_equality",
      ifelse(!gap$fpca_reject_equal & !gap$proposed_reject_equal,
             "both_insufficient_evidence", "methods_differ")
    )
    cbind(gap, pooled[, setdiff(names(pooled), c("lambda_pool")), drop = FALSE])
  }))
}

compact_real_component_inference <- function(table) {
  output <- table[, c(
    "spike_index", "hat_lambda", "fpca_alpha_se", "fpca_alpha_ci_lower",
    "fpca_alpha_ci_upper", "hat_alpha", "se_alpha", "ci_alpha_lower",
    "ci_alpha_upper", "hat_Delta", "se_Delta", "ci_Delta_lower",
    "ci_Delta_upper", "phase_lower_bound", "hat_r2", "se_r2",
    "ci_r2_lower", "ci_r2_upper", "phase_valid",
    "fpca_simple_gap_wald_valid", "fpca_equality_rejected",
    "proposed_equality_rejected", "proposed_gap_local_inference_valid",
    "prespecified_cluster_member"
  ), drop = FALSE]
  names(output) <- c(
    "rank", "sample_eigenvalue", "fpca_eigenvalue_se",
    "fpca_eigenvalue_ci_lower", "fpca_eigenvalue_ci_upper",
    "proposed_alpha", "proposed_alpha_se", "proposed_alpha_ci_lower",
    "proposed_alpha_ci_upper", "hat_Delta", "Delta_se", "Delta_ci_lower",
    "Delta_ci_upper", "Delta_one_sided_lower", "hat_r2", "r2_se",
    "r2_ci_lower", "r2_ci_upper", "phase_valid",
    "fpca_gap_wald_regular", "fpca_equality_rejected",
    "proposed_equality_rejected", "proposed_gap_local_inference_valid",
    "prespecified_cluster_member"
  )
  output
}

compact_real_gap_inference <- function(gaps) {
  output <- gaps[, c(
    "pair", "raw_relative_gap_legacy", "sample_gap_sym", "fpca_gap_sym",
    "fpca_gap_wald_se", "fpca_gap_wald_ci_raw_lower",
    "fpca_gap_wald_ci_raw_upper", "fpca_gap_wald_ci_lower",
    "fpca_gap_wald_ci_upper", "fpca_gap_wald_regular",
    "fpca_p_equal_general", "fpca_anderson_Q", "fpca_anderson_p_equal",
    "corrected_relative_gap_legacy", "proposed_gap_sym", "T_proposed",
    "proposed_p_equal", "proposed_gap_ci_raw_lower",
    "proposed_gap_ci_raw_upper", "proposed_gap_ci_raw_empty",
    "proposed_gap_upper95", "proposed_gap_upper95_boundary_fallback",
    "proposed_gap_ci_boundary_lower",
    "proposed_gap_ci_boundary_upper", "Delta_pool", "r2_pool",
    "T_sample_HC", "p_sample_HC", "phase_valid",
    "fpca_equality_rejected",
    "proposed_equality_rejected", "proposed_gap_local_inference_valid",
    "prespecified_cluster_member", "equality_evidence"
  ), drop = FALSE]
  output$fpca_gap_ci_status <- ifelse(
    output$fpca_gap_wald_regular,
    "regular distinct-root interval",
    "not valid under near-tie"
  )
  nonregular <- !output$fpca_gap_wald_regular
  output$fpca_gap_wald_ci_lower[nonregular] <- NA_real_
  output$fpca_gap_wald_ci_upper[nonregular] <- NA_real_
  output
}

gap_k0_sensitivity <- function(fit, K0_values, J, confidence_level = 0.95,
                               equality_draws = 100000L,
                               equality_seed = 82471L) {
  bind_rows_base(lapply(sort(unique(as.integer(K0_values))), function(K0) {
    result <- adjacent_gap_inference(
      fit, K0, min(J, K0), confidence_level,
      equality_draws = equality_draws, equality_seed = equality_seed
    )
    result$K0_sensitivity <- K0
    result
  }))
}

prespecified_cluster_inference <- function(gaps, candidate_pairs) {
  candidate_pairs <- unique(as.character(candidate_pairs))
  parsed <- lapply(candidate_pairs, function(pair) {
    as.integer(strsplit(sub("^PC", "", pair), "-PC", fixed = TRUE)[[1]])
  })
  members <- unlist(parsed, use.names = FALSE)
  if (anyDuplicated(members)) stop("Prespecified candidate clusters must not overlap.")
  keep <- gaps$pair %in% candidate_pairs
  columns <- c(
    "pair", "multiplicity", "lambda_pool", "pooled_alpha", "pooled_alpha_se",
    "pooled_alpha_ci_lower", "pooled_alpha_ci_upper", "pooled_Delta",
    "pooled_Delta_se", "pooled_Delta_ci_lower", "pooled_Delta_ci_upper",
    "pooled_Delta_one_sided_lower", "pooled_r2", "pooled_r2_se",
    "pooled_r2_ci_lower", "pooled_r2_ci_upper", "pooled_inference_valid"
  )
  output <- gaps[keep, columns, drop = FALSE]
  output$candidate_status <- "prespecified_candidate_cluster"
  output
}

k0_sensitivity <- function(fit, K0_values, J, confidence_level = 0.95,
                           certification_threshold = 0.20,
                           minimum_relative_gap = 0.05) {
  rows <- lapply(sort(unique(as.integer(K0_values))), function(K0) {
    if (K0 >= length(fit$values)) return(NULL)
    result <- infer_real_components(
      fit, K0, min(J, K0), confidence_level, certification_threshold,
      minimum_relative_gap
    )$table
    result$K0_sensitivity <- K0
    result
  })
  bind_rows_base(rows)
}

write_real_inference_table <- function(data, path, dataset_label) {
  direction_columns <- grep("^corrected_", names(data), value = TRUE)
  columns <- c(
    "spike_index", "hat_lambda", "hat_alpha", "hat_Delta", "phase_lower_bound",
    "hat_r2", direction_columns, "relative_neighbor_gap", "phase_valid",
    "prespecified_cluster_member"
  )
  display <- data[, columns, drop = FALSE]
  direction_labels <- gsub("_", "-", sub("^corrected_", "Proposed ", direction_columns), fixed = TRUE)
  names(display) <- c(
    "PC", "$\\hat\\lambda$", "$\\hat\\alpha$", "$\\hat\\Delta$",
    "Phase lower", "$\\hat r^2$", direction_labels, "Relative gap", "Phase valid",
    "Candidate cluster"
  )
  display$PC <- as.character(as.integer(display$PC))
  write_latex_table(
    display, path, digits = 3L,
    caption = paste(dataset_label, "dimension-free PCA inference"),
    label = paste0("tab:", gsub("[^a-z0-9]+", "-", tolower(dataset_label)), "-inference")
  )
}

refresh_saved_separation_artifact <- function(path, minimum_relative_gap = 0.05) {
  data <- read_required_csv(path)
  data <- refresh_component_separation_flags(data, minimum_relative_gap)
  write_csv_atomic(data, path)
  invisible(data)
}

augment_real_data_run_artifacts <- function(dataset, run_id, artifacts,
                                            root = PROJECT_ROOT) {
  paths <- real_data_paths(root)
  registry_path <- file.path(paths$manifest, "run_registry.csv")
  registry <- read_required_csv(registry_path)
  row_index <- which(registry$dataset == dataset & registry$run_id == run_id)
  if (length(row_index) != 1L) stop("Could not uniquely augment real-data run ", run_id)
  for (name in names(artifacts)) {
    if (!name %in% names(registry)) registry[[name]] <- NA_character_
    registry[[name]][row_index] <- relative_to_root(artifacts[[name]], root)
  }
  write_csv_atomic(registry, registry_path)
  latest_path <- file.path(paths$manifest, paste0(dataset, "__latest.csv"))
  latest <- read_required_csv(latest_path)
  if (nrow(latest) == 1L && identical(latest$run_id, run_id)) {
    for (name in names(artifacts)) {
      latest[[name]] <- relative_to_root(artifacts[[name]], root)
    }
    write_csv_atomic(latest, latest_path)
  }
  invisible(artifacts)
}

# Use one compact manuscript canvas per statistical object.
real_data_publication_par <- function(mar = c(5.1, 5.3, 1.0, 1.0)) {
  par(
    mar = mar, cex = 1.08, cex.axis = 1.14, cex.lab = 1.25,
    mgp = c(3.2, 0.9, 0), tcl = -0.28,
    family = "sans"
  )
}

append_real_figure_manifest <- function(row, root = PROJECT_ROOT) {
  path <- file.path(real_data_paths(root)$output, "figure_output_manifest.csv")
  manifest <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE) else data.frame()
  if (nrow(manifest)) manifest <- manifest[manifest$filename != row$filename, , drop = FALSE]
  manifest <- if (nrow(manifest)) bind_rows_base(list(manifest, row)) else row
  write_csv_atomic(manifest, path)
  refresh_figure_output_summary(root)
  invisible(path)
}

save_real_publication_figure <- function(stem, plot_data, plotter, root, dataset,
                                        width = 6.4, height = 5.3, png_dpi = 200L,
                                        setting = "", quantity = "", comparison = "",
                                        source_run_id = "") {
  paths <- real_data_paths(root, dataset)
  data_path <- file.path(paths$output, paste0(stem, ".csv"))
  write_csv_atomic(plot_data, data_path)
  pdf_path <- file.path(paths$manuscript, paste0(stem, ".pdf"))
  png_path <- file.path(paths$manuscript, paste0(stem, ".png"))
  render <- function() {
    old <- real_data_publication_par(); on.exit(par(old), add = TRUE)
    plotter(read_required_csv(data_path))
  }
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  tryCatch(render(), finally = grDevices::dev.off())
  grDevices::png(png_path, width = width, height = height, units = "in", res = png_dpi)
  tryCatch(render(), finally = grDevices::dev.off())
  append_real_figure_manifest(data.frame(
    dataset = dataset, filename = relative_to_root(pdf_path, root),
    plot_data_file = relative_to_root(data_path, root), setting = setting,
    quantity = quantity, comparison = comparison,
    method_scope = figure_method_scope(quantity, comparison), source_run_id = source_run_id,
    stringsAsFactors = FALSE
  ), root)
  invisible(c(pdf = pdf_path, png = png_path, data = data_path))
}
