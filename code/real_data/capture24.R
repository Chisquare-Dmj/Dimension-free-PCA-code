# CAPTURE-24 free-living wrist-accelerometry large-domain application.
#
# Each independent observation is one participant's wake-aligned 24-hour
# activity function. Raw triaxial acceleration is converted to orientation-
# invariant ENMO in mg, averaged over non-overlapping 5-second epochs, and then
# transformed as log(1 + ENMO). Participants are not individually centered or
# scaled. PCA centers across participants and uses covariance divisor n.

CAPTURE24_RECORD_URL <- "https://ora.ox.ac.uk/objects/uuid:99d7c092-d865-4a19-b096-cc16440cd001"
CAPTURE24_DOWNLOAD_URL <- paste0(CAPTURE24_RECORD_URL, "/files/rpr76f381b")
CAPTURE24_LICENSE <- "Creative Commons Attribution 4.0 International (CC BY 4.0)"
CAPTURE24_ARCHIVE_BYTES <- 6902652480

default_capture24_config <- function() {
  list(
    # Nested domains measure how spectral complexity grows as longer prefixes
    # of the circularly rotated, wake-aligned daily function are included.
    domain_hours = c(6, 12, 18, 24),
    epoch_seconds = 5L,
    expected_sample_rate_hz = 100L,
    expected_participants = 151L,
    preprocessing_version = 2L,
    sleep_gap_tolerance_minutes = 20,
    minimum_main_sleep_hours = 1,
    preprocessing_cores = 10L,
    preprocessing_chunk_rows = 500000L,
    coverage_sensitivity_minutes = 30,
    # Reviewer-facing wake-anchor sensitivity. This changes only cohort
    # inclusion and reuses the deterministic five-second preprocessing cache.
    wake_anchor_sensitivity_hours = 3,
    K0 = 8L,
    K0_values = c(6L, 8L, 10L),
    J = 6L,
    confidence_level = 0.95,
    certification_threshold = 0.20,
    minimum_relative_gap = 0.05,
    eigengap_definition = "symmetric-v1",
    fpca_equality_draws = 100000L,
    fpca_equality_seed = 82471L,
    master_seed = MASTER_SEED,
    ncores = 10L,
    root = PROJECT_ROOT,
    use_cache = TRUE
  )
}

validate_capture24_config <- function(config) {
  if (!identical(as.numeric(config$domain_hours), c(6, 12, 18, 24))) {
    stop("CAPTURE-24 domain_hours must remain c(6, 12, 18, 24) for the prespecified pilot.")
  }
  if (config$epoch_seconds != 5L || 86400 %% config$epoch_seconds != 0L) {
    stop("CAPTURE-24 uses prespecified non-overlapping 5-second ENMO epochs.")
  }
  if (config$expected_participants != 151L) {
    stop("CAPTURE-24 expected_participants must remain 151 for the public QC cohort.")
  }
  if (config$preprocessing_version != 2L) {
    stop("CAPTURE-24 preprocessing_version must match the active timestamp-aware implementation.")
  }
  if (config$minimum_main_sleep_hours < 1) {
    stop("minimum_main_sleep_hours must be at least one hour for a stable wake anchor.")
  }
  if (config$wake_anchor_sensitivity_hours < config$minimum_main_sleep_hours) {
    stop("wake_anchor_sensitivity_hours cannot be below minimum_main_sleep_hours.")
  }
  if (config$J > config$K0 || !config$K0 %in% config$K0_values) {
    stop("J cannot exceed K0, and K0_values must contain the primary K0.")
  }
  if (!identical(config$eigengap_definition, "symmetric-v1") ||
      config$fpca_equality_draws < 100000L) {
    stop("CAPTURE-24 requires symmetric eigengaps and at least 100000 FPCA null draws.")
  }
  if (config$preprocessing_cores < 1L || config$preprocessing_chunk_rows < 500L) {
    stop("preprocessing_cores and preprocessing_chunk_rows must be positive.")
  }
  invisible(config)
}

capture24_source_dir <- function(config) {
  file.path(real_data_paths(config$root, "capture24")$data, "source")
}

capture24_archive_path <- function(config) {
  file.path(capture24_source_dir(config), "capture24.zip")
}

download_capture24 <- function(config = list()) {
  cfg <- validate_capture24_config(merge_config(default_capture24_config(), config))
  source_dir <- capture24_source_dir(cfg)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- capture24_archive_path(cfg)
  complete <- file.exists(destination) && identical(as.numeric(file.info(destination)$size),
                                                    as.numeric(CAPTURE24_ARCHIVE_BYTES))
  if (!complete) {
    if (file.exists(destination)) unlink(destination)
    temporary <- paste0(destination, ".part")
    if (file.exists(temporary)) unlink(temporary)
    pipeline_log(sprintf("CAPTURE-24 | DOWNLOAD START | expected_bytes=%s", CAPTURE24_ARCHIVE_BYTES))
    status <- system2(
      "curl",
      c("-fL", "--retry", "5", "--retry-delay", "5",
        "--connect-timeout", "30", "-o", shQuote(temporary), shQuote(CAPTURE24_DOWNLOAD_URL))
    )
    if (!identical(status, 0L) || !file.exists(temporary)) {
      stop("CAPTURE-24 archive download failed: ", CAPTURE24_DOWNLOAD_URL)
    }
    observed_bytes <- as.numeric(file.info(temporary)$size)
    if (!identical(observed_bytes, as.numeric(CAPTURE24_ARCHIVE_BYTES))) {
      stop("CAPTURE-24 archive has ", observed_bytes, " bytes; expected ", CAPTURE24_ARCHIVE_BYTES, ".")
    }
    if (!file.rename(temporary, destination)) stop("Could not finalize ", destination)
  }
  test_status <- system2("unzip", c("-tqq", shQuote(destination)), stdout = FALSE, stderr = FALSE)
  if (!identical(test_status, 0L)) stop("CAPTURE-24 ZIP integrity test failed.")
  manifest <- data.frame(
    source = CAPTURE24_RECORD_URL,
    download_url = CAPTURE24_DOWNLOAD_URL,
    license = CAPTURE24_LICENSE,
    file = basename(destination),
    bytes = as.numeric(file.info(destination)$size),
    sha256 = sha256_file(destination),
    verified_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(manifest, file.path(source_dir, "download_manifest.csv"))
  message("CAPTURE-24 download complete and ZIP integrity verified: ", destination)
  invisible(destination)
}

capture24_preprocessed_dir <- function(config) {
  key <- paste0("cfg", config_fingerprint(list(
    preprocessing_version = config$preprocessing_version,
    epoch_seconds = config$epoch_seconds,
    sleep_gap_tolerance_minutes = config$sleep_gap_tolerance_minutes
  )))
  file.path(real_data_paths(config$root, "capture24")$data, "processed", key)
}

read_capture24_binary <- function(path, n, m, size, what) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  values <- readBin(connection, what = what, n = n * m, size = size, endian = "little")
  if (length(values) != n * m) stop("Unexpected binary matrix size in ", path)
  matrix(values, nrow = n, ncol = m, byrow = TRUE)
}

prepare_capture24 <- function(config = list()) {
  cfg <- validate_capture24_config(merge_config(default_capture24_config(), config))
  required_real_data_packages("jsonlite")
  archive <- capture24_archive_path(cfg)
  if (!file.exists(archive) || file.info(archive)$size != CAPTURE24_ARCHIVE_BYTES) {
    stop("The complete CAPTURE-24 archive is missing. Run dataset=capture24 action=download first.")
  }
  processed_dir <- capture24_preprocessed_dir(cfg)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(processed_dir, "preprocessing_manifest.json")
  qc_path <- file.path(processed_dir, "subject_qc.csv")
  curve_path <- file.path(processed_dir, "wake_aligned_log1p_enmo_mg_float32.bin")
  coverage_path <- file.path(processed_dir, "wake_aligned_coverage_uint8.bin")
  required <- c(manifest_path, qc_path, curve_path, coverage_path)
  if (!cfg$use_cache || !all(file.exists(required))) {
    if (!nzchar(Sys.which("python3"))) stop("python3 is required for CAPTURE-24 streaming preprocessing.")
    python_packages <- system2(
      "python3", c("-c", shQuote("import numpy, pandas")),
      stdout = FALSE, stderr = FALSE
    )
    if (!identical(python_packages, 0L)) {
      stop("CAPTURE-24 preprocessing requires the Python packages numpy and pandas.")
    }
    if (dir.exists(processed_dir)) {
      unlink(list.files(processed_dir, full.names = TRUE), recursive = TRUE)
    }
    helper <- file.path(cfg$root, "code", "real_data", "capture24_preprocess.py")
    status <- system2(
      "python3",
      c(
        shQuote(helper), "--archive", shQuote(archive),
        "--output-dir", shQuote(processed_dir),
        "--epoch-seconds", cfg$epoch_seconds,
        "--sleep-gap-minutes", cfg$sleep_gap_tolerance_minutes,
        "--chunk-rows", cfg$preprocessing_chunk_rows,
        "--n-jobs", min(cfg$ncores, cfg$preprocessing_cores)
      )
    )
    if (!identical(status, 0L) || !all(file.exists(required))) {
      stop("CAPTURE-24 streaming preprocessing failed.")
    }
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  n <- as.integer(manifest$participants)
  m <- as.integer(manifest$grid_size)
  if (n != cfg$expected_participants || m != 86400L / cfg$epoch_seconds) {
    stop("CAPTURE-24 preprocessed dimensions do not match the prespecified 151 x 17280 design.")
  }
  raw_curve_all <- read_capture24_binary(curve_path, n, m, 4L, "numeric")
  coverage_all <- read_capture24_binary(coverage_path, n, m, 1L, "integer") > 0L
  subject_qc <- read.csv(qc_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(subject_qc$participant_id, sprintf("P%03d", seq_len(n)))) {
    stop("CAPTURE-24 participant files are not the expected P001--P151 sequence.")
  }
  public_metadata <- read.csv(
    file.path(processed_dir, "metadata.csv"), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata_index <- match(subject_qc$participant_id, public_metadata$pid)
  if (anyNA(metadata_index)) stop("CAPTURE-24 metadata does not cover every participant file.")
  subject_qc$age_group <- public_metadata$age[metadata_index]
  subject_qc$sex <- public_metadata$sex[metadata_index]
  subject_qc$female <- as.integer(subject_qc$sex == "F")
  subject_qc$wake_alignment_valid <- subject_qc$main_sleep_hours >= cfg$minimum_main_sleep_hours
  subject_qc$analysis_included <- subject_qc$wake_alignment_valid
  subject_qc$exclusion_reason <- ifelse(
    subject_qc$analysis_included, "included", "main_sleep_shorter_than_minimum"
  )
  keep <- subject_qc$analysis_included
  if (sum(keep) < 2L) stop("Fewer than two CAPTURE-24 participants have a valid wake anchor.")
  raw_curve <- raw_curve_all[keep, , drop = FALSE]
  coverage <- coverage_all[keep, , drop = FALSE]
  raw_curve[!coverage] <- NA_real_
  pointwise_mean <- colMeans(raw_curve, na.rm = TRUE)
  if (any(!is.finite(pointwise_mean))) stop("At least one wake-aligned epoch is unobserved for every participant.")
  filled <- raw_curve
  missing <- which(is.na(filled), arr.ind = TRUE)
  if (nrow(missing)) filled[missing] <- pointwise_mean[missing[, 2L]]
  analysis_qc <- subject_qc[keep, , drop = FALSE]
  dt_hours <- cfg$epoch_seconds / 3600
  list(
    weighted_X = filled * sqrt(dt_hours),
    raw_curve = raw_curve,
    filled_curve = filled,
    coverage = coverage,
    pointwise_mean = pointwise_mean,
    subject_qc = analysis_qc,
    all_subject_qc = subject_qc,
    participant_id = analysis_qc$participant_id,
    dt_hours = dt_hours,
    grid_size = m,
    processed_dir = processed_dir,
    imputation = "pointwise cohort mean; centered missing residual is zero"
  )
}

capture24_domain_fit <- function(weighted_X, hours, config) {
  m <- as.integer(hours * 3600 / config$epoch_seconds)
  X <- weighted_X[, seq_len(m), drop = FALSE]
  fit <- gram_pca(X, vectors = TRUE, center = TRUE)
  inference_result <- infer_real_components(
    fit, config$K0, config$J, config$confidence_level,
    config$certification_threshold, config$minimum_relative_gap,
    equality_draws = config$fpca_equality_draws,
    equality_seed = config$fpca_equality_seed
  )
  inference <- inference_result$table
  gaps <- inference_result$gaps
  total <- sum(fit$values)
  inference$domain_hours <- hours
  inference$grid_size <- m
  inference$n_subjects <- nrow(X)
  inference$sample_effective_rank <- spectral_effective_rank(fit$values)
  inference$effective_rank_over_n <- inference$sample_effective_rank / nrow(X)
  inference$pc1_variance_percent <- 100 * fit$values[1] / total
  inference$cumulative_pc1_5_percent <- 100 * sum(head(fit$values, 5L)) / total
  gaps$domain_hours <- hours
  gaps$grid_size <- m
  gaps$n_subjects <- nrow(X)
  spectrum_count <- min(30L, length(fit$values) - 1L)
  spectrum <- data.frame(
    domain_hours = hours,
    sample_pc = seq_len(spectrum_count),
    sample_eigenvalue = fit$values[seq_len(spectrum_count)],
    eigengap_to_next = fit$values[seq_len(spectrum_count)] - fit$values[seq_len(spectrum_count) + 1L],
    relative_eigengap_to_next =
      (fit$values[seq_len(spectrum_count)] - fit$values[seq_len(spectrum_count) + 1L]) /
      fit$values[seq_len(spectrum_count)]
  )
  list(fit = fit, inference = inference, gaps = gaps, spectrum = spectrum)
}

summarize_capture24_pilot <- function(domain_inference) {
  bind_rows_base(lapply(sort(unique(domain_inference$domain_hours)), function(hours) {
    x <- domain_inference[domain_inference$domain_hours == hours, , drop = FALSE]
    get_pc <- function(column, pc) x[[column]][match(pc, x$spike_index)]
    data.frame(
      domain_hours = hours,
      n_subjects = x$n_subjects[1],
      grid_size = x$grid_size[1],
      sample_effective_rank = x$sample_effective_rank[1],
      effective_rank_over_n = x$effective_rank_over_n[1],
      pc1_variance_percent = x$pc1_variance_percent[1],
      cumulative_pc1_5_percent = x$cumulative_pc1_5_percent[1],
      hat_Delta_pc1 = get_pc("hat_Delta", 1L),
      hat_r2_pc1 = get_pc("hat_r2", 1L),
      hat_Delta_pc2 = get_pc("hat_Delta", 2L),
      hat_r2_pc2 = get_pc("hat_r2", 2L)
    )
  }))
}

assess_capture24_pilot <- function(pilot) {
  ordered <- pilot[order(pilot$domain_hours), , drop = FALSE]
  initial <- ordered[1, , drop = FALSE]
  final <- ordered[nrow(ordered), , drop = FALSE]
  ratio_pass <- final$effective_rank_over_n >= 0.15
  rank_pass <- final$sample_effective_rank >= 8
  data.frame(
    final_domain_hours = final$domain_hours,
    n_subjects = final$n_subjects,
    initial_effective_rank = initial$sample_effective_rank,
    final_effective_rank = final$sample_effective_rank,
    effective_rank_growth_ratio = final$sample_effective_rank / initial$sample_effective_rank,
    final_effective_rank_over_n = final$effective_rank_over_n,
    final_pc1_variance_percent = final$pc1_variance_percent,
    target_ratio_min = 0.15,
    target_effective_rank_min = 8,
    implied_pc1_variance_max_percent = 100 / (0.15 * final$n_subjects),
    effective_rank_ratio_at_least_0.15 = ratio_pass,
    effective_rank_at_least_eight = rank_pass,
    overall_pass = ratio_pass && rank_pass,
    decision = if (ratio_pass && rank_pass) "proceed_to_full_analysis" else "do_not_proceed_to_full_analysis",
    stringsAsFactors = FALSE
  )
}

capture24_coverage_sensitivity <- function(prepared, config) {
  missing_minutes <- rowSums(!prepared$coverage) * config$epoch_seconds / 60
  keep <- missing_minutes <= config$coverage_sensitivity_minutes
  if (sum(keep) < 2L) return(data.frame())
  raw <- prepared$raw_curve[keep, , drop = FALSE]
  means <- colMeans(raw, na.rm = TRUE)
  if (any(!is.finite(means))) return(data.frame())
  index <- which(is.na(raw), arr.ind = TRUE)
  if (nrow(index)) raw[index] <- means[index[, 2L]]
  X <- raw * sqrt(prepared$dt_hours)
  bind_rows_base(lapply(config$domain_hours, function(hours) {
    m <- as.integer(hours * 3600 / config$epoch_seconds)
    fit <- gram_pca(X[, seq_len(m), drop = FALSE], vectors = FALSE, center = TRUE)
    effective_rank <- spectral_effective_rank(fit$values)
    data.frame(
      maximum_missing_minutes = config$coverage_sensitivity_minutes,
      domain_hours = hours,
      n_subjects = nrow(X),
      grid_size = m,
      sample_effective_rank = effective_rank,
      effective_rank_over_n = effective_rank / nrow(X),
      pc1_variance_percent = 100 * fit$values[1] / sum(fit$values),
      cumulative_pc1_5_percent = 100 * sum(head(fit$values, 5L)) / sum(fit$values)
    )
  }))
}

capture24_wake_anchor_sensitivity <- function(prepared, config) {
  keep <- prepared$subject_qc$main_sleep_hours >= config$wake_anchor_sensitivity_hours
  if (sum(keep) < 2L) stop("Wake-anchor sensitivity retained fewer than two participants.")
  raw <- prepared$raw_curve[keep, , drop = FALSE]
  means <- colMeans(raw, na.rm = TRUE)
  if (any(!is.finite(means))) stop("Wake-anchor sensitivity has an entirely unobserved epoch.")
  index <- which(is.na(raw), arr.ind = TRUE)
  if (nrow(index)) raw[index] <- means[index[, 2L]]
  weighted_X <- raw * sqrt(prepared$dt_hours)
  bind_rows_base(lapply(config$domain_hours, function(hours) {
    result <- capture24_domain_fit(weighted_X, hours, config)
    inference <- result$inference
    get_pc <- function(column, pc) inference[[column]][match(pc, inference$spike_index)]
    data.frame(
      minimum_main_sleep_hours = config$wake_anchor_sensitivity_hours,
      domain_hours = hours,
      n_subjects = nrow(weighted_X),
      grid_size = as.integer(hours * 3600 / config$epoch_seconds),
      sample_effective_rank = inference$sample_effective_rank[1],
      effective_rank_over_n = inference$effective_rank_over_n[1],
      pc1_variance_percent = inference$pc1_variance_percent[1],
      hat_r2_pc1 = get_pc("hat_r2", 1L),
      hat_r2_pc2 = get_pc("hat_r2", 2L),
      stringsAsFactors = FALSE
    )
  }))
}

capture24_cluster_diagnostics <- function(table) {
  pairs <- list(c(3L, 4L), c(5L, 6L))
  bind_rows_base(lapply(pairs, function(pair) {
    x <- table[match(pair, table$spike_index), , drop = FALSE]
    if (anyNA(x$spike_index)) return(NULL)
    data.frame(
      empirical_cluster = paste0("PC", pair[1], "-PC", pair[2]),
      first_sample_eigenvalue = x$hat_lambda[1],
      second_sample_eigenvalue = x$hat_lambda[2],
      pair_gap = x$hat_lambda[1] - x$hat_lambda[2],
      pair_relative_gap = (x$hat_lambda[1] - x$hat_lambda[2]) / x$hat_lambda[1],
      mean_individual_reliability = mean(x$hat_r2),
      both_individually_regular = all(x$regular_inference_valid),
      interpretation = "interpret_empirical_eigenspace_not_individual_directions",
      stringsAsFactors = FALSE
    )
  }))
}

capture24_wake_aligned_direction <- function(epoch_seconds) {
  hour <- (seq_len(86400L / epoch_seconds) - 0.5) * epoch_seconds / 3600
  # Contrast the first 16 hours with the final 8 hours of the circularly
  # wake-aligned day while keeping the direction orthogonal to the constant.
  value <- ifelse(hour < 16, 1, -2)
  normalize_direction(value * sqrt(epoch_seconds / 3600))
}

capture24_score_analysis <- function(fit, oriented_vectors, prepared, J) {
  raw_vectors <- lapply(seq_len(J), function(j) feature_eigenvector(fit, j))
  signs <- vapply(seq_len(J), function(j) sign(sum(raw_vectors[[j]] * oriented_vectors[[j]])), numeric(1))
  signs[signs == 0] <- 1
  scores <- sweep(fit$gram_vectors[, seq_len(J), drop = FALSE], 2L,
                  sqrt(fit$n * fit$values[seq_len(J)]) * signs, "*")
  score_table <- data.frame(
    participant_id = prepared$participant_id,
    age_group = prepared$subject_qc$age_group,
    sex = prepared$subject_qc$sex,
    mean_enmo_mg = prepared$subject_qc$mean_enmo_mg_observed,
    missing_24h_minutes = prepared$subject_qc$missing_24h_minutes,
    scores,
    check.names = FALSE
  )
  names(score_table)[(ncol(score_table) - J + 1L):ncol(score_table)] <- paste0("PC", seq_len(J), "_score")
  association <- bind_rows_base(lapply(seq_len(J), function(j) {
    score <- scores[, j]
    age_test <- kruskal.test(score ~ factor(prepared$subject_qc$age_group))
    sex_test <- t.test(score ~ factor(prepared$subject_qc$sex))
    female_mean <- mean(score[prepared$subject_qc$female == 1L])
    male_mean <- mean(score[prepared$subject_qc$female == 0L])
    data.frame(
      spike_index = j,
      correlation_mean_enmo_mg = cor(scores[, j], prepared$subject_qc$mean_enmo_mg_observed),
      correlation_missing_minutes = cor(scores[, j], prepared$subject_qc$missing_24h_minutes),
      age_group_kruskal_p_value = unname(age_test$p.value),
      female_minus_male_score = female_mean - male_mean,
      sex_t_test_p_value = unname(sex_test$p.value),
      stringsAsFactors = FALSE
    )
  }))
  list(scores = score_table, association = association)
}

write_capture24_pilot_table <- function(data, path) {
  display <- data
  names(display) <- c(
    "Hours", "n", "m", "Effective rank", "Effective rank / n",
    "PC1 variance (%)", "PC1--5 cumulative (%)",
    "Delta PC1", "r2 PC1", "Delta PC2", "r2 PC2"
  )
  write_latex_table(display, path, digits = 3L,
                    caption = "CAPTURE-24 large-domain pilot diagnostics",
                    label = "tab:capture24-pilot")
}

write_capture24_main_inference_table <- function(data, path) {
  display <- data[data$spike_index %in% 1:2,
                  c("spike_index", "hat_alpha", "hat_Delta", "hat_r2"), drop = FALSE]
  names(display) <- c("PC", "$\\hat\\alpha$", "$\\hat\\Delta$", "$\\hat r^2$")
  display$PC <- as.character(as.integer(display$PC))
  write_latex_table(
    display, path, digits = 3L,
    caption = "CAPTURE-24 leading component estimates",
    label = "tab:capture24-leading-components"
  )
}

write_capture24_wake_anchor_table <- function(data, path) {
  display <- data[, c(
    "domain_hours", "n_subjects", "sample_effective_rank",
    "effective_rank_over_n", "pc1_variance_percent", "hat_r2_pc1", "hat_r2_pc2"
  ), drop = FALSE]
  names(display) <- c("Hours", "n", "Effective rank", "Effective rank / n",
                      "PC1 variance (%)", "r2 PC1", "r2 PC2")
  write_latex_table(
    display, path, digits = 3L,
    caption = "CAPTURE-24 three-hour wake-anchor sensitivity",
    label = "tab:capture24-wake-anchor-sensitivity"
  )
}

write_capture24_cluster_table <- function(data, path) {
  display <- data[, c(
    "empirical_cluster", "pair_relative_gap", "mean_individual_reliability",
    "both_individually_regular"
  ), drop = FALSE]
  names(display) <- c("Cluster", "Relative gap", "Mean individual r2", "Both regular")
  write_latex_table(
    display, path, digits = 3L,
    caption = "CAPTURE-24 near-cluster diagnostics",
    label = "tab:capture24-near-clusters"
  )
}

capture24_fpca_vs_proposed <- function(domains, domain_spectrum) {
  leading <- domains[domains$spike_index %in% 1:2, , drop = FALSE]
  totals <- vapply(leading$domain_hours, function(hours) {
    # domain_spectrum intentionally stores only the leading 30 eigenvalues.
    # Recover the complete empirical trace from the already saved full-spectrum
    # PC1 PVE rather than treating that truncated spectrum as the denominator.
    pc1 <- domains[domains$domain_hours == hours & domains$spike_index == 1L, ]
    pc1$hat_lambda[1] / (pc1$pc1_variance_percent[1] / 100)
  }, numeric(1))
  data.frame(
    domain_hours = leading$domain_hours,
    PC = paste0("PC", leading$spike_index),
    lambda_fpca = leading$hat_lambda,
    PVE_fpca = leading$hat_lambda / totals,
    PVE_fpca_percent = 100 * leading$hat_lambda / totals,
    alpha_proposed = leading$hat_alpha,
    Delta_proposed = leading$hat_Delta,
    r2_proposed = leading$hat_r2,
    neighbor_gap = leading$relative_neighbor_gap,
    status = leading$inference_status,
    eigenvalue_correction = leading$hat_lambda - leading$hat_alpha,
    relative_eigenvalue_correction = (leading$hat_lambda - leading$hat_alpha) /
      leading$hat_alpha,
    eigenvalue_ratio = leading$hat_lambda / leading$hat_alpha,
    stringsAsFactors = FALSE
  )
}

refresh_capture24_supplementary_outputs <- function(root = PROJECT_ROOT, run_id = NULL) {
  run_id <- resolved_real_data_run_id("capture24", root, run_id)
  inference_path <- resolve_real_data_run("capture24", "inference", root, run_id)
  table <- read_required_csv(inference_path)
  cfg <- default_capture24_config()
  cfg$root <- root
  prepared <- prepare_capture24(cfg)
  wake_anchor_sensitivity <- capture24_wake_anchor_sensitivity(prepared, cfg)
  leading_components <- table[table$spike_index %in% 1:2, , drop = FALSE]
  cluster_diagnostics <- capture24_cluster_diagnostics(table)
  domains <- read_required_csv(resolve_real_data_run("capture24", "domain_expansion", root, run_id))
  domain_spectrum <- read_required_csv(resolve_real_data_run("capture24", "domain_spectrum", root, run_id))
  fpca_comparison <- capture24_fpca_vs_proposed(domains, domain_spectrum)
  prefix <- sub("__spectral_inference[.]csv$", "", inference_path)
  figure_dir <- real_data_paths(root, "capture24")$manuscript
  artifacts <- list(
    wake_anchor_sensitivity = paste0(prefix, "__wake_anchor_sensitivity.csv"),
    wake_anchor_table = paste0(prefix, "__wake_anchor_sensitivity.tex"),
    main_inference = paste0(prefix, "__leading_regular_components.csv"),
    main_table = paste0(prefix, "__leading_regular_components.tex"),
    cluster_diagnostics = paste0(prefix, "__near_cluster_diagnostics.csv"),
    cluster_table = paste0(prefix, "__near_cluster_diagnostics.tex"),
    fpca_comparison = paste0(prefix, "__fpca_vs_proposed.csv")
  )
  write_csv_atomic(wake_anchor_sensitivity, artifacts$wake_anchor_sensitivity)
  write_csv_atomic(leading_components, artifacts$main_inference)
  write_csv_atomic(cluster_diagnostics, artifacts$cluster_diagnostics)
  write_csv_atomic(fpca_comparison, artifacts$fpca_comparison)
  write_csv_atomic(fpca_comparison, file.path(dirname(inference_path), "capture24_fpca_vs_proposed.csv"))
  write_capture24_wake_anchor_table(wake_anchor_sensitivity, artifacts$wake_anchor_table)
  write_capture24_main_inference_table(leading_components, artifacts$main_table)
  write_capture24_cluster_table(cluster_diagnostics, artifacts$cluster_table)
  metadata_path <- resolve_real_data_run("capture24", "metadata", root, run_id)
  metadata <- jsonlite::read_json(metadata_path, simplifyVector = FALSE)
  metadata$wake_anchor_sensitivity <- list(
    minimum_main_sleep_hours = cfg$wake_anchor_sensitivity_hours,
    n_subjects = unique(wake_anchor_sensitivity$n_subjects),
    artifact = relative_to_root(artifacts$wake_anchor_sensitivity, root)
  )
  metadata$separation_diagnostic <- paste(
    "each adjacent gap is divided by the larger eigenvalue in its pair;",
    "both members therefore receive the same pair-gap value"
  )
  write_json_atomic(metadata, metadata_path)
  augment_real_data_run_artifacts("capture24", run_id, artifacts, root)
  invisible(artifacts)
}

run_capture24_analysis <- function(config = list()) {
  cfg <- validate_capture24_config(merge_config(default_capture24_config(), config))
  stage <- Sys.time(); pipeline_log("CAPTURE-24 | PREPROCESS START")
  prepared <- prepare_capture24(cfg)
  X <- prepared$weighted_X
  pipeline_log(sprintf("CAPTURE-24 | PREPROCESS COMPLETE | elapsed=%s | n=%d | m=%d",
                       format_elapsed(stage), nrow(X), ncol(X)))

  stage <- Sys.time(); pipeline_log("CAPTURE-24 | NESTED-DOMAIN PCA START")
  domain_results <- lapply(cfg$domain_hours, function(hours) capture24_domain_fit(X, hours, cfg))
  domains <- bind_rows_base(lapply(domain_results, `[[`, "inference"))
  domain_gaps <- bind_rows_base(lapply(domain_results, `[[`, "gaps"))
  domains$prespecified_cluster_member <-
    domains$domain_hours == 24 & domains$spike_index %in% 3:6
  domain_gaps$prespecified_cluster_member <-
    domain_gaps$domain_hours == 24 & domain_gaps$pair %in% c("PC3-PC4", "PC5-PC6")
  domain_spectrum <- bind_rows_base(lapply(domain_results, `[[`, "spectrum"))
  pilot <- summarize_capture24_pilot(domains)
  pilot_decision <- assess_capture24_pilot(pilot)
  coverage_sensitivity <- capture24_coverage_sensitivity(prepared, cfg)
  wake_anchor_sensitivity <- capture24_wake_anchor_sensitivity(prepared, cfg)
  full_fit <- domain_results[[length(domain_results)]]$fit
  direction <- capture24_wake_aligned_direction(cfg$epoch_seconds)
  inference <- infer_real_components(
    full_fit, cfg$K0, cfg$J, cfg$confidence_level,
    cfg$certification_threshold, cfg$minimum_relative_gap,
    directions = list(wake_aligned = direction), orientation_direction = direction,
    equality_draws = cfg$fpca_equality_draws,
    equality_seed = cfg$fpca_equality_seed
  )
  table <- inference$table
  primary_gaps <- inference$gaps
  capture_candidate_pairs <- c("PC3-PC4", "PC5-PC6")
  table$prespecified_cluster_member <- table$spike_index %in% 3:6
  primary_gaps$prespecified_cluster_member <- primary_gaps$pair %in% capture_candidate_pairs
  table$dataset <- "CAPTURE_24_free_living_wrist_accelerometry"
  table$n_subjects <- nrow(X)
  table$domain_hours <- 24
  table$grid_size <- ncol(X)
  table$K0_primary <- cfg$K0
  pipeline_log(sprintf("CAPTURE-24 | NESTED-DOMAIN PCA COMPLETE | elapsed=%s", format_elapsed(stage)))

  stage <- Sys.time(); pipeline_log("CAPTURE-24 | K0 SENSITIVITY AND ASSOCIATIONS START")
  sensitivity <- bind_rows_base(lapply(seq_along(domain_results), function(index) {
    result <- k0_sensitivity(
      domain_results[[index]]$fit, cfg$K0_values, cfg$J,
      cfg$confidence_level, cfg$certification_threshold, cfg$minimum_relative_gap
    )
    result$domain_hours <- cfg$domain_hours[index]
    result
  }))
  gap_sensitivity <- bind_rows_base(lapply(seq_along(domain_results), function(index) {
    result <- gap_k0_sensitivity(
      domain_results[[index]]$fit, cfg$K0_values, cfg$J, cfg$confidence_level
    )
    result$domain_hours <- cfg$domain_hours[index]
    result
  }))
  scores <- capture24_score_analysis(full_fit, inference$vectors, prepared, nrow(table))
  functions <- data.frame(
    wake_aligned_time_hours = (seq_len(ncol(X)) - 0.5) * cfg$epoch_seconds / 3600,
    mean_log1p_enmo_mg = prepared$pointwise_mean
  )
  for (j in seq_len(min(3L, nrow(table)))) {
    functions[[paste0("PC", j, "_function")]] <- inference$vectors[[j]] / sqrt(prepared$dt_hours)
  }
  pipeline_log(sprintf("CAPTURE-24 | K0 SENSITIVITY AND ASSOCIATIONS COMPLETE | elapsed=%s",
                       format_elapsed(stage)))

  leading_components <- table[table$spike_index %in% 1:2, , drop = FALSE]
  cluster_diagnostics <- capture24_cluster_diagnostics(table)
  compact_components <- compact_real_component_inference(table)
  compact_components$domain_hours <- 24
  all_domain_components <- bind_rows_base(lapply(cfg$domain_hours, function(hours) {
    x <- compact_real_component_inference(domains[domains$domain_hours == hours, , drop = FALSE])
    x$domain_hours <- hours
    x
  }))
  compact_gaps <- compact_real_gap_inference(primary_gaps)
  cluster_inference <- prespecified_cluster_inference(primary_gaps, capture_candidate_pairs)
  fpca_comparison <- capture24_fpca_vs_proposed(domains, domain_spectrum)

  cfg$n_subjects <- nrow(X); cfg$grid_size <- ncol(X)
  run_id <- make_run_id(
    "real_capture24", cfg,
    c("n_subjects", "grid_size", "domain_hours", "epoch_seconds",
      "sleep_gap_tolerance_minutes", "minimum_main_sleep_hours",
      "K0", "J", "eigengap_definition", "fpca_equality_draws",
      "fpca_equality_seed")
  )
  paths <- real_data_paths(cfg$root, "capture24")
  prefix <- file.path(paths$output, run_id)
  artifacts <- list(
    inference = paste0(prefix, "__spectral_inference.csv"),
    domain_expansion = paste0(prefix, "__domain_expansion.csv"),
    domain_spectrum = paste0(prefix, "__domain_spectrum.csv"),
    pilot_summary = paste0(prefix, "__pilot_summary.csv"),
    pilot_decision = paste0(prefix, "__pilot_decision.csv"),
    coverage_sensitivity = paste0(prefix, "__coverage_sensitivity.csv"),
    wake_anchor_sensitivity = paste0(prefix, "__wake_anchor_sensitivity.csv"),
    k0_sensitivity = paste0(prefix, "__k0_sensitivity.csv"),
    subject_qc = paste0(prefix, "__subject_qc.csv"),
    score_association = paste0(prefix, "__score_association.csv"),
    subject_scores = paste0(prefix, "__subject_scores.csv"),
    principal_functions = paste0(prefix, "__principal_functions.csv"),
    config = paste0(prefix, "__config.csv"),
    metadata = paste0(prefix, "__metadata.json"),
    table = paste0(prefix, "__spectral_inference.tex"),
    pilot_table = paste0(prefix, "__pilot_summary.tex"),
    main_inference = paste0(prefix, "__leading_regular_components.csv"),
    main_table = paste0(prefix, "__leading_regular_components.tex"),
    cluster_diagnostics = paste0(prefix, "__near_cluster_diagnostics.csv"),
    cluster_table = paste0(prefix, "__near_cluster_diagnostics.tex"),
    wake_anchor_table = paste0(prefix, "__wake_anchor_sensitivity.tex"),
    fpca_comparison = paste0(prefix, "__fpca_vs_proposed.csv")
  )
  artifacts$pc1_pc6_inference <- paste0(prefix, "__pc1_pc6_inference.csv")
  artifacts$all_domain_pc1_pc6_inference <- paste0(prefix, "__all_domain_pc1_pc6_inference.csv")
  artifacts$adjacent_gap_inference <- paste0(prefix, "__adjacent_gap_inference.csv")
  artifacts$adjacent_gap_k0_sensitivity <- paste0(prefix, "__adjacent_gap_k0_sensitivity.csv")
  artifacts$cluster_inference <- paste0(prefix, "__prespecified_candidate_cluster_inference.csv")
  write_csv_atomic(table, artifacts$inference)
  write_csv_atomic(domains, artifacts$domain_expansion)
  write_csv_atomic(domain_spectrum, artifacts$domain_spectrum)
  write_csv_atomic(pilot, artifacts$pilot_summary)
  write_csv_atomic(pilot_decision, artifacts$pilot_decision)
  write_csv_atomic(coverage_sensitivity, artifacts$coverage_sensitivity)
  write_csv_atomic(wake_anchor_sensitivity, artifacts$wake_anchor_sensitivity)
  write_csv_atomic(sensitivity, artifacts$k0_sensitivity)
  write_csv_atomic(prepared$all_subject_qc, artifacts$subject_qc)
  write_csv_atomic(scores$association, artifacts$score_association)
  write_csv_atomic(scores$scores, artifacts$subject_scores)
  write_csv_atomic(functions, artifacts$principal_functions)
  write_csv_atomic(leading_components, artifacts$main_inference)
  write_csv_atomic(cluster_diagnostics, artifacts$cluster_diagnostics)
  write_csv_atomic(fpca_comparison, artifacts$fpca_comparison)
  write_csv_atomic(compact_components, artifacts$pc1_pc6_inference)
  write_csv_atomic(all_domain_components, artifacts$all_domain_pc1_pc6_inference)
  write_csv_atomic(compact_gaps, artifacts$adjacent_gap_inference)
  write_csv_atomic(gap_sensitivity, artifacts$adjacent_gap_k0_sensitivity)
  write_csv_atomic(cluster_inference, artifacts$cluster_inference)
  write_csv_atomic(compact_components, file.path(paths$output, "capture24_pc1_pc6_inference.csv"))
  write_csv_atomic(all_domain_components, file.path(paths$output, "capture24_all_domain_pc1_pc6_inference.csv"))
  write_csv_atomic(compact_gaps, file.path(paths$output, "capture24_adjacent_gap_inference.csv"))
  write_csv_atomic(gap_sensitivity, file.path(paths$output, "capture24_adjacent_gap_k0_sensitivity.csv"))
  write_csv_atomic(cluster_inference, file.path(paths$output, "capture24_candidate_cluster_inference.csv"))
  write_csv_atomic(fpca_comparison, file.path(paths$output, "capture24_fpca_vs_proposed.csv"))
  write_csv_atomic(table, file.path(paths$output, "capture24_spectral_inference.csv"))
  write_csv_atomic(domains, file.path(paths$output, "capture24_domain_analysis.csv"))
  write_csv_atomic(sensitivity, file.path(paths$output, "capture24_k0_sensitivity.csv"))
  write_csv_atomic(coverage_sensitivity,
                   file.path(paths$output, "capture24_missingness_sensitivity.csv"))
  write_csv_atomic(wake_anchor_sensitivity,
                   file.path(paths$output, "capture24_sleep_duration_sensitivity.csv"))
  config_table <- as_config_table(cfg, "real_capture24"); config_table$run_id <- run_id
  write_csv_atomic(config_table, artifacts$config)
  metadata <- list(
    dataset = "CAPTURE-24", source = CAPTURE24_RECORD_URL, license = CAPTURE24_LICENSE,
    independent_unit = "participant; public post-QC files with a valid annotation-defined wake anchor",
    public_participant_files = cfg$expected_participants,
    n_subjects = nrow(X), excluded_wake_anchor = sum(!prepared$all_subject_qc$analysis_included),
    epoch_seconds = cfg$epoch_seconds,
    transformation = "5-second mean ENMO in mg followed by log1p",
    enmo = "max(sqrt(x^2+y^2+z^2)-1,0)*1000 mg",
    aligned_at = "end of longest annotation-defined sleep episode; circular 24-hour rotation",
    sleep_gap_tolerance_minutes = cfg$sleep_gap_tolerance_minutes,
    minimum_main_sleep_hours = cfg$minimum_main_sleep_hours,
    subject_centered = FALSE, subject_scaled = FALSE,
    across_subject_centered = TRUE, covariance_divisor = "n",
    missing_handling = prepared$imputation,
    pilot_thresholds = list(effective_rank_min = 8, effective_rank_over_n_min = 0.15),
    pilot_decision = as.list(pilot_decision[1, , drop = FALSE]),
    wake_anchor_sensitivity = paste0(
      "minimum main-sleep episode of ", cfg$wake_anchor_sensitivity_hours,
      " hours"
    ),
    analysis_role = paste(
      "formal large-domain application; the conservative pilot screen is retained",
      "for transparency and is not treated as a theorem assumption"
    )
  )
  write_json_atomic(metadata, artifacts$metadata)
  write_real_inference_table(table, artifacts$table, "CAPTURE-24 wake-aligned ENMO")
  write_capture24_pilot_table(pilot, artifacts$pilot_table)
  write_capture24_main_inference_table(leading_components, artifacts$main_table)
  write_capture24_cluster_table(cluster_diagnostics, artifacts$cluster_table)
  write_capture24_wake_anchor_table(wake_anchor_sensitivity, artifacts$wake_anchor_table)
  register_real_data_run("capture24", run_id, artifacts, cfg$root)
  plot_capture24_results(cfg$root, run_id)
  message("CAPTURE-24 analysis complete: ", run_id)
  invisible(list(table = table, pilot = pilot, pilot_decision = pilot_decision,
                 coverage_sensitivity = coverage_sensitivity,
                 wake_anchor_sensitivity = wake_anchor_sensitivity,
                 subject_qc = prepared$subject_qc, artifacts = artifacts, run_id = run_id))
}

plot_capture24_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  run_id <- resolved_real_data_run_id("capture24", root, run_id)
  domains <- read_required_csv(resolve_real_data_run("capture24", "domain_expansion", root, run_id))
  spectrum <- read_required_csv(resolve_real_data_run("capture24", "domain_spectrum", root, run_id))
  functions <- read_required_csv(resolve_real_data_run("capture24", "principal_functions", root, run_id))
  summary <- aggregate(cbind(sample_effective_rank, pc1_variance_percent) ~ domain_hours,
                       domains, function(x) x[1])
  comparison <- capture24_fpca_vs_proposed(domains, spectrum)
  paths <- real_data_paths(root, "capture24")
  write_csv_atomic(comparison, file.path(paths$output, "capture24_fpca_vs_proposed.csv"))
  outputs <- list()
  full <- spectrum[spectrum$domain_hours == 24, ]; full <- full[seq_len(min(25L, nrow(full))), ]
  outputs$spectrum <- save_real_publication_figure(
    "capture24_spectrum_24h", full, function(x) {
      plot(x$sample_pc, x$sample_eigenvalue, type = "b", pch = 16, log = "y", col = "grey35",
           xlab = "Sample PC", ylab = "Sample eigenvalue")
      points(x$sample_pc[1:2], x$sample_eigenvalue[1:2], pch = 16, col = "#0072B2", cex = 1.25)
      legend("topright", c("PC1-PC2", "Remaining spectrum"),
             col = c("#0072B2", "grey35"), pch = 16, bty = "n")
    }, root, "capture24", source_run_id = run_id, setting = "24-hour wake-aligned curve",
    quantity = "sample spectrum", comparison = "PC1-PC2 and remaining sample spectrum")
  outputs$effective_rank <- save_real_publication_figure(
    "capture24_effective_rank", summary, function(x) {
      plot(x$domain_hours, x$sample_effective_rank, type = "b", pch = 16, col = "#0072B2", lwd = 2,
           xlab = "Wake-aligned domain length (hours)", ylab = "Effective rank")
    }, root, "capture24", source_run_id = run_id, setting = "6, 12, 18, and 24 hours",
    quantity = "effective rank", comparison = "none")
  reliability <- domains[domains$spike_index %in% 1:2,
                         c("domain_hours", "spike_index", "hat_r2"), drop = FALSE]
  outputs$reliability <- save_real_publication_figure(
    "capture24_pc_reliability", reliability, function(x) {
      colors <- c("#0072B2", "#D55E00")
      plot(NA, xlim = range(x$domain_hours), ylim = range(x$hat_r2),
           xlab = "Wake-aligned domain length (hours)", ylab = "Principal-component reliability")
      for (j in 1:2) { z <- x[x$spike_index == j, ]; lines(z$domain_hours, z$hat_r2,
        type = "b", pch = 15 + j, col = colors[j], lwd = 2) }
      legend("bottomleft", paste0("PC", 1:2), col = colors, pch = 16:17, lty = 1, bty = "n")
    }, root, "capture24", source_run_id = run_id, setting = "6, 12, 18, and 24 hours",
    quantity = "r2 estimate", comparison = "PC1 and PC2")
  outputs$pve <- save_real_publication_figure(
    "capture24_pc1_pve", summary, function(x) {
      plot(x$domain_hours, x$pc1_variance_percent, type = "b", pch = 16, col = "#CC79A7", lwd = 2,
           xlab = "Wake-aligned domain length (hours)", ylab = "PC1 variance explained (%)")
    }, root, "capture24", source_run_id = run_id, setting = "6, 12, 18, and 24 hours",
    quantity = "sample FPCA PVE", comparison = "none")
  diagnostics <- domains[
    domains$domain_hours == 24 & domains$spike_index %in% 1:6,
    c("spike_index", "hat_Delta", "phase_lower_bound", "hat_r2"), drop = FALSE
  ]
  outputs$phase <- save_real_publication_figure(
    "capture24_pc1_pc6_phase_reliability", diagnostics, function(x) {
      limits <- range(c(x$hat_Delta, x$phase_lower_bound, x$hat_r2), finite = TRUE)
      padding <- max(0.01, 0.08 * diff(limits))
      plot(x$spike_index, x$hat_Delta, type = "b", pch = 16, col = "#0072B2",
           ylim = limits + c(-padding, padding), xlab = "Sample PC", ylab = "Estimate")
      lines(x$spike_index, x$phase_lower_bound, type = "b", pch = 17,
            lty = 2, col = "#D55E00")
      lines(x$spike_index, x$hat_r2, type = "b", pch = 15, col = "#009E73")
      legend("bottomleft", c("Phase margin", "Phase lower bound", "Reliability"),
             col = c("#0072B2", "#D55E00", "#009E73"), pch = c(16, 17, 15),
             lty = c(1, 2, 1), bty = "n")
    }, root, "capture24", source_run_id = run_id,
    setting = "24-hour wake-aligned curve; primary K0=8",
    quantity = "phase margin and reliability estimates", comparison = "sample PCs 1--6")
  for (j in 1:2) {
    z <- comparison[comparison$PC == paste0("PC", j), ]
    outputs[[paste0("pc", j, "_eigenvalue")]] <- save_real_publication_figure(
      paste0("capture24_pc", j, "_eigenvalue_comparison"), z, function(x) {
        ylim <- range(c(x$lambda_fpca, x$alpha_proposed))
        plot(x$domain_hours, x$lambda_fpca, type = "b", pch = 16, col = "#D55E00", lwd = 2,
             ylim = ylim, xlab = "Wake-aligned domain length (hours)", ylab = "Eigenvalue")
        lines(x$domain_hours, x$alpha_proposed, type = "b", pch = 17, col = "#0072B2", lwd = 2)
        legend("topleft", c("Sample FPCA eigenvalue", "Proposed spike estimate"),
               col = c("#D55E00", "#0072B2"), pch = 16:17, lty = 1, bty = "n")
      }, root, "capture24", source_run_id = run_id, setting = paste0("PC", j, "; 6, 12, 18, and 24 hours"),
      quantity = "sample eigenvalue and corrected spike estimate", comparison = "FPCA and Proposed")
  }
  function_names <- c("mean_log1p_enmo_mg", "PC1_function", "PC2_function", "PC3_function")
  file_names <- c("capture24_mean_function", "capture24_pc1_function", "capture24_pc2_function", "capture24_pc3_function")
  ylabels <- c("Mean log(1 + ENMO mg)", rep("Principal function", 3L))
  colors <- c("grey25", "#0072B2", "#D55E00", "#009E73")
  for (i in seq_along(function_names)) {
    z <- functions[, c("wake_aligned_time_hours", function_names[i])]; names(z)[2] <- "value"
    outputs[[file_names[i]]] <- save_real_publication_figure(
      file_names[i], z, function(x) {
        plot(x$wake_aligned_time_hours, x$value, type = "l", lwd = 2, col = colors[i],
             xlab = "Wake-aligned time (hours)", ylab = ylabels[i])
        if (i > 1L) abline(h = 0, lty = 3)
      }, root, "capture24", source_run_id = run_id, setting = "24-hour wake-aligned function",
      quantity = function_names[i], comparison = "none")
  }
  invisible(outputs)
}
