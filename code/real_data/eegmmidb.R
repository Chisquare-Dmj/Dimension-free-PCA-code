# PhysioNet EEG Motor Movement/Imagery Database (EEGMMIDB) application.
#
# One subject is one independent 64-dimensional functional observation. Run R04
# contains imagined left/right fist events: T1 is left and T2 is right. The
# analysis uses right-minus-left smoothed log 8--30 Hz power contrasts.

EEGMMIDB_BASE_URL <- "https://physionet.org/files/eegmmidb/1.0.0"
EEGMMIDB_EXCLUDED_SUBJECTS <- c(88L, 92L, 100L, 104L)

default_eegmmidb_config <- function() {
  list(
    run = 4L,
    subject_ids = setdiff(seq_len(109L), EEGMMIDB_EXCLUDED_SUBJECTS),
    band_hz = c(8, 30),
    smoothing_seconds = 0.25,
    epoch_seconds = 4,
    output_hz = 20,
    K0 = 10L,
    K0_values = c(8L, 10L, 12L),
    J = 6L,
    confidence_level = 0.95,
    certification_threshold = 0.20,
    minimum_relative_gap = 0.05,
    master_seed = MASTER_SEED,
    ncores = 1L,
    root = PROJECT_ROOT,
    use_cache = TRUE
  )
}

validate_eegmmidb_config <- function(config) {
  if (!identical(as.integer(config$run), 4L)) stop("The EEG application is defined for imagined-fist run R04.")
  if (length(config$band_hz) != 2L || config$band_hz[1] <= 0 || config$band_hz[2] <= config$band_hz[1]) {
    stop("band_hz must contain increasing positive lower and upper frequencies.")
  }
  if (config$output_hz <= 0 || config$epoch_seconds <= 0 || config$smoothing_seconds <= 0) {
    stop("output_hz, epoch_seconds, and smoothing_seconds must be positive.")
  }
  if (config$J > config$K0) stop("J cannot exceed K0.")
  if (!config$K0 %in% config$K0_values) stop("K0_values must contain the primary K0.")
  invisible(config)
}

eeg_subject_label <- function(subject_id) sprintf("S%03d", as.integer(subject_id))

eeg_edf_path <- function(subject_id, config) {
  paths <- real_data_paths(config$root, "eegmmidb")
  file.path(paths$data, "raw", sprintf("%sR%02d.edf", eeg_subject_label(subject_id), config$run))
}

download_eegmmidb <- function(config = list()) {
  cfg <- validate_eegmmidb_config(merge_config(default_eegmmidb_config(), config))
  paths <- real_data_paths(cfg$root, "eegmmidb")
  dir.create(file.path(paths$data, "raw"), recursive = TRUE, showWarnings = FALSE)
  worker <- function(subject_id) {
    label <- eeg_subject_label(subject_id)
    filename <- sprintf("%sR%02d.edf", label, cfg$run)
    destination <- file.path(paths$data, "raw", filename)
    if (file.exists(destination) && file.info(destination)$size > 1000) return(destination)
    temporary <- paste0(destination, ".part")
    url <- sprintf("%s/%s/%s", EEGMMIDB_BASE_URL, label, filename)
    status <- tryCatch({
      utils::download.file(url, temporary, mode = "wb", quiet = TRUE, method = "libcurl")
      0L
    }, error = function(error) 1L)
    if (status != 0L || !file.exists(temporary) || file.info(temporary)$size <= 1000) {
      unlink(temporary)
      stop("Failed to download ", url)
    }
    if (!file.rename(temporary, destination)) stop("Could not finalize ", destination)
    destination
  }
  files <- unlist(parallel_map(cfg$subject_ids, worker, cfg$ncores), use.names = FALSE)
  manifest <- data.frame(
    subject_id = cfg$subject_ids,
    file = basename(files),
    bytes = file.info(files)$size,
    md5 = unname(tools::md5sum(files)),
    source_url = sprintf(
      "%s/%s/%sR%02d.edf", EEGMMIDB_BASE_URL,
      vapply(cfg$subject_ids, eeg_subject_label, character(1)),
      vapply(cfg$subject_ids, eeg_subject_label, character(1)), cfg$run
    )
  )
  manifest_path <- file.path(paths$data, "raw", "download_manifest.csv")
  write_csv_atomic(manifest, manifest_path)
  message("EEGMMIDB download complete: ", length(files), " R04 files.")
  invisible(manifest_path)
}

analytic_signal <- function(x) {
  n <- length(x)
  spectrum <- fft(x)
  multiplier <- numeric(n)
  multiplier[1] <- 1
  if (n %% 2L == 0L) {
    multiplier[2:(n / 2L)] <- 2
    multiplier[n / 2L + 1L] <- 1
  } else {
    multiplier[2:((n + 1L) / 2L)] <- 2
  }
  fft(spectrum * multiplier, inverse = TRUE) / n
}

centered_moving_average <- function(x, width) {
  width <- max(1L, as.integer(width))
  if (width == 1L) return(as.numeric(x))
  smoothed <- as.numeric(stats::filter(x, rep(1 / width, width), sides = 2L))
  good <- which(is.finite(smoothed))
  if (!length(good)) stop("Smoothing window is longer than the signal.")
  smoothed[seq_len(good[1] - 1L)] <- smoothed[good[1]]
  if (good[length(good)] < length(smoothed)) {
    smoothed[(good[length(good)] + 1L):length(smoothed)] <- smoothed[good[length(good)]]
  }
  smoothed
}

band_log_power <- function(signal_matrix, sampling_rate, band_hz, smoothing_seconds) {
  required_real_data_packages("signal")
  normalized_band <- band_hz / (sampling_rate / 2)
  if (normalized_band[2] >= 1) stop("The upper band frequency must be below Nyquist.")
  filter <- signal::butter(4L, normalized_band, type = "pass")
  smoothing_width <- max(1L, round(smoothing_seconds * sampling_rate))
  output <- matrix(NA_real_, nrow(signal_matrix), ncol(signal_matrix))
  for (channel in seq_len(nrow(signal_matrix))) {
    filtered <- signal::filtfilt(filter, signal_matrix[channel, ])
    power <- Mod(analytic_signal(filtered))^2
    power <- centered_moving_average(power, smoothing_width)
    floor_value <- stats::median(power, na.rm = TRUE) * 1e-6 + 1e-18
    output[channel, ] <- log(power + floor_value)
  }
  output
}

normalize_eeg_channel <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

read_eeg_subject_contrast <- function(path, config) {
  required_real_data_packages(c("edfReader", "signal"))
  header <- edfReader::readEdfHeader(path)
  ordinary <- which(!header$sHeaders$isAnnotation)
  if (length(ordinary) < 64L) stop("Expected at least 64 ordinary EEG channels in ", basename(path))
  ordinary <- ordinary[seq_len(64L)]
  rates <- header$sHeaders$sRate[ordinary]
  if (length(unique(rates)) != 1L) stop("EEG channels do not share one sampling rate in ", basename(path))
  sampling_rate <- rates[1]
  signals <- edfReader::readEdfSignals(
    header, signals = c(ordinary, "Annotations"),
    recordStarts = FALSE, mergeASignals = TRUE
  )
  signal_matrix <- do.call(rbind, lapply(signals[seq_along(ordinary)], `[[`, "signal"))
  channel_names <- vapply(signals[seq_along(ordinary)], `[[`, character(1), "label")
  annotations_index <- which(vapply(signals, function(x) inherits(x, "ebdfASignal"), logical(1)))[1]
  if (is.na(annotations_index)) stop("EDF+ annotations are missing from ", basename(path))
  annotations <- signals[[annotations_index]]$annotations
  annotations$annotation <- trimws(annotations$annotation)

  log_power <- band_log_power(signal_matrix, sampling_rate, config$band_hz, config$smoothing_seconds)
  epoch_length <- as.integer(round(config$epoch_seconds * sampling_rate))
  extract_epochs <- function(label) {
    onsets <- annotations$onset[annotations$annotation == label]
    epochs <- lapply(onsets, function(onset) {
      start <- as.integer(round(onset * sampling_rate)) + 1L
      stop <- start + epoch_length - 1L
      if (start < 1L || stop > ncol(log_power)) return(NULL)
      log_power[, start:stop, drop = FALSE]
    })
    Filter(Negate(is.null), epochs)
  }
  left <- extract_epochs("T1")
  right <- extract_epochs("T2")
  if (length(left) < 2L || length(right) < 2L) {
    stop("Too few complete T1/T2 epochs in ", basename(path),
         ": left=", length(left), ", right=", length(right))
  }
  mean_epoch <- function(epochs) Reduce(`+`, epochs) / length(epochs)
  contrast <- mean_epoch(right) - mean_epoch(left)
  downsample_step <- as.integer(round(sampling_rate / config$output_hz))
  if (abs(sampling_rate / downsample_step - config$output_hz) > 1e-8) {
    stop("output_hz must divide the EDF sampling rate exactly.")
  }
  usable <- floor(ncol(contrast) / downsample_step) * downsample_step
  contrast <- contrast[, seq_len(usable), drop = FALSE]
  downsampled <- t(vapply(seq_len(nrow(contrast)), function(channel) {
    colMeans(matrix(contrast[channel, ], nrow = downsample_step))
  }, numeric(usable / downsample_step)))
  list(
    contrast = downsampled,
    channel_names = channel_names,
    sampling_rate = sampling_rate,
    output_hz = sampling_rate / downsample_step,
    left_epochs = length(left),
    right_epochs = length(right)
  )
}

eeg_preprocessing_cache_id <- function(config) {
  config_fingerprint(config[c("run", "band_hz", "smoothing_seconds", "epoch_seconds", "output_hz")])
}

prepare_eegmmidb <- function(config = list()) {
  cfg <- validate_eegmmidb_config(merge_config(default_eegmmidb_config(), config))
  paths <- real_data_paths(cfg$root, "eegmmidb")
  cache_id <- eeg_preprocessing_cache_id(cfg)
  prepared_path <- file.path(paths$data, "prepared", paste0("eegmmidb_r04__", cache_id, ".rds"))
  if (isTRUE(cfg$use_cache) && file.exists(prepared_path)) return(readRDS(prepared_path))
  edf_files <- vapply(cfg$subject_ids, eeg_edf_path, character(1), config = cfg)
  missing <- edf_files[!file.exists(edf_files)]
  if (length(missing)) {
    stop(length(missing), " EEG files are missing. Run dataset=eeg action=download first.")
  }
  worker <- function(index) {
    subject_id <- cfg$subject_ids[index]
    tryCatch({
      value <- read_eeg_subject_contrast(edf_files[index], cfg)
      list(ok = TRUE, subject_id = subject_id, value = value, error = NA_character_)
    }, error = function(error) {
      list(ok = FALSE, subject_id = subject_id, value = NULL, error = conditionMessage(error))
    })
  }
  processed <- parallel_map(seq_along(edf_files), worker, cfg$ncores)
  log <- data.frame(
    subject_id = vapply(processed, `[[`, integer(1), "subject_id"),
    included = vapply(processed, `[[`, logical(1), "ok"),
    error = vapply(processed, `[[`, character(1), "error"),
    stringsAsFactors = FALSE
  )
  usable <- processed[log$included]
  if (length(usable) < 30L) stop("Only ", length(usable), " usable EEG subjects remain after preprocessing.")
  reference_names <- normalize_eeg_channel(usable[[1]]$value$channel_names)
  reference_shape <- dim(usable[[1]]$value$contrast)
  for (item in usable) {
    if (!identical(dim(item$value$contrast), reference_shape)) stop("EEG contrasts have inconsistent dimensions.")
    if (!identical(normalize_eeg_channel(item$value$channel_names), reference_names)) {
      stop("EEG channel order is inconsistent across subjects.")
    }
  }
  contrasts <- lapply(usable, function(item) item$value$contrast)
  unweighted_X <- t(vapply(contrasts, function(x) as.vector(t(x)), numeric(prod(reference_shape))))
  dt <- 1 / usable[[1]]$value$output_hz
  prepared <- list(
    weighted_X = unweighted_X * sqrt(dt),
    unweighted_X = unweighted_X,
    subject_ids = vapply(usable, `[[`, integer(1), "subject_id"),
    channel_names = usable[[1]]$value$channel_names,
    normalized_channel_names = reference_names,
    n_channels = reference_shape[1],
    n_time_points = reference_shape[2],
    output_hz = usable[[1]]$value$output_hz,
    dt = dt,
    preprocessing_log = log,
    cache_id = cache_id
  )
  saveRDS(prepared, prepared_path, compress = "xz")
  prepared$prepared_path <- prepared_path
  prepared
}

eeg_lateralization_direction <- function(prepared) {
  template <- matrix(0, prepared$n_channels, prepared$n_time_points)
  active <- seq_len(prepared$n_time_points) / prepared$output_hz
  active <- active >= 0.5 & active <= 3.5
  assign_channel <- function(channel, value) {
    index <- match(normalize_eeg_channel(channel), prepared$normalized_channel_names)
    if (is.na(index)) stop("Required lateralization channel is missing: ", channel)
    template[index, active] <<- value
  }
  for (channel in c("FC3", "C3", "CP3")) assign_channel(channel, -1)
  for (channel in c("FC4", "C4", "CP4")) assign_channel(channel, 1)
  normalize_direction(as.vector(t(template)) * sqrt(prepared$dt))
}

eeg_lateralization_comparison <- function(table, pc = 2L) {
  point <- table[table$spike_index == pc, , drop = FALSE]
  data.frame(
    method = c("FPCA", "Proposed"),
    estimate = c(point$raw_lateralization[1], point$corrected_lateralization[1]),
    stringsAsFactors = FALSE
  )
}

run_eegmmidb_analysis <- function(config = list()) {
  cfg <- validate_eegmmidb_config(merge_config(default_eegmmidb_config(), config))
  prepared <- prepare_eegmmidb(cfg)
  n <- nrow(prepared$weighted_X)
  direction <- eeg_lateralization_direction(prepared)
  fit <- gram_pca(prepared$weighted_X, vectors = TRUE, center = TRUE)
  inference <- infer_real_components(
    fit, cfg$K0, cfg$J, cfg$confidence_level, cfg$certification_threshold,
    cfg$minimum_relative_gap,
    directions = list(lateralization = direction),
    orientation_direction = direction
  )
  table <- inference$table
  table$dataset <- "EEGMMIDB_R04"
  table$n_subjects <- n
  table$p_channels <- prepared$n_channels
  table$n_time_points <- prepared$n_time_points
  table$K0_primary <- cfg$K0

  sensitivity <- k0_sensitivity(
    fit, cfg$K0_values, min(cfg$J, 4L), cfg$confidence_level,
    cfg$certification_threshold, cfg$minimum_relative_gap
  )
  fpca_comparison <- eeg_lateralization_comparison(table, 2L)

  cfg$n_subjects <- n
  cfg$p_channels <- prepared$n_channels
  cfg$n_time_points <- prepared$n_time_points
  run_id <- make_run_id(
    "real_eegmmidb", cfg,
    c("run", "n_subjects", "p_channels", "n_time_points", "band_hz", "output_hz", "K0", "J")
  )
  paths <- real_data_paths(cfg$root, "eegmmidb")
  prefix <- file.path(paths$output, run_id)
  artifacts <- list(
    inference = paste0(prefix, "__spectral_inference.csv"),
    spectrum = paste0(prefix, "__sample_spectrum.csv"),
    k0_sensitivity = paste0(prefix, "__k0_sensitivity.csv"),
    preprocessing_log = paste0(prefix, "__preprocessing_log.csv"),
    config = paste0(prefix, "__config.csv"),
    metadata = paste0(prefix, "__metadata.json"),
    table = paste0(prefix, "__spectral_inference.tex")
  )
  artifacts$fpca_comparison <- paste0(prefix, "__pc2_lateralization_fpca_vs_proposed.csv")
  write_csv_atomic(table, artifacts$inference)
  write_csv_atomic(data.frame(sample_pc = seq_along(fit$values), sample_eigenvalue = fit$values), artifacts$spectrum)
  write_csv_atomic(sensitivity, artifacts$k0_sensitivity)
  write_csv_atomic(prepared$preprocessing_log, artifacts$preprocessing_log)
  write_csv_atomic(fpca_comparison, artifacts$fpca_comparison)
  write_csv_atomic(fpca_comparison, file.path(paths$output, "eeg_pc2_lateralization_fpca_vs_proposed.csv"))
  write_csv_atomic(fpca_comparison, file.path(paths$output, "eeg_fpca_vs_proposed.csv"))
  write_csv_atomic(table, file.path(paths$output, "eegmmidb_spectral_inference.csv"))
  write_csv_atomic(sensitivity, file.path(paths$output, "eegmmidb_k0_sensitivity.csv"))
  config_table <- as_config_table(cfg, "real_eegmmidb"); config_table$run_id <- run_id
  write_csv_atomic(config_table, artifacts$config)
  metadata <- list(
    dataset = "PhysioNet EEGMMIDB v1.0.0", run = "R04 imagined left/right fist",
    source = EEGMMIDB_BASE_URL, n_subjects = n,
    excluded_subjects = EEGMMIDB_EXCLUDED_SUBJECTS,
    p_channels = prepared$n_channels, n_time_points = prepared$n_time_points,
    output_hz = prepared$output_hz, sample_effective_rank = spectral_effective_rank(fit$values),
    centered = TRUE, covariance_divisor = "n", primary_K0 = cfg$K0
  )
  write_json_atomic(metadata, artifacts$metadata)
  write_real_inference_table(table, artifacts$table, "EEGMMIDB R04")
  register_real_data_run("eegmmidb", run_id, artifacts, cfg$root)
  plot_eegmmidb_results(cfg$root, run_id)
  message("EEGMMIDB analysis complete: ", run_id)
  invisible(list(table = table, fit = fit, vectors = inference$vectors, artifacts = artifacts, run_id = run_id))
}

plot_eegmmidb_results <- function(root = PROJECT_ROOT, run_id = NULL) {
  run_id <- resolved_real_data_run_id("eegmmidb", root, run_id)
  table <- read_required_csv(resolve_real_data_run("eegmmidb", "inference", root, run_id))
  spectrum <- read_required_csv(resolve_real_data_run("eegmmidb", "spectrum", root, run_id))
  outputs <- list()
  spectrum_plot <- spectrum[seq_len(min(25L, nrow(spectrum))), ]
  outputs$spectrum <- save_real_publication_figure(
    "eegmmidb_spectrum", spectrum_plot, function(x) {
      plot(x$sample_pc, x$sample_eigenvalue, type = "b", pch = 16,
           xlab = "Sample PC", ylab = "Sample eigenvalue")
    }, root, "eegmmidb", source_run_id = run_id, setting = "EEGMMIDB R04; n=105",
    quantity = "sample spectrum", comparison = "none")
  phase <- table[, c("spike_index", "hat_Delta", "phase_lower_bound", "hat_r2")]
  outputs$phase <- save_real_publication_figure(
    "eegmmidb_phase_reliability", phase, function(x) {
      limits <- range(c(x$hat_Delta, x$phase_lower_bound, x$hat_r2), finite = TRUE)
      padding <- max(0.01, 0.08 * diff(limits))
      plot(x$spike_index, x$hat_Delta, type = "b", pch = 16, col = "#0072B2",
           ylim = limits + c(-padding, padding), xlab = "Sample PC", ylab = "Estimate")
      lines(x$spike_index, x$phase_lower_bound, type = "b", pch = 17, lty = 2, col = "#D55E00")
      lines(x$spike_index, x$hat_r2, type = "b", pch = 15, col = "#009E73")
      legend("bottomleft", c("Phase margin", "Phase lower bound", "Principal-component reliability"),
             col = c("#0072B2", "#D55E00", "#009E73"), pch = c(16, 17, 15),
             lty = c(1, 2, 1), bty = "n")
    }, root, "eegmmidb", source_run_id = run_id, setting = "EEGMMIDB R04; primary K0=10",
    quantity = "phase margin and reliability estimates", comparison = "sample PCs 1--6")
  invisible(outputs)
}
