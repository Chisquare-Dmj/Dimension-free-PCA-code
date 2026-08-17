# Shared configuration, reproducibility, parallelism, and file-output helpers.

PROJECT_ROOT <- normalizePath(
  Sys.getenv("DIMENSION_FREE_PCA_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

MASTER_SEED <- 20260810L
DEFAULT_CONFIDENCE_LEVEL <- 0.95

merge_config <- function(defaults, overrides = list()) {
  if (length(overrides) == 0L) {
    return(defaults)
  }
  modifyList(defaults, overrides, keep.null = TRUE)
}

experiment_config <- function(defaults, pipeline_config, experiment_name) {
  global_names <- setdiff(names(pipeline_config), names(pipeline_config)[vapply(pipeline_config, is.list, logical(1))])
  # A command-line parameter is applied only when that experiment declares it.
  # This lets grouped selections accept parameters without polluting other runs.
  global_overrides <- pipeline_config[intersect(global_names, names(defaults))]
  local_overrides <- pipeline_config[[experiment_name]]
  if (is.null(local_overrides)) {
    local_overrides <- list()
  }
  merge_config(merge_config(defaults, global_overrides), local_overrides)
}

initialize_output_dirs <- function(root = PROJECT_ROOT) {
  paths <- list(
    replicate = file.path(root, "output", "data", "replicate"),
    truth = file.path(root, "output", "data", "truth"),
    summary = file.path(root, "output", "data", "summary"),
    manifest = file.path(root, "output", "data", "manifest"),
    figure = file.path(root, "manuscript", "simulation"),
    table = file.path(root, "output", "table"),
    log = file.path(root, "output", "log")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  invisible(paths)
}

result_path <- function(kind, filename, root = PROJECT_ROOT) {
  dirs <- initialize_output_dirs(root)
  if (!kind %in% names(dirs)) {
    stop("Unknown output kind: ", kind)
  }
  file.path(dirs[[kind]], filename)
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    stop("Could not atomically write ", path)
  }
  invisible(path)
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required saved result is missing: ", path, "\nRun the corresponding simulation mode first.")
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

as_config_table <- function(config, experiment) {
  encode <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (is.atomic(x)) {
      if (!is.null(names(x)) && any(nzchar(names(x)))) {
        return(paste(paste0(names(x), "=", x), collapse = ";"))
      }
      return(paste(x, collapse = ";"))
    }
    paste(capture.output(dput(x)), collapse = "")
  }
  data.frame(
    experiment = experiment,
    parameter = names(config),
    value = vapply(config, encode, character(1)),
    stringsAsFactors = FALSE
  )
}

sanitize_run_token <- function(x) {
  x <- gsub("[.]", "p", as.character(x))
  x <- gsub("[^A-Za-z0-9-]+", "-", x)
  gsub("(^-+|-+$)", "", x)
}

compact_config_value <- function(x) {
  if (is.null(x)) return("none")
  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    return(paste(paste0(names(x), "x", x), collapse = "-"))
  }
  if (length(x) > 6L) return(paste0(x[1], "-to-", x[length(x)], "-k", length(x)))
  paste(x, collapse = "-")
}

short_text_hash <- function(x) {
  ints <- utf8ToInt(paste(x, collapse = "|"))
  weights <- (seq_along(ints) %% 1009L) + 1L
  sprintf("%08x", as.integer(sum((ints * weights) %% 2147483629) %% 2147483629))
}

config_fingerprint <- function(config) {
  # Runtime placement and worker count do not change the scientific result.
  scientific_config <- config[setdiff(names(config), c("root", "ncores"))]
  table <- as_config_table(scientific_config, experiment = "fingerprint")
  short_text_hash(paste(table$parameter, table$value, sep = "="))
}

source_output_stem <- function(label, run_ids) {
  paste0(label, "__src-", short_text_hash(run_ids))
}

make_run_id <- function(experiment, config, fields) {
  missing <- setdiff(fields, names(config))
  if (length(missing)) stop("Cannot build run ID; missing config fields: ", paste(missing, collapse = ", "))
  pieces <- vapply(fields, function(field) {
    paste0(sanitize_run_token(field), sanitize_run_token(compact_config_value(config[[field]])))
  }, character(1))
  readable <- paste(c(experiment, pieces), collapse = "_")
  # The readable prefix identifies the main design at a glance. The fingerprint
  # prevents collisions when less visible settings such as seeds or scenario
  # subsets differ. The cap leaves room for artifact suffixes on common filesystems.
  if (nchar(readable) > 190L) readable <- substr(readable, 1L, 190L)
  run_id <- paste0(readable, "_cfg", config_fingerprint(config))
  run_id
}

run_artifact_path <- function(kind, run_id, artifact, extension = "csv", root = PROJECT_ROOT) {
  result_path(kind, paste0(run_id, "__", artifact, ".", extension), root)
}

save_config_table <- function(config, experiment, run_id, root = PROJECT_ROOT) {
  path <- run_artifact_path("truth", run_id, "config", "csv", root)
  table <- as_config_table(config, experiment)
  table$run_id <- run_id
  write_csv_atomic(table, path)
  path
}

relative_to_root <- function(path, root) {
  root <- paste0(normalizePath(root, winslash = "/", mustWork = FALSE), "/")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (startsWith(path, root)) substring(path, nchar(root) + 1L) else path
}

register_experiment_run <- function(experiment, run_id, artifacts, root = PROJECT_ROOT) {
  row <- data.frame(
    experiment = experiment,
    run_id = run_id,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
  for (name in names(artifacts)) row[[name]] <- relative_to_root(artifacts[[name]], root)
  registry_path <- result_path("manifest", "run_registry.csv", root)
  registry <- if (file.exists(registry_path)) read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  if (nrow(registry)) registry <- registry[!(registry$experiment == experiment & registry$run_id == run_id), , drop = FALSE]
  registry <- if (nrow(registry)) bind_rows_base(list(registry, row)) else row
  write_csv_atomic(registry, registry_path)
  write_csv_atomic(row, result_path("manifest", paste0(experiment, "__latest.csv"), root))
  invisible(row)
}

resolve_run_artifact <- function(experiment, artifact, root = PROJECT_ROOT, run_id = NULL) {
  manifest_path <- if (is.null(run_id)) {
    result_path("manifest", paste0(experiment, "__latest.csv"), root)
  } else {
    result_path("manifest", "run_registry.csv", root)
  }
  manifest <- read_required_csv(manifest_path)
  if (!is.null(run_id)) manifest <- manifest[manifest$experiment == experiment & manifest$run_id == run_id, , drop = FALSE]
  if (nrow(manifest) != 1L) stop("Could not uniquely resolve run for ", experiment, if (!is.null(run_id)) paste0(": ", run_id) else ".")
  if (!artifact %in% names(manifest) || is.na(manifest[[artifact]]) || !nzchar(manifest[[artifact]])) {
    stop("Artifact '", artifact, "' is not registered for run ", manifest$run_id)
  }
  path <- manifest[[artifact]]
  if (!grepl("^/", path)) path <- file.path(root, path)
  if (!file.exists(path)) stop("Registered artifact is missing: ", path)
  path
}

resolved_run_id <- function(experiment, root = PROJECT_ROOT, run_id = NULL) {
  manifest_path <- if (is.null(run_id)) result_path("manifest", paste0(experiment, "__latest.csv"), root) else result_path("manifest", "run_registry.csv", root)
  manifest <- read_required_csv(manifest_path)
  if (!is.null(run_id)) manifest <- manifest[manifest$experiment == experiment & manifest$run_id == run_id, , drop = FALSE]
  if (nrow(manifest) != 1L) stop("Could not uniquely resolve run ID for ", experiment)
  manifest$run_id
}

replication_seed <- function(experiment_id, replication_id, master_seed = MASTER_SEED) {
  as.integer(master_seed + 100000L * as.integer(experiment_id) + as.integer(replication_id))
}

with_fixed_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expression)
}

parallel_map <- function(indices, FUN, ncores = 1L) {
  ncores <- max(1L, min(as.integer(ncores), length(indices)))
  if (ncores == 1L || .Platform$OS.type == "windows") {
    return(lapply(indices, FUN))
  }
  parallel::mclapply(indices, FUN, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = FALSE)
}

bind_rows_base <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (length(x) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(x, names), use.names = FALSE))
  normalized <- lapply(x, function(item) {
    missing <- setdiff(columns, names(item))
    for (name in missing) item[[name]] <- NA
    item[columns]
  })
  out <- do.call(rbind, normalized)
  rownames(out) <- NULL
  out
}

figure_method_scope <- function(quantity, comparison) {
  text <- tolower(paste(quantity, comparison))
  has_fpca <- grepl(
    "fpca|sample spectr|sample eigenvalue|principal function|individual alignment|effective rank|variance explained",
    text
  )
  has_proposed <- grepl(
    "proposed|phase margin|reliability|wald|studentized|k0-sensitivity|cumulant|population-spike",
    text
  )
  if (has_fpca && has_proposed) return("Both")
  if (has_fpca) return("FPCA")
  "Proposed"
}

refresh_figure_output_summary <- function(root = PROJECT_ROOT) {
  specifications <- list(
    simulation = file.path(root, "output", "data", "summary",
                           "simulation_figure_output_manifest.csv"),
    real_data = file.path(root, "output", "real_data", "figure_output_manifest.csv")
  )
  pieces <- lapply(names(specifications), function(source) {
    path <- specifications[[source]]
    if (!file.exists(path)) return(NULL)
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    x$output_group <- source
    if (!"dataset" %in% names(x)) x$dataset <- if (source == "simulation") "simulation" else NA_character_
    if (!"source_run_ids" %in% names(x) && "source_run_id" %in% names(x)) {
      x$source_run_ids <- x$source_run_id
    }
    if (!"method_scope" %in% names(x)) {
      x$method_scope <- mapply(figure_method_scope, x$quantity, x$comparison)
    }
    x$analysis_setting <- x$setting
    x[, c("output_group", "dataset", "filename", "plot_data_file", "analysis_setting",
          "quantity", "method_scope", "comparison", "source_run_ids"), drop = FALSE]
  })
  summary <- bind_rows_base(pieces)
  path <- file.path(root, "output", "figure_output_summary.csv")
  write_csv_atomic(summary, path)
  invisible(path)
}

generate_standard_scores <- function(distribution, nrow, ncol) {
  distribution <- tolower(distribution)
  if (distribution %in% c("gaussian", "gau", "normal")) {
    return(matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol))
  }
  if (distribution %in% c("uniform", "unif")) {
    return(matrix(runif(nrow * ncol, -sqrt(3), sqrt(3)), nrow = nrow, ncol = ncol))
  }
  if (distribution %in% c("t12", "t")) {
    return(matrix(rt(nrow * ncol, df = 12) * sqrt(5 / 6), nrow = nrow, ncol = ncol))
  }
  stop("Unsupported score distribution: ", distribution)
}

score_fourth_cumulant <- function(distribution) {
  distribution <- tolower(distribution)
  if (distribution %in% c("gaussian", "gau", "normal")) return(0)
  if (distribution %in% c("uniform", "unif")) return(-6 / 5)
  if (distribution %in% c("t12", "t")) return(6 / (12 - 4))
  stop("Unsupported score distribution: ", distribution)
}

matrix_sqrt_psd <- function(A, tolerance = 1e-12) {
  eig <- eigen(A, symmetric = TRUE)
  values <- pmax(eig$values, tolerance * max(1, eig$values[1]))
  sweep(eig$vectors, 2L, sqrt(values), "*") %*% t(eig$vectors)
}

block_diagonal <- function(blocks) {
  row_sizes <- vapply(blocks, nrow, integer(1))
  col_sizes <- vapply(blocks, ncol, integer(1))
  out <- matrix(0, sum(row_sizes), sum(col_sizes))
  row_start <- 1L
  col_start <- 1L
  for (i in seq_along(blocks)) {
    rows <- row_start:(row_start + row_sizes[i] - 1L)
    cols <- col_start:(col_start + col_sizes[i] - 1L)
    out[rows, cols] <- blocks[[i]]
    row_start <- max(rows) + 1L
    col_start <- max(cols) + 1L
  }
  out
}

fixed_haar_matrix <- function(dimension, seed) {
  with_fixed_seed(seed, {
    qr_result <- qr(matrix(rnorm(dimension * dimension), dimension, dimension))
    Q <- qr.Q(qr_result)
    signs <- sign(diag(qr.R(qr_result)))
    signs[signs == 0] <- 1
    sweep(Q, 2L, signs, "*")
  })
}

fixed_orthonormal_directions <- function(dimension, count, seed) {
  with_fixed_seed(seed, qr.Q(qr(matrix(rnorm(dimension * count), dimension, count))))
}

progress_message <- function(label, current, total) {
  cat(sprintf("[%s] %d/%d (%.1f%%)\n", label, current, total, 100 * current / total))
  flush.console()
}

validate_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer.")
  }
  as.integer(x)
}

STANDARD_REPLICATE_COLUMNS <- c(
  "experiment", "run_id", "scenario", "score_distribution", "n", "p", "T", "N", "m",
  "replication", "replications_planned", "spike_index", "M", "K0", "decay",
  "master_seed", "population_seed", "confidence_level", "certification_threshold",
  "population_phase_status", "regular_inference_valid",
  "true_alpha", "true_psi", "true_Delta", "true_r2",
  "hat_lambda", "hat_alpha", "hat_Delta", "hat_r2",
  "se_alpha", "se_Delta", "se_r2",
  "ci_alpha_lower", "ci_alpha_upper", "ci_Delta_lower", "ci_Delta_upper",
  "ci_r2_lower", "ci_r2_upper", "cover_alpha", "cover_Delta", "cover_r2",
  "actual_signal_overlap", "theta_true", "theta_raw", "theta_corrected",
  "Q_true", "Q_raw", "Q_corrected", "individual_alignment_repeated"
)

standardize_replicate_output <- function(data) {
  for (column in setdiff(STANDARD_REPLICATE_COLUMNS, names(data))) data[[column]] <- NA
  data[c(STANDARD_REPLICATE_COLUMNS, setdiff(names(data), STANDARD_REPLICATE_COLUMNS))]
}

add_run_metadata <- function(data, config, run_id, replications_planned,
                             population_seed = NA_integer_) {
  data$run_id <- run_id
  data$replications_planned <- replications_planned
  data$M <- if (!is.null(config$M)) config$M else NA_integer_
  data$decay <- if (!is.null(config$decay)) config$decay else NA_real_
  data$master_seed <- config$master_seed
  data$population_seed <- population_seed
  data$confidence_level <- config$confidence_level
  data$certification_threshold <- config$certification_threshold
  data
}

validate_inference_config <- function(config, require_dimensions = character()) {
  for (name in c("M", "K0", require_dimensions)) {
    if (!is.null(config[[name]])) validate_positive_integer(config[[name]], name)
  }
  M <- config[["M"]]
  K0 <- config[["K0"]]
  n <- config[["n"]]
  replications <- config[["replications"]]
  confidence_level <- config[["confidence_level"]]
  if (!is.null(M) && !is.null(K0) && K0 < M) stop("K0 must be at least M.")
  if (!is.null(n) && !is.null(K0) && K0 >= n) stop("K0 must be smaller than n.")
  if (!is.null(replications)) validate_positive_integer(replications, "replications")
  if (!is.null(confidence_level) && (confidence_level <= 0 || confidence_level >= 1)) stop("confidence_level must lie strictly between 0 and 1.")
  invisible(config)
}
