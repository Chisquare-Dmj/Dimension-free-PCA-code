# Dimension-free PCA estimators from the deflated empirical spectrum.

gram_pca <- function(X, vectors = TRUE, center = FALSE, tolerance = 1e-10) {
  # Simulations use center=FALSE because observations are generated with mean
  # zero. Real data should normally use center=TRUE. The covariance divisor is
  # always n, including after centering, to match the theoretical normalization.
  if (!is.logical(center) || length(center) != 1L || is.na(center)) stop("center must be TRUE or FALSE.")
  column_means <- if (center) colMeans(X) else rep(0, ncol(X))
  if (center) X <- sweep(X, 2L, column_means, "-")
  n <- nrow(X)
  gram <- tcrossprod(X) / n
  eig <- eigen(gram, symmetric = TRUE, only.values = !vectors)
  values <- pmax(eig$values, 0)
  result <- list(values = values, gram_vectors = if (vectors) eig$vectors else NULL,
                 n = n, X = X, centered = center, column_means = column_means)
  class(result) <- "dfpca_fit"
  result
}

feature_eigenvector <- function(fit, index, tolerance = 1e-10) {
  lambda <- fit$values[index]
  if (!is.finite(lambda) || lambda <= tolerance * max(1, fit$values[1])) {
    stop("Sample eigenvalue ", index, " is numerically zero; its feature direction is undefined.")
  }
  as.numeric(crossprod(fit$X, fit$gram_vectors[, index]) / sqrt(fit$n * lambda))
}

feature_projection <- function(fit, index, direction) {
  lambda <- fit$values[index]
  if (!is.finite(lambda) || lambda <= 0) return(NA_real_)
  scores <- as.numeric(fit$X %*% direction)
  sum(fit$gram_vectors[, index] * scores) / sqrt(fit$n * lambda)
}

population_spike_truth <- function(population_eigenvalues, M, n, spike_indices = seq_len(M)) {
  population_eigenvalues <- sort(as.numeric(population_eigenvalues), decreasing = TRUE)
  if (M >= length(population_eigenvalues)) stop("M must be smaller than the population spectrum length.")
  bulk <- population_eigenvalues[(M + 1L):length(population_eigenvalues)]
  rows <- lapply(spike_indices, function(j) {
    alpha <- population_eigenvalues[j]
    gaps <- alpha - bulk
    psi <- alpha * (1 + sum(bulk / gaps) / n)
    Delta <- 1 - sum((bulk / gaps)^2) / n
    r2 <- alpha / psi * Delta
    data.frame(spike_index = j, true_alpha = alpha, true_psi = psi, true_Delta = Delta, true_r2 = r2)
  })
  bind_rows_base(rows)
}

deflated_spectral_moments <- function(sample_eigenvalues, target_index, K0, n, orders = 1:4) {
  values <- sort(as.numeric(sample_eigenvalues), decreasing = TRUE)
  if (K0 >= length(values)) stop("K0 must be smaller than the number of sample eigenvalues.")
  z <- values[target_index]
  deflated_spectral_moments_at(sample_eigenvalues, z, K0, n, orders)
}

deflated_spectral_moments_at <- function(sample_eigenvalues, z, K0, n, orders = 1:4) {
  values <- sort(as.numeric(sample_eigenvalues), decreasing = TRUE)
  if (K0 >= length(values)) stop("K0 must be smaller than the number of sample eigenvalues.")
  if (!is.finite(z)) stop("The transform evaluation point must be finite.")
  bulk <- values[(K0 + 1L):length(values)]
  gaps <- bulk - z
  if (any(abs(gaps) <= .Machine$double.eps * max(1, abs(z)))) {
    stop("A deflated bulk eigenvalue is numerically tied with the target sample eigenvalue.")
  }
  moments <- vapply(orders, function(s) sum(gaps^(-s)) / n, numeric(1))
  names(moments) <- paste0("m", orders)
  moments
}

fpca_scores <- function(fit, target_index, tolerance = 1e-10) {
  if (is.null(fit$gram_vectors)) stop("Classical FPCA inference requires Gram eigenvectors.")
  lambda <- fit$values[target_index]
  scores <- sqrt(fit$n * lambda) * fit$gram_vectors[, target_index]
  score_identity_error <- abs(mean(scores^2) - lambda)
  if (score_identity_error > tolerance * max(1, abs(lambda))) {
    stop("Gram-score identity failed for sample PC ", target_index, ".")
  }
  attr(scores, "identity_error") <- score_identity_error
  scores
}

# Classical fixed-covariance FPCA inference follows the empirical-score
# expansion of Dauxois, Pousse and Romain (1982) and Hall and Hosseini-Nasab
# (2006). It is intentionally independent of the dimension-free correction.
fpca_eigenvalue_interval <- function(fit, target_index,
                                     confidence_level = DEFAULT_CONFIDENCE_LEVEL,
                                     tolerance = 1e-10) {
  lambda <- fit$values[target_index]
  scores <- fpca_scores(fit, target_index, tolerance)
  score_variance <- stats::var(scores^2)
  standard_error <- sqrt(score_variance / fit$n)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  data.frame(
    fpca_alpha_hat = lambda,
    fpca_var_lambda = score_variance,
    fpca_alpha_se = standard_error,
    fpca_alpha_ci_lower = lambda - z * standard_error,
    fpca_alpha_ci_upper = lambda + z * standard_error,
    fpca_alpha_ci_length = 2 * z * standard_error,
    fpca_score_identity_error = attr(scores, "identity_error")
  )
}

symmetric_relative_gap <- function(upper, lower) {
  2 * (upper - lower) / (upper + lower)
}

legacy_asymmetric_relative_gap <- function(upper, lower) {
  (upper - lower) / upper
}

.fpca_null_draw_cache <- new.env(parent = emptyenv())

fpca_null_standard_draws <- function(B, seed) {
  key <- paste(as.integer(B), as.integer(seed), sep = "_")
  if (!exists(key, envir = .fpca_null_draw_cache, inherits = FALSE)) {
    draws <- with_fixed_seed(seed, matrix(stats::rnorm(2L * B), ncol = 2L))
    assign(key, draws, envir = .fpca_null_draw_cache)
  }
  get(key, envir = .fpca_null_draw_cache, inherits = FALSE)
}

# Classical repeated-root inference. This implements the fixed-covariance
# Gaussian perturbation-block limit and Anderson's Gaussian likelihood-ratio
# benchmark; it does not call or approximate the proposed high-complexity law.
fpca_pair_inference <- function(fit, j,
                                confidence_level = DEFAULT_CONFIDENCE_LEVEL,
                                equality_draws = 100000L,
                                equality_seed = 82471L,
                                known_distinct = NULL,
                                tolerance = 1e-10) {
  j <- as.integer(j)
  if (j < 1L || j + 1L > length(fit$values)) {
    stop("j must identify an adjacent empirical FPCA pair.")
  }
  equality_draws <- as.integer(equality_draws)
  if (equality_draws < 100000L) stop("The general FPCA equality test requires at least 100000 draws.")
  lambda_upper <- fit$values[j]
  lambda_lower <- fit$values[j + 1L]
  lambda_pool <- (lambda_upper + lambda_lower) / 2
  score_upper <- fpca_scores(fit, j, tolerance)
  score_lower <- fpca_scores(fit, j + 1L, tolerance)
  score_squares <- cbind(score_upper^2, score_lower^2)
  score_covariance <- stats::cov(score_squares)
  gradient <- c(
    4 * lambda_lower / (lambda_upper + lambda_lower)^2,
    -4 * lambda_upper / (lambda_upper + lambda_lower)^2
  )
  gap_variance <- as.numeric(crossprod(gradient, score_covariance %*% gradient)) / fit$n
  gap_se <- sqrt(max(0, gap_variance))
  gap <- symmetric_relative_gap(lambda_upper, lambda_lower)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  ci_raw <- gap + c(-1, 1) * z * gap_se

  block_contributions <- cbind(
    score_upper^2 - lambda_pool,
    score_upper * score_lower,
    score_lower^2 - lambda_pool
  )
  block_covariance <- stats::cov(block_contributions)
  # The eigengap of a 2x2 symmetric block is
  # sqrt((Z11-Z22)^2 + (2 Z12)^2), so only this 2D transform is needed.
  transform_matrix <- rbind(c(1, 0, -1), c(0, 2, 0))
  spacing_covariance <- transform_matrix %*% block_covariance %*% t(transform_matrix)
  decomposition <- eigen(spacing_covariance, symmetric = TRUE)
  square_root <- decomposition$vectors %*%
    diag(sqrt(pmax(decomposition$values, 0)), 2L) %*% t(decomposition$vectors)
  null_coordinates <- fpca_null_standard_draws(equality_draws, equality_seed) %*% square_root
  null_spacing <- sqrt(rowSums(null_coordinates^2)) / lambda_pool
  observed_spacing <- sqrt(fit$n) * gap
  p_general <- (1 + sum(null_spacing >= observed_spacing)) / (equality_draws + 1)

  anderson_Q <- fit$n * (
    2 * log(lambda_pool) - log(lambda_upper) - log(lambda_lower)
  )
  alpha <- 1 - confidence_level
  regular <- if (is.null(known_distinct)) p_general < alpha else isTRUE(known_distinct)
  data.frame(
    fpca_gap_sym = gap,
    fpca_gap_legacy_asymmetric = legacy_asymmetric_relative_gap(lambda_upper, lambda_lower),
    fpca_gap_wald_se = gap_se,
    fpca_gap_wald_ci_raw_lower = ci_raw[1],
    fpca_gap_wald_ci_raw_upper = ci_raw[2],
    fpca_gap_wald_ci_lower = max(0, ci_raw[1]),
    fpca_gap_wald_ci_upper = min(2, ci_raw[2]),
    fpca_gap_wald_regular = regular,
    fpca_p_equal_general = p_general,
    fpca_equality_draws = equality_draws,
    fpca_anderson_Q = anderson_Q,
    fpca_anderson_p_equal = stats::pchisq(anderson_Q, df = 2, lower.tail = FALSE),
    fpca_pair_score_identity_error = max(
      attr(score_upper, "identity_error"), attr(score_lower, "identity_error")
    )
  )
}

invert_noncentral_chisq_gap_ci <- function(observed_squared_statistic,
                                           effective_nscale,
                                           confidence_level = DEFAULT_CONFIDENCE_LEVEL,
                                           maximum_ncp = 1e8) {
  x <- as.numeric(observed_squared_statistic)
  nscale <- as.numeric(effective_nscale)
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("observed_squared_statistic must be one finite nonnegative number.")
  }
  if (length(nscale) != 1L || !is.finite(nscale) || nscale <= 0) {
    stop("effective_nscale must be one finite positive number.")
  }
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.")
  }
  alpha <- 1 - confidence_level
  lower_tail <- alpha / 2
  upper_tail <- 1 - alpha / 2
  central_cdf <- stats::pchisq(x, df = 2, ncp = 0)

  # For a sufficiently small observation, no nonnegative noncentrality is in
  # the equal-tail confidence set. Preserve that empty-set information rather
  # than incorrectly turning it into an interval containing zero.
  if (central_cdf < lower_tail) {
    return(data.frame(
      gap_ci_lower = NA_real_, gap_ci_upper = NA_real_,
      ncp_ci_lower = NA_real_, ncp_ci_upper = NA_real_,
      ci_empty = TRUE, ci_contains_zero = FALSE
    ))
  }

  solve_ncp <- function(target_probability) {
    objective <- function(ncp) stats::pchisq(x, df = 2, ncp = ncp) - target_probability
    upper <- 1
    while (objective(upper) > 0 && upper < maximum_ncp) upper <- min(maximum_ncp, 2 * upper)
    if (objective(upper) > 0) return(Inf)
    stats::uniroot(objective, interval = c(0, upper), tol = 1e-10)$root
  }
  ncp_lower <- if (central_cdf <= upper_tail) 0 else solve_ncp(upper_tail)
  ncp_upper <- solve_ncp(lower_tail)
  gap_from_ncp <- function(ncp) if (is.infinite(ncp)) Inf else 2 * sqrt(max(0, ncp) / nscale)
  data.frame(
    gap_ci_lower = gap_from_ncp(ncp_lower),
    gap_ci_upper = gap_from_ncp(ncp_upper),
    ncp_ci_lower = ncp_lower,
    ncp_ci_upper = ncp_upper,
    ci_empty = FALSE,
    ci_contains_zero = ncp_lower == 0
  )
}

invert_noncentral_chisq_gap_upper <- function(observed_squared_statistic,
                                              effective_nscale,
                                              confidence_level = DEFAULT_CONFIDENCE_LEVEL,
                                              maximum_ncp = 1e8) {
  x <- as.numeric(observed_squared_statistic)
  nscale <- as.numeric(effective_nscale)
  if (!is.finite(x) || length(x) != 1L || x < 0 ||
      !is.finite(nscale) || length(nscale) != 1L || nscale <= 0) {
    stop("The observed statistic and effective sample scale must be finite and valid.")
  }
  alpha <- 1 - confidence_level
  objective <- function(ncp) stats::pchisq(x, df = 2, ncp = ncp) - alpha
  point_gap <- 2 * sqrt(x / nscale)
  if (objective(0) <= 0) {
    # The literal lower-tail inversion is empty for an unusually tiny radial
    # observation. A Gaussian confidence ball gives a conservative, finite
    # boundary-aware upper bound instead of the misleading value zero.
    radius <- sqrt(stats::qchisq(confidence_level, df = 2))
    return(2 * (sqrt(x) + radius) / sqrt(nscale))
  }
  upper <- 1
  while (objective(upper) > 0 && upper < maximum_ncp) upper <- min(maximum_ncp, 2 * upper)
  ncp_upper <- if (objective(upper) > 0) Inf else {
    stats::uniroot(objective, c(0, upper), tol = 1e-10)$root
  }
  if (is.infinite(ncp_upper)) Inf else max(point_gap, 2 * sqrt(ncp_upper / nscale))
}

multiplicity2_gap_inference <- function(sample_eigenvalues, j, n, K0,
                                        confidence_level = DEFAULT_CONFIDENCE_LEVEL) {
  values <- sort(as.numeric(sample_eigenvalues), decreasing = TRUE)
  j <- as.integer(j)
  if (j < 1L || j + 1L > length(values)) stop("j must identify an adjacent sample-eigenvalue pair.")
  if (K0 < j + 1L) stop("K0 must exclude both members of the candidate pair.")
  lambda_upper <- values[j]
  lambda_lower <- values[j + 1L]
  lambda_pool <- (lambda_upper + lambda_lower) / 2
  m_upper <- deflated_spectral_moments_at(values, lambda_upper, K0, n, 1L)
  m_lower <- deflated_spectral_moments_at(values, lambda_lower, K0, n, 1L)
  m_pool <- deflated_spectral_moments_at(values, lambda_pool, K0, n, 1:4)
  alpha_upper <- -1 / m_upper[["m1"]]
  alpha_lower <- -1 / m_lower[["m1"]]
  alpha_pool <- -1 / m_pool[["m1"]]
  Delta_pool <- m_pool[["m1"]]^2 / m_pool[["m2"]]
  r2_pool <- -m_pool[["m1"]] / (lambda_pool * m_pool[["m2"]])
  sample_gap_sym <- symmetric_relative_gap(lambda_upper, lambda_lower)
  proposed_gap_sym <- symmetric_relative_gap(alpha_upper, alpha_lower)
  sample_gap_legacy <- legacy_asymmetric_relative_gap(lambda_upper, lambda_lower)
  proposed_gap_legacy <- legacy_asymmetric_relative_gap(alpha_upper, alpha_lower)
  valid <- all(is.finite(c(alpha_upper, alpha_lower, alpha_pool, Delta_pool, r2_pool))) &&
    alpha_upper > 0 && alpha_lower > 0 && Delta_pool > 0 && r2_pool > 0

  if (!valid) {
    return(data.frame(
      sample_lambda_upper = lambda_upper, sample_lambda_lower = lambda_lower,
      sample_gap_sym = sample_gap_sym, raw_relative_gap = sample_gap_legacy,
      raw_relative_gap_legacy = sample_gap_legacy, lambda_pool = lambda_pool,
      alpha_hat_upper = alpha_upper, alpha_hat_lower = alpha_lower,
      alpha_pool = alpha_pool, proposed_gap_sym = proposed_gap_sym,
      corrected_relative_gap = proposed_gap_legacy,
      corrected_relative_gap_legacy = proposed_gap_legacy,
      Delta_pool = Delta_pool, r2_pool = r2_pool,
      T_sample_HC = NA_real_, p_sample_HC = NA_real_,
      T_proposed = NA_real_, T2_proposed = NA_real_, p_proposed = NA_real_,
      proposed_gap_ci_raw_lower = NA_real_, proposed_gap_ci_raw_upper = NA_real_,
      proposed_gap_ci_raw_empty = NA, proposed_gap_upper95 = NA_real_,
      proposed_gap_upper95_boundary_fallback = NA,
      proposed_gap_ci_boundary_lower = NA_real_,
      proposed_gap_ci_boundary_upper = NA_real_,
      proposed_gap_local_inference_valid = FALSE, inference_valid = FALSE
    ))
  }

  T_sample_HC <- sqrt(n * Delta_pool) * sample_gap_sym / (2 * r2_pool)
  T_proposed <- sqrt(n * Delta_pool) * proposed_gap_sym / 2
  proposed_ci <- invert_noncentral_chisq_gap_ci(T_proposed^2, n * Delta_pool, confidence_level)
  proposed_upper <- invert_noncentral_chisq_gap_upper(
    T_proposed^2, n * Delta_pool, confidence_level
  )
  upper_boundary_fallback <-
    stats::pchisq(T_proposed^2, df = 2, ncp = 0) <= 1 - confidence_level
  boundary_lower <- if (proposed_ci$ci_empty) 0 else proposed_ci$gap_ci_lower
  boundary_upper <- if (proposed_ci$ci_empty) proposed_upper else proposed_ci$gap_ci_upper
  data.frame(
    sample_lambda_upper = lambda_upper, sample_lambda_lower = lambda_lower,
    sample_gap_sym = sample_gap_sym, raw_relative_gap = sample_gap_legacy,
    raw_relative_gap_legacy = sample_gap_legacy, lambda_pool = lambda_pool,
    alpha_hat_upper = alpha_upper, alpha_hat_lower = alpha_lower,
    alpha_pool = alpha_pool, proposed_gap_sym = proposed_gap_sym,
    corrected_relative_gap = proposed_gap_legacy,
    corrected_relative_gap_legacy = proposed_gap_legacy,
    Delta_pool = Delta_pool, r2_pool = r2_pool,
    T_sample_HC = T_sample_HC,
    p_sample_HC = stats::pchisq(T_sample_HC^2, df = 2, lower.tail = FALSE),
    T_proposed = T_proposed, T2_proposed = T_proposed^2,
    p_proposed = stats::pchisq(T_proposed^2, df = 2, lower.tail = FALSE),
    proposed_gap_ci_raw_lower = proposed_ci$gap_ci_lower,
    proposed_gap_ci_raw_upper = proposed_ci$gap_ci_upper,
    proposed_gap_ci_raw_empty = proposed_ci$ci_empty,
    proposed_gap_upper95 = proposed_upper,
    proposed_gap_upper95_boundary_fallback = upper_boundary_fallback,
    proposed_gap_ci_boundary_lower = boundary_lower,
    proposed_gap_ci_boundary_upper = boundary_upper,
    proposed_gap_local_inference_valid = TRUE,
    inference_valid = TRUE
  )
}

pooled_cluster_inference <- function(sample_eigenvalues, cluster_indices, K0, n,
                                     confidence_level = DEFAULT_CONFIDENCE_LEVEL) {
  indices <- sort(unique(as.integer(cluster_indices)))
  if (length(indices) < 2L || any(diff(indices) != 1L)) stop("cluster_indices must be consecutive and contain at least two PCs.")
  values <- sort(as.numeric(sample_eigenvalues), decreasing = TRUE)
  location <- mean(values[indices])
  m <- deflated_spectral_moments_at(values, location, K0, n, 1:4)
  alpha <- -1 / m[["m1"]]
  Delta <- m[["m1"]]^2 / m[["m2"]]
  r2 <- -m[["m1"]] / (location * m[["m2"]])
  d_Delta <- 2 * m[["m1"]] - 2 * m[["m1"]]^2 * m[["m3"]] / m[["m2"]]^2
  d_r2 <- r2 * (m[["m2"]] / m[["m1"]] - 1 / location - 2 * m[["m3"]] / m[["m2"]])
  multiplicity <- length(indices)
  se_alpha <- sqrt(2 * alpha^2 / (n * Delta)) / sqrt(multiplicity)
  se_Delta <- sqrt(2 * d_Delta^2 / (n * m[["m2"]])) / sqrt(multiplicity)
  se_r2 <- sqrt(2 * d_r2^2 / (n * m[["m2"]])) / sqrt(multiplicity)
  z_two <- qnorm(1 - (1 - confidence_level) / 2)
  z_one <- qnorm(confidence_level)
  data.frame(
    cluster = paste0("PC", min(indices), "-PC", max(indices)),
    multiplicity = multiplicity, lambda_pool = location,
    pooled_alpha = alpha, pooled_alpha_se = se_alpha,
    pooled_alpha_ci_lower = alpha - z_two * se_alpha,
    pooled_alpha_ci_upper = alpha + z_two * se_alpha,
    pooled_Delta = Delta, pooled_Delta_se = se_Delta,
    pooled_Delta_ci_lower = Delta - z_two * se_Delta,
    pooled_Delta_ci_upper = Delta + z_two * se_Delta,
    pooled_Delta_one_sided_lower = Delta - z_one * se_Delta,
    pooled_r2 = r2, pooled_r2_se = se_r2,
    pooled_r2_ci_lower = r2 - z_two * se_r2,
    pooled_r2_ci_upper = r2 + z_two * se_r2,
    pooled_inference_valid = all(is.finite(c(alpha, Delta, r2, se_alpha, se_Delta, se_r2))) && Delta > 0 && r2 > 0
  )
}

simulate_goe_spacing <- function(m, adjacent_rank = 1L, B = 10000L, seed = MASTER_SEED) {
  m <- as.integer(m); adjacent_rank <- as.integer(adjacent_rank); B <- as.integer(B)
  if (m < 2L) stop("m must be at least two.")
  if (adjacent_rank < 1L || adjacent_rank >= m) stop("adjacent_rank must lie between 1 and m-1.")
  if (B < 1L) stop("B must be positive.")
  set.seed(seed)
  vapply(seq_len(B), function(iteration) {
    G <- matrix(0, m, m)
    diag(G) <- rnorm(m, sd = sqrt(2))
    upper <- upper.tri(G)
    G[upper] <- rnorm(sum(upper))
    G <- G + t(G)
    diag(G) <- diag(G) / 2
    eta <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    eta[adjacent_rank] - eta[adjacent_rank + 1L]
  }, numeric(1))
}

run_fpca_repeated_root_qa <- function(replications = 150L, n = 400L,
                                      equality_draws = 100000L,
                                      seed = 761L) {
  p_values <- with_fixed_seed(seed, vapply(seq_len(replications), function(iteration) {
    fit <- gram_pca(matrix(stats::rnorm(n * 2L), n, 2L), vectors = TRUE)
    test <- fpca_pair_inference(
      fit, 1L, equality_draws = equality_draws, equality_seed = 82471L
    )
    c(general = test$fpca_p_equal_general,
      anderson = test$fpca_anderson_p_equal)
  }, numeric(2L)))
  rejection <- rowMeans(p_values < 0.05)
  agreement <- stats::cor(p_values["general", ], p_values["anderson", ])
  stopifnot(
    abs(rejection[["general"]] - 0.05) < 0.04,
    abs(rejection[["anderson"]] - 0.05) < 0.04,
    agreement > 0.90
  )
  data.frame(
    replications = replications, n = n,
    fpca_general_type1 = rejection[["general"]],
    fpca_anderson_type1 = rejection[["anderson"]],
    p_value_correlation = agreement
  )
}

run_eigengap_smoke_checks <- function() {
  probabilities <- c(0.25, 0.50, 0.75)
  half_spacing <- simulate_goe_spacing(2L, 1L, 20000L, seed = 7301L) / 2
  empirical <- unname(quantile(half_spacing, probabilities))
  rayleigh <- sqrt(-2 * log(1 - probabilities))
  stopifnot(max(abs(empirical - rayleigh)) < 0.035)

  statistic_grid <- c(0.1, 1, 2, 5, 10)
  stopifnot(max(abs(
    pchisq(statistic_grid, df = 2, lower.tail = FALSE) - exp(-statistic_grid / 2)
  )) < 1e-14)

  central_acceptance <- function(x) {
    qchisq(0.025, df = 2) <= x && x <= qchisq(0.975, df = 2)
  }
  for (x in c(0.01, 0.1, 1, 3, 8, 12)) {
    interval <- invert_noncentral_chisq_gap_ci(x, 100, 0.95)
    stopifnot(identical(isTRUE(interval$ci_contains_zero), central_acceptance(x)))
  }

  # A deterministic large-n repeated-outlier spectrum checks the first-order
  # equivalence between correction before gap construction and correct
  # high-complexity calibration of the raw sample gap.
  n <- 1000000L
  spectrum <- c(6 + 1 / sqrt(n), 6 - 1 / sqrt(n), rep(1, 2000L))
  check <- multiplicity2_gap_inference(spectrum, 1L, n, 2L, 0.95)
  stopifnot(abs(check$p_sample_HC - check$p_proposed) < 0.01)

  # The local-alternative parameterization has an exact symmetric-gap identity.
  for (n_value in c(150, 300, 600)) for (delta in c(0, 0.5, 1, 2, 3)) {
    upper <- 6 * (1 + delta / (2 * sqrt(n_value)))
    lower <- 6 * (1 - delta / (2 * sqrt(n_value)))
    stopifnot(abs(symmetric_relative_gap(upper, lower) - delta / sqrt(n_value)) < 1e-12)
  }

  # Tiny observed gaps may make the literal equal-tail set empty, but the
  # boundary-adjusted magnitude interval and upper confidence bound stay finite.
  boundary <- multiplicity2_gap_inference(c(6, 6, rep(1, 100)), 1L, 1000L, 2L)
  stopifnot(
    isTRUE(boundary$proposed_gap_ci_raw_empty),
    all(is.finite(c(boundary$proposed_gap_upper95,
                    boundary$proposed_gap_ci_boundary_lower,
                    boundary$proposed_gap_ci_boundary_upper)))
  )

  # Classical empirical-score normalization and Gaussian variance QA.
  gaussian_fit <- with_fixed_seed(991L, gram_pca(matrix(rnorm(600L * 2L), 600L, 2L)))
  gaussian_interval <- fpca_eigenvalue_interval(gaussian_fit, 1L)
  gaussian_reference <- gaussian_fit$values[1] * sqrt(2 / gaussian_fit$n)
  stopifnot(
    gaussian_interval$fpca_score_identity_error < 1e-10,
    abs(gaussian_interval$fpca_alpha_se / gaussian_reference - 1) < 0.10
  )

  # Classical functions must remain independent of high-complexity quantities.
  classical_source <- paste(
    deparse(body(fpca_eigenvalue_interval)),
    deparse(body(fpca_pair_inference)), collapse = "\n"
  )
  stopifnot(!grepl("Delta|r2", classical_source))
  run_fpca_repeated_root_qa()
  invisible(TRUE)
}

dimension_free_inference <- function(sample_eigenvalues, target_index, K0, n,
                                     confidence_level = DEFAULT_CONFIDENCE_LEVEL,
                                     certification_threshold = 0.20) {
  m <- deflated_spectral_moments(sample_eigenvalues, target_index, K0, n, 1:4)
  lambda <- sort(as.numeric(sample_eigenvalues), decreasing = TRUE)[target_index]
  alpha <- -1 / m[["m1"]]
  Delta <- m[["m1"]]^2 / m[["m2"]]
  r2 <- -m[["m1"]] / (lambda * m[["m2"]])

  d_alpha <- 1 / Delta
  d_Delta <- 2 * m[["m1"]] - 2 * m[["m1"]]^2 * m[["m3"]] / m[["m2"]]^2
  d_r2 <- r2 * (m[["m2"]] / m[["m1"]] - 1 / lambda - 2 * m[["m3"]] / m[["m2"]])

  se_alpha <- if (Delta > 0) sqrt(2 * alpha^2 / (n * Delta)) else NA_real_
  se_Delta <- sqrt(2 * d_Delta^2 / (n * m[["m2"]]))
  se_r2 <- sqrt(2 * d_r2^2 / (n * m[["m2"]]))
  z_two <- qnorm(1 - (1 - confidence_level) / 2)
  z_one <- qnorm(confidence_level)

  data.frame(
    spike_index = target_index,
    K0 = K0,
    hat_lambda = lambda,
    mhat1 = m[["m1"]], mhat2 = m[["m2"]], mhat3 = m[["m3"]], mhat4 = m[["m4"]],
    hat_alpha = alpha, hat_Delta = Delta, hat_r2 = r2,
    d_alpha = d_alpha, d_Delta = d_Delta, d_r2 = d_r2,
    se_alpha = se_alpha, se_Delta = se_Delta, se_r2 = se_r2,
    ci_alpha_lower = alpha - z_two * se_alpha,
    ci_alpha_upper = alpha + z_two * se_alpha,
    ci_Delta_lower = Delta - z_two * se_Delta,
    ci_Delta_upper = Delta + z_two * se_Delta,
    ci_r2_lower = r2 - z_two * se_r2,
    ci_r2_upper = r2 + z_two * se_r2,
    ci_r2_display_lower = max(0, r2 - z_two * se_r2),
    ci_r2_display_upper = min(1, r2 + z_two * se_r2),
    phase_lower_bound = Delta - z_one * se_Delta,
    phase_certified = as.integer(Delta - z_one * se_Delta > certification_threshold)
  )
}

attach_truth_and_coverage <- function(estimate, truth) {
  out <- cbind(estimate, truth[setdiff(names(truth), "spike_index")])
  out$population_phase_status <- ifelse(out$true_Delta > 0, "supercritical", "subcritical")
  out$regular_inference_valid <- TRUE
  out$cover_alpha <- as.integer(out$ci_alpha_lower <= out$true_alpha & out$true_alpha <= out$ci_alpha_upper)
  out$cover_Delta <- as.integer(out$ci_Delta_lower <= out$true_Delta & out$true_Delta <= out$ci_Delta_upper)
  out$cover_r2 <- as.integer(out$ci_r2_lower <= out$true_r2 & out$true_r2 <= out$ci_r2_upper)
  out
}

cluster_reliability <- function(sample_eigenvalues, cluster_indices, K0, n) {
  estimates <- lapply(cluster_indices, function(j) {
    dimension_free_inference(sample_eigenvalues, j, K0, n)
  })
  mean(vapply(estimates, function(x) x$hat_r2, numeric(1)))
}

population_derivatives <- function(population_eigenvalues, M, n, spike_index) {
  truth <- population_spike_truth(population_eigenvalues, M, n, spike_index)
  bulk <- sort(population_eigenvalues, decreasing = TRUE)[-(seq_len(M))]
  alpha <- truth$true_alpha
  psi <- truth$true_psi
  Delta <- truth$true_Delta
  phi_second <- 2 * sum(bulk^2 / (alpha - bulk)^3) / n
  d_alpha <- 1 / truth$true_Delta
  d_Delta <- phi_second / Delta
  dr2_dalpha <- Delta / psi + alpha * phi_second / psi - alpha * Delta^2 / psi^2
  d_r2 <- dr2_dalpha / Delta
  c(d_alpha = d_alpha, d_Delta = d_Delta, d_r2 = d_r2)
}

orient_projection <- function(projection, anchor_projection) {
  ifelse(anchor_projection < 0, -projection, projection)
}
