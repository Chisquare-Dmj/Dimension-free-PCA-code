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
  bulk <- values[(K0 + 1L):length(values)]
  gaps <- bulk - z
  if (any(abs(gaps) <= .Machine$double.eps * max(1, abs(z)))) {
    stop("A deflated bulk eigenvalue is numerically tied with the target sample eigenvalue.")
  }
  moments <- vapply(orders, function(s) sum(gaps^(-s)) / n, numeric(1))
  names(moments) <- paste0("m", orders)
  moments
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
