# Population models and observation generators used by all experiments.

calibrate_population_spike <- function(population_values, target_index, signal_count,
                                       sample_size, target_Delta) {
  bulk_values <- population_values[-seq_len(signal_count)]
  phase_equation <- function(alpha) {
    1 - sum((bulk_values / (alpha - bulk_values))^2) / sample_size - target_Delta
  }
  lower <- max(bulk_values) * (1 + 1e-10)
  upper <- max(population_values[target_index], lower + 1)
  while (phase_equation(upper) < 0) upper <- upper * 2
  uniroot(phase_equation, c(lower, upper), tol = 1e-12)$root
}

compound_symmetry_block <- function(size = 10L, correlation = 0.55) {
  correlation * matrix(1, size, size) + (1 - correlation) * diag(size)
}

build_panel_population <- function(scenario, p = 20L, T = 30L, seed,
                                   sample_size = 300L,
                                   p2_leading_multipliers = c(5, 4),
                                   target_Delta3 = -0.10) {
  scenario <- match.arg(scenario, c("P1_independent", "P2_block_AR1"))
  d <- p * T
  if (scenario == "P1_independent") {
    population <- with_fixed_seed(seed, {
      A <- diag(runif(p, 0.5, 1.5))
      B <- diag(runif(T, 0.5, 1.5))
      bulk <- kronecker(B, A)
      directions <- fixed_orthonormal_directions(d, 3L, seed + 1L)
      theta <- c(6, 4, 2)
      initial_sigma <- bulk + directions %*% diag(theta, 3L) %*% t(directions)
      initial_eig <- eigen(initial_sigma, symmetric = TRUE)
      final_values <- initial_eig$values
      original_alpha3 <- final_values[3]
      final_values[3] <- calibrate_population_spike(
        final_values, 3L, 3L, sample_size, target_Delta3
      )
      sigma <- sweep(initial_eig$vectors, 2L, final_values, "*") %*% t(initial_eig$vectors)
      list(
        A = A, B = B, bulk = bulk, directions = initial_eig$vectors[, 1:3],
        theta = theta, sigma = sigma,
        calibrated_multiplier3 = final_values[3] / original_alpha3
      )
    })
  } else {
    A <- outer(seq_len(p), seq_len(p), function(i, j) 0.6^abs(i - j))
    block <- compound_symmetry_block(10L, 0.55)
    if (T %% 10L != 0L) stop("P2_block_AR1 requires T to be divisible by 10.")
    B <- block_diagonal(replicate(T / 10, block, simplify = FALSE))
    bulk <- kronecker(B, A)
    eig <- eigen(bulk, symmetric = TRUE)
    final_values <- eig$values
    final_values[1:2] <- final_values[1:2] * p2_leading_multipliers
    original_alpha3 <- final_values[3]
    final_values[3] <- calibrate_population_spike(
      final_values, 3L, 3L, sample_size, target_Delta3
    )
    sigma <- sweep(eig$vectors, 2L, final_values, "*") %*% t(eig$vectors)
    population <- list(
      A = A, B = B, bulk = bulk, directions = eig$vectors[, 1:3],
      theta = final_values[1:3], sigma = sigma,
      calibrated_multiplier3 = final_values[3] / original_alpha3
    )
  }
  eig_sigma <- eigen(population$sigma, symmetric = TRUE)
  population$values <- eig_sigma$values
  population$vectors <- eig_sigma$vectors
  population$sqrt <- sweep(eig_sigma$vectors, 2L, sqrt(pmax(eig_sigma$values, 0)), "*") %*% t(eig_sigma$vectors)
  population$scenario <- scenario
  population$p <- p
  population$T <- T
  population$dimension <- d
  population
}

sample_panel_observations <- function(n, population, distribution = "Gaussian") {
  generate_standard_scores(distribution, n, population$dimension) %*% population$sqrt
}

build_vector_calibration_population <- function(n, p, Delta) {
  c_ratio <- (p - 1) / n
  alpha <- 1 + sqrt(c_ratio / (1 - Delta))
  values <- c(alpha, rep(1, p - 1L))
  list(
    values = values,
    sqrt_values = sqrt(values),
    alpha = alpha,
    psi = alpha * (1 + c_ratio / (alpha - 1)),
    Delta = Delta,
    r2 = alpha / (alpha * (1 + c_ratio / (alpha - 1))) * Delta,
    n = n, p = p
  )
}

sample_vector_observations <- function(n, population, distribution = "Gaussian") {
  sweep(generate_standard_scores(distribution, n, population$p), 2L, population$sqrt_values, "*")
}

build_cross_covariance <- function(scenario, p, seed) {
  scenario <- match.arg(scenario, c("F1_diagonal", "F2_block", "F3_Haar"))
  if (scenario == "F1_diagonal") {
    values <- c(6, 4, 3, rep(1, p - 3L))
    vectors <- diag(p)
  } else if (scenario == "F2_block") {
    block <- compound_symmetry_block(10L, 0.55)
    if (p %% 10L != 0L) stop("F2_block requires p to be divisible by 10.")
    base <- block_diagonal(replicate(p / 10, block, simplify = FALSE))
    eig <- eigen(base, symmetric = TRUE)
    values <- eig$values
    values[1:3] <- values[1:3] * c(5, 4, 3)
    vectors <- eig$vectors
  } else {
    values <- c(6, 4, 3, rep(1, p - 3L))
    vectors <- fixed_haar_matrix(p, seed)
  }
  list(
    scenario = scenario,
    values = values,
    vectors = vectors,
    sqrt = sweep(vectors, 2L, sqrt(values), "*") %*% t(vectors)
  )
}

build_functional_population <- function(scenario, p = 100L, N = 200L, decay = 1.1, seed = 201L) {
  cross <- build_cross_covariance(scenario, p, seed)
  mu <- seq_len(N)^(-decay)
  spectrum_grid <- expand.grid(coordinate = seq_len(p), basis = seq_len(N))
  spectrum_grid$value <- cross$values[spectrum_grid$coordinate] * mu[spectrum_grid$basis]
  spectrum_grid <- spectrum_grid[order(spectrum_grid$value, decreasing = TRUE), ]
  rownames(spectrum_grid) <- NULL
  list(
    scenario = scenario, p = p, N = N, decay = decay,
    cross = cross, mu = mu,
    values = spectrum_grid$value,
    spectrum_map = spectrum_grid
  )
}

sample_functional_coefficients <- function(n, population, distribution = "Gaussian") {
  # The returned n x (p*N) matrix contains exact orthonormal-basis coefficients.
  # Columns are ordered by basis block: all p coordinates for basis 1, then basis 2, etc.
  blocks <- vector("list", population$N)
  for (k in seq_len(population$N)) {
    scores <- generate_standard_scores(distribution, n, population$p)
    blocks[[k]] <- (scores %*% population$cross$sqrt) * sqrt(population$mu[k])
  }
  do.call(cbind, blocks)
}

functional_population_direction <- function(population, spike_index) {
  map <- population$spectrum_map[spike_index, ]
  direction <- numeric(population$p * population$N)
  block_start <- (map$basis - 1L) * population$p
  direction[block_start + seq_len(population$p)] <- population$cross$vectors[, map$coordinate]
  direction
}

large_domain_spectrum <- function(n, repeated = FALSE, N = 2L * n, T = n / 3, decay = 1.1) {
  T <- as.integer(T)
  N <- as.integer(N)
  values <- numeric(N)
  if (repeated) {
    values[1:3] <- c(6, 6, 3)
  } else {
    values[1:3] <- c(6, 4, 3)
  }
  if (N >= 4L) values[4:N] <- (floor((seq.int(4L, N) - 4L) / T) + 1)^(-decay)
  values
}

sample_large_domain_coefficients <- function(n, eigenvalues, distribution = "Gaussian") {
  sweep(generate_standard_scores(distribution, n, length(eigenvalues)), 2L, sqrt(eigenvalues), "*")
}

large_domain_target <- function(N) {
  f <- numeric(N)
  f[c(1, 2, 13)] <- c(0.6, 0.3, sqrt(0.55))
  f
}

large_domain_basis_matrix <- function(N, T, m) {
  grid <- seq_len(m) * T / m
  sqrt(2 / T) * cos(outer(seq_len(N), grid, function(j, t) j * pi * t / T))
}

coefficient_to_grid <- function(coefficients, T, m) {
  basis <- large_domain_basis_matrix(ncol(coefficients), T, m)
  weight <- T / m
  coefficients %*% basis * sqrt(weight)
}

direction_to_grid <- function(direction, T, m) {
  basis <- large_domain_basis_matrix(length(direction), T, m)
  as.numeric(direction %*% basis) * sqrt(T / m)
}
