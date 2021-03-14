#' @export
solve_weighted_LASSOgating = function(X, gamma, weights, alpha, mu, nu2, verbose = FALSE) {
  # fits a weighted multivariate Gaussian problem with lasso-regularized maximum likelihood by using a
  # coordinate ascent algorithm. The regularization concerns the mean vector and the cov. matrix is
  # supposed to be diagonal.

  iter <- 0
  converge <- FALSE
  max_iter <- 300
  threshold <- 1e-7

  n <- nrow(X)
  p <- ncol(X)

  tau <- weights

  # Objective function for the initial model
  R <- diag(nu2)
  Q_old <- sum(tau * log(alpha)) + sum(tau * mvgaussian_pdf(X, mu, R, "diagonal")$log_fxi) - gamma * sum(abs(mu))
  while (!converge && (iter < max_iter)) {
    for (j in 1:p) {
      Xj <- X[, j]
      ## Gating Network parameters update
      # coordinate ascent for lasso to updata Mu's
      u <- sum(tau * Xj) / sum(tau)
      eta <- gamma * nu2[j] / sum(tau)
      mu[j] <- sign(u) * (max(abs(u) - eta, 0)) # Soft-thresholding operator
      # the covariance matrix elements
      z <- (Xj - ones(n, 1) * mu[j]) * sqrt(tau)
      nu2[j] <- (t(z) %*% z) / sum(tau) + .Machine$double.eps
    }
    R <- diag(nu2) # cov matrix update

    # Objective function
    Q <- sum(tau * log(alpha)) + sum(tau * mvgaussian_pdf(X, mu, R, "diagonal")$log_fxi) - gamma * sum(abs(mu))

    # Convergence test
    converge  <- abs(Q - Q_old) <= threshold
    if (verbose) {
      message("Coordinate Ascent for Gating Net: Iteration: ", iter, " | Q_gating:: ", Q_old)
    }
    Q_old <- Q
    iter <- iter + 1
  }

  return(list(mu = mu, nu2 = nu2, R = R))
}
