#' @export
solve_weighted_LASSOexpert = function(X, y, lambda, weights, beta0, beta, sigma2, verbose = FALSE) {
  # fits a weighted Gaussian regression problem with lasso-regularized maximum likelihood by using a
  # coordinate ascent algorithm

  iter <- 0
  converge <- FALSE
  max_iter <- 300
  threshold <- 1e-7

  n <- nrow(X)
  p <- ncol(X)

  tau <- weights

  # Objective function for the initial model
  log_Phi_y <- -0.5 * log(2 * pi) - 0.5 * log(sigma2) - 0.5 * ((y - beta0 * ones(n, 1) - X %*% beta) ^ 2) / sigma2
  Q_old <- sum(tau * log_Phi_y) - lambda * sum(abs(beta))
  while (!converge && (iter < max_iter)) {
    for (j in 1:p) {
      Xj <- X[, j]
      ## Expert Network parameters update
      # coordinate ascent for lasso to update the reg coefficients Beta's
      Rkj <- y - beta0 * ones(n, 1) -  X %*% beta + beta[j] * Xj
      u <- (t(tau * Xj) %*% Xj) ^ (-1) * t(tau * Xj) %*% Rkj
      eta <- lambda * (t(tau * Xj) %*% Xj) ^ (-1) * sigma2
      beta[j] <- sign(u) * (max(abs(u) - eta, 0)) # Soft-thresholding operator
    }

    beta0 <- sum(tau * (y - X %*% beta)) / sum(tau) # Update the intercept
    sigma2 <- sum(tau * ((y - beta0 * ones(n, 1) - X %*% beta) ^ 2)) / sum(tau)  + .Machine$double.eps# Update the regressor variance

    # Objective function
    log_Phi_y <- -0.5 * log(2 * pi) - 0.5 * log(sigma2) - 0.5 * ((y - beta0 * ones(n, 1) - X %*% beta) ^ 2) / sigma2
    Q <- sum(tau * log_Phi_y) - lambda * sum(abs(beta))

    # convergence test
    converge  <- abs(Q - Q_old) <= threshold
    if (verbose) {
      message("Coordinate Ascent for Experts Net: Iteration: ", iter, " | Q_gating:: "  , Q_old)
    }
    Q_old <- Q
    iter <- iter + 1
  }

  return(list(beta = beta, sigma2 = sigma2, beta0 = beta0))
}
