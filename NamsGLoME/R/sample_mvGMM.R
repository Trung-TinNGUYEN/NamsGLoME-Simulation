#' @export
sample_mvGMM = function(gmm, n) {
  # Create an n-sample from a multivariate Gaussian mixture model

  p <- ncol(gmm$means)
  X <- matrix(data = 0, nrow = n, ncol = p)

  Z <- stats::rmultinom(n = n, size = 1, prob = gmm$weights)
  klas <- apply(X = Z, MARGIN = 2, FUN = which.max)
  for (i in 1:n) {
    X[i, ] <- MASS::mvrnorm(mu = gmm$means[klas[i], ], Sigma = gmm$covariances[, , klas[i]])
  }

  return(list(X = X, klas = klas))
}
