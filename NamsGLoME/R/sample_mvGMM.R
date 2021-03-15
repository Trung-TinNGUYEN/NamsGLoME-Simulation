#' @export
sample_mvGMM = function(gmm, n) {
  # create an n-sample from (truncated) Gaussian mixture model.

  p <- ncol(gmm$means)
  X <- matrix(data = 0, nrow = n, ncol = p)

  Z <- stats::rmultinom(n = n, size = 1, prob = gmm$weights)
  klas <- apply(X = Z, MARGIN = 2, FUN = which.max)
  for (i in 1:n) {
    # Truncated Gaussian
    #X[i, ] <- truncnorm::rtruncn   orm(1, a = 0, b = 1, mean = gmm$means[klas[i], ], sd = gmm$covariances[, , klas[i]])
    #X[i, ] <- stats::rnorm(1, mean = gmm$means[klas[i], ], sd = gmm$covariances[, , klas[i]])
    X[i, ] <- MASS::mvrnorm(mu = gmm$means[klas[i], ], Sigma = gmm$covariances[, , klas[i]])
  }

  return(list(X = X, klas = klas))
}
