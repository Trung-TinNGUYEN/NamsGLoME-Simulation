#' @export
mvgaussian_pdf = function(X, mu, sigma, covtype = c("full", "diagonal")) {
  ###########################################################################
  #  Clacul d'une densite gaussienne  pour un echantillon donn??
  # Entrees:
  #
  #           X    : Tableau d'individus [nxd]
  #           mu [1xd]; centre
  #           sigma # mat de variance [dxd]
  # Sorties:
  #
  #            fxi : (2*pi*|sigma|)^(-.5*d)*exp{-0.5*(xi - mu)
  #            sigma^-1(xi-mu)'}  pour i=1,...n. vecteur de dim n
  #            log_fxi:  log(fxi) pour i=1,...n. vecteur de dim n pour
  #
  # #########################################################################

  covtype <- match.arg(covtype)

  sigma <- as.matrix(sigma) # Make sure sigma(1x1) is still a matrix for diag

  if (covtype == "diagonal") {
    diagelements <- diag(sigma) + .Machine$double.eps
    invSigma <- diag(as.matrix(1 / diagelements))
    detSigma <- prod(diagelements)
  } else if (covtype == "full") {
    detSigma <- det(sigma)
    invSigma <- solve(sigma, tol = 0)
  }

  n <- nrow(X)
  d <- ncol(X)

  z <- ((X - ones(n, 1) %*% mu) %*% invSigma) * (X - ones(n, 1) %*% mu)
  mahalanobis <- apply(z, 1, sum)
  log_fxi <- -(d / 2) * log(2 * pi) - 0.5 * log(detSigma) - 0.5 * mahalanobis

  denom <- (2 * pi) ^ (d / 2) * (detSigma) ^ (1 / 2)
  fxi <-  exp(-0.5 * mahalanobis) / denom

  return(list(log_fxi = log_fxi, fxi = fxi))
}
