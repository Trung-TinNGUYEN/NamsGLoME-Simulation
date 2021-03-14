#' @export
sample_MoGATE = function(alpha, mu, R, beta0, beta, sigma, n) {
  # Draws an n-sample from a mixture of gaussian-gated experts model

  K <- length(alpha)

  # Sample the predictors
  gmm <- list()
  gmm$weights <- alpha
  gmm$means <- mu
  gmm$covariances <- R
  mvGMM <- sample_mvGMM(gmm, n)
  X <- mvGMM$X
  klasx <- mvGMM$klas

  # Sample the responses

  # Calculate the gating net probabilites
  gatingProb <- posterior_mvGMM(X, gmm)$tau
  y <- zeros(n, 1)
  Z <- zeros(n, K)
  klasy <- zeros(K, 1)
  for (i in 1:n) {
    Zik <- stats::rmultinom(n = 1, size = 1, prob = gatingProb[i,])
    Z[i,] <- Zik
    zi <- which(Zik == 1)
    klasy[i] <- zi

    y[i] <- beta0[zi] + X[i, ] %*% beta[, zi] + stats::rnorm(n = 1, mean = 0, sd = sigma[zi])
  }

  # # Statistics (means, variances)
  # # E[yi|zi=k]
  # Ey_k = Xbeta*beta;
  # # E[yi]
  # Ey = sum(Piik.*Ey_k,2);
  # # Var[yi|zi=k]
  # Vary_k = Sigmak.^2;
  # # Var[yi]
  # Vary = sum(Piik.*(Ey_k.^2 + ones(n,1)*Vary_k),2) - Ey.^2;
  #
  # stats.Ey_k = Ey_k;
  # stats.Ey = Ey;
  # stats.Vary_k = Vary_k;
  # stats.Vary = Vary;

  stats <- list()
  stats$klasx <- klasx
  stats$klasy <- klasy

  return(list(X = X, y = y, stats = stats))
}
