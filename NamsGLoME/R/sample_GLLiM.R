#' @export
sample_GLLiM = function(alpha, mu, R, beta0, beta, sigma, n, ny) {
  # Output: Draws an n-sample from a mixture of gaussian-gated experts model with linear mean experts (GLLiM).
    # X (nxp), y (nx1): sample data set (X_i,y_i)_{i=1,...n}
    # Xny(n*nyxp),yny(n*nyx1): repeat n_y times the sample X for Monte Carlo method
    #                         then combine by rows
    # stats$klasx <- klasx: label for laten variable used to sample X
    # stats$klasy <- klasy: label for laten variable used to sample y
    # gllim_pdf (n*nyx1): value of  pdf of gllim model at (xi,yi), i.e, f(yi|xi, Psi)


  # Input:
    # alpha: prior probabilities for the weights of Gaussian gating functions
    # mu: mean vectors of Gaussian gating functions
    # R: Covariance matrices of Gaussian gating functions
    # beta0, beta: coefficients of the linear means of Gaussian experts.
    # sigma: Covariance matrices of Gaussian experts
    # n: sample size
    # ny: repeat n_y times the sample X for Monte Carlo method then combine by rows

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
  gatingProb <- posterior_mvGMM(X, gmm)$tau # nxK


  y <- zeros(n, 1)
  Z <- zeros(n, K)
  klasy <- zeros(K, 1)

  for (i in 1:n) {
    Zik <- stats::rmultinom(n = 1, size = 1, prob = gatingProb[i,]) # Kx1
    Z[i,] <- Zik
    zi <- which(Zik == 1)
    klasy[i] <- zi

    y[i] <- beta0[zi] + X[i, ] %*% beta[, zi] + stats::rnorm(n = 1, mean = 0, sd = sigma[zi])
  }

  # Repeat n_y times the sample X for Monte Carlo method then combine by rows
  Xny <- do.call(rbind, replicate(ny, X, simplify=FALSE))
  klasx_ny <- do.call(rbind, replicate(ny, klasx, simplify=FALSE))

  gatingProb_ny <- posterior_mvGMM(Xny, gmm)$tau # (n*ny)xK

  yny <- zeros(n*ny, 1)
  Zny <- zeros(n*ny, K)
  klasy_ny <- zeros(K, 1)

  gllim_pYxz <- matrix(0,n*ny,K)

  for (i in 1:(n*ny)) {
    Zik_ny <- stats::rmultinom(n = 1, size = 1, prob = gatingProb_ny[i,]) # Kx1
    Zny[i,] <- Zik_ny
    zi_ny <- which(Zik_ny == 1)
    klasy_ny[i] <- zi_ny

    yny[i] <- beta0[zi_ny] + Xny[i, ] %*% beta[, zi_ny] + stats::rnorm(n = 1, mean = 0, sd = sigma[zi_ny])
  }

  for (k in 1:K){
    gllim_pYxz[,k] <- exp(loggausspdf(t(yny),t(matrix(beta0[k] + Xny%*% beta[, k],ncol=1)),
                                      matrix(sigma[k]^2,ncol = 1))) # (n*ny)xK
  }
  # Summation of log-likelihood value over all sample.
  gllim_pdf <- sum(log(rowSums(gllim_pYxz*gatingProb_ny)))

  # Vector of pdf values over all sample s_0(Y_i|X_i), i = 1,..., n*ny.
  gllim_pdf_vec <- rowSums(gllim_pYxz*gatingProb_ny)

  # Summary statistics
  stats <- list()
  stats$klasx <- klasx
  stats$klasy <- klasy

  stats$klasx_ny <- klasx_ny
  stats$klasy_ny <- klasy_ny


  return(list(X = X, y = y, Xny = Xny, yny = yny, stats = stats,
              gllim_pdf = gllim_pdf,gllim_pdf_vec = gllim_pdf_vec))
}
