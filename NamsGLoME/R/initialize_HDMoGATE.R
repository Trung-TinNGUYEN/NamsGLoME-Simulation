#' @export
initialize_HDMoGATE = function(X, y, K) {
  # Initializes a mixture of high-dimensional gaussian-gated mixture-of-experts by the EM-Lasso algorithm

  n <- nrow(X)
  p <- ncol(X)

  # Gating Net parameters
  outkmeans <- stats::kmeans(x = X, centers = K)

  klas <- outkmeans$cluster
  mu <- outkmeans$centers

  R <- array(data = 0, dim = c(p, p, K))
  alpha <- matrix(data = 0, nrow = 1, ncol = K)

  for (k in 1:K) {
    Xk <- X[klas == k, ]
    nk <- nrow(Xk)
    alpha[k] <- nk / n
    z <- (Xk - ones(nk, 1) %*% mu[k, ])
    R[, , k] <- diag(apply(X = z ^ 2, MARGIN = 2, FUN = sum)) / nk
  }
  #     nu2=zeros(K,p);
  #     for k=1:K
  #         nu2(k,:) = diag(R(:,:,k));
  #         R(:,:,k) = diag(nu2(k,:));
  #     end
  # Expert Net parameters
  #klas = randi(K,n,1);
  beta0 <- matrix(data = 0, nrow = 1, ncol = K)
  beta <- matrix(data = 0, nrow = p, ncol = K)
  sigma2 <- matrix(data = 0, nrow = 1, ncol = K)

  for (k in 1:K) {
    Xk <- X[klas == k, ]
    yk <- y[klas == k]
    nk <- length(yk)

    # The regression coefficients
    beta[, k] <- solve(t(Xk) %*% Xk, tol = 0) %*% t(Xk) %*% yk

    beta0[k] <- sum(yk - Xk %*% beta[, k])

    # the variances sigma2k
    sigma2[k] <- sum((yk - beta0[k] * ones(nk, 1) - Xk %*% beta[, k]) ^ 2) / nk
  }

  return(list(alpha = alpha, mu = mu, R = R, beta0 = beta0, beta = beta, sigma2 = sigma2))
}
