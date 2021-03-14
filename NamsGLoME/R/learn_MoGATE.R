#' @export
learn_MoGATE = function(X, y, K, nb_EM_runs = 1, verbose = FALSE) {
  # Fits a mixture of gaussian-gated mixture-of-experts by the EM algorithm

  n <- nrow(X)
  p <- ncol(X)

  if (dim(y)[1] != n) {
    y <- t(y)
  }

  max_iter_EM <- 1500
  threshold <- 1e-6

  loglik <- -Inf
  best_loglik <- -Inf
  stored_cputime <- c()
  EM_try <- 1

  while (EM_try <= nb_EM_runs) {
    if (nb_EM_runs > 1 && verbose) {
      message("EM try number: ", EM_try, "\n")
    }

    time <- Sys.time()
    # EM Initialisation
    initialization <- initialize_HDMoGATE(X, y, K)
    alpha <- initialization$alpha
    mu <- initialization$mu
    R <- initialization$R
    beta0 <- initialization$beta0
    beta <- initialization$beta
    sigma2 <- initialization$sigma2
    nu2 <- matrix(data = 0, nrow = K, ncol = p)
    for (k in 1:K) {
      nu2[k,] <- diag(R[, , k])
    }

    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf
    stored_loglik <- c()

    ## EM ####
    while (!converge && (iter < max_iter_EM)) {
      iter <- iter + 1

      ## E-Step
      log_alpha_Phi_xy <- matrix(data = 0, nrow = n, ncol = K)
      for (k in 1:K) {
        # Gating network conditional density
        log_Phi_x <- mvgaussian_pdf(X, mu[k,], R[, , k], covtype = "diagonal")$log_fxi

        # Expert Network conditional density
        log_Phi_y <- -0.5 * log(2 * pi) - 0.5 * log(sigma2[k]) - 0.5 * ((y - beta0[k] * ones(n, 1) - X %*% beta[, k]) ^ 2) / sigma2[k]
        #log_Phi_y <-  mvgaussian_pdf(y, beta0(k)*ones(n,1) + X*beta[, k], sigma2(k));

        # HDMoGATE conditional density
        log_alpha_Phi_xy[, k] <- log(alpha[k]) * ones(n, 1) + log_Phi_x + log_Phi_y

      }
      #log_tau  <- log_normalize(log_piik_fik);
      log_sum_alpha_PhixPhiy <- logsumexp(log_alpha_Phi_xy, 1)
      log_tau <- log_alpha_Phi_xy - log_sum_alpha_PhixPhiy %*% ones(1, K)

      tau <- exp(log_tau)
      tau <- tau / (apply(X = tau, MARGIN = 1, FUN = sum) %*% ones(1, K))

      ## M-Step
      for (k in 1:K) {
        tauk <- tau[, k]

        # Gating Network

        # The mixing proportions
        alpha[k] <- sum(tauk) / n

        # The Gaussian means
        mu[k,] <- apply(X = X * (tauk %*% ones(1, p)), MARGIN = 2, sum) / sum(tauk)
        z <- (X - ones(n, 1) %*% mu[k,]) * (sqrt(tauk) %*% ones(1, p))

        # The Gaussian cov matrices
        # R[, , k] <- (t(z) %*% z) / sum(tauk)
        R[, , k] <- diag(apply(X = z^2, MARGIN = 2, FUN = sum) / sum(tauk)) + .Machine$double.eps
        nu2[k, ] <- diag(R[, , k])

        # Expert Network
        Xk <- X * (sqrt(tauk %*% ones(1, p)))

        yk <- y * sqrt(tauk)

        # Update the regression coefficients
        beta[, k] <- solve(t(Xk) %*% Xk, tol = 0) %*% t(Xk) %*% yk
        beta0[k] <- sum(tauk * (y - X %*% beta[, k])) / sum(tauk)

        # Update the variances sigma2
        sigma2[k] <- sum(tauk * ((y - beta0[k] * ones(n, 1) - X %*% beta[, k]) ^ 2)) / sum(tauk) + .Machine$double.eps
      }

      # Observed-data log-likelihood
      loglik <- sum(log_sum_alpha_PhixPhiy)

      if (verbose) {
        message("EM - MoE: Iteration: ", iter, " | log-likelihood: "  , loglik)
      }

      converge <- abs((loglik - prev_loglik) / prev_loglik) <= threshold
      if (is.na(converge)) {
        converge <- FALSE
      } # Basically for the first iteration when prev_loglik is -inf
      prev_loglik <- loglik
      stored_loglik <- c(stored_loglik, loglik)
    }
    EM_try <- EM_try + 1
    stored_cputime <- c(stored_cputime, Sys.time() - time)

    # Results

    # Model params
    param <- list()
    param$alpha <- alpha
    param$mu <- mu
    param$R <- R
    param$nu2 <- nu2
    param$beta0 <- beta0
    param$beta <- beta
    param$sigma2 <- sigma2

    solution <- list()
    solution$param <- param

    # Model stats

    # Piik <- Piik(1:m,:);#To Be UPDATED
    # tau <- tau(1:m,:);
    # solution.Piik <- Piik;
    solution$stats <- list()
    solution$stats$tau <- tau
    solution$stats$log_Phi_x <- log_Phi_x
    solution$stats$log_Phi_y <- log_Phi_y
    solution$stats$log_alpha_Phi_xy <- log_alpha_Phi_xy
    solution$stats$ml <- loglik
    solution$stats$stored_loglik <- stored_loglik

    # Parameter vector of the estimated model
    # rk <- c()# zeros(K*p*(p+1)/2, 1);
    # for (k in 1:K) {
    #   rk <- c(rk, t(apply(lower.tri(R[, , K], diag = TRUE), 1, function(x) x[x != 0])))
    # }

    # Parameter vector of the estimated model
    psi <- c(as.vector(alpha), as.vector(mu), as.vector(nu2), as.vector(beta0), as.vector(beta), as.vector(sigma2))
    solution$stats$psi <- psi

    # Classsification pour EM : MAP(piik) (cas particulier ici to ensure a convex segmentation of the curve(s).
    outMAP <- MAP(tau)
    klas <- outMAP$klas
    Zik <- outMAP$Z
    solution$stats$klas <- klas

    # Statistics (means, variances) To Be UPDATED
    #         # E[yi|zi=k]
    #         Ey_k <- X(1:m,:)*Beta;
    #         solution.Ey_k <- Ey_k;
    #         # E[yi]
    #         Ey <- sum(Piik.*Ey_k,2);
    #         solution.Ey <- Ey;
    #
    #         # Var[yi|zi=k]
    #         Vary_k <- sigma2;
    #         solution.Vary_k <- Vary_k;
    #
    #         # Var[yi]
    #         Vary <- sum(Piik.*(Ey_k.^2 + ones(m,1)*Vary_k),2) - Ey.^2;
    #         solution.Vary <- Vary;


    # BIC AIC and ICL
    df <- length(psi)
    solution$stats$df <- df

    solution$stats$BIC <- solution$stats$ml - (df * log(n) / 2)
    solution$stats$AIC <- solution$stats$ml - df
    ## CL: complete-data loglikelihood
    zik_log_alpha_Phi_xy <- Zik * log_alpha_Phi_xy
    sum_zik_log_fik <- apply(X = zik_log_alpha_Phi_xy, MARGIN = 1, sum)
    comp_loglik <- sum(sum_zik_log_fik)
    solution$stats$CL <- comp_loglik
    solution$stats$ICL <- solution$stats$CL - (df * log(n) / 2)
    ##
    if (nb_EM_runs > 1 && verbose) {
      message("Max value of the log-likelihood: ",
              solution$stats$ml,
              "\n\n")
    }
    if (loglik > best_loglik) {
      best_solution <- solution
      best_loglik <- loglik
    }
  }

  solution <- best_solution

  if (nb_EM_runs > 1 && verbose) {
    message("Max value of the log-likelihood: ", solution$stats$ml, "\n")
  }

  solution$stats$cputime <- mean(stored_cputime)
  solution$stats$stored_cputime <- stored_cputime

  return(solution)
}
