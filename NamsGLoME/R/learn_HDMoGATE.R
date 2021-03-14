#' @export
learn_HDMoGATE = function(X, y, K, lambda, gamma, nb_EM_runs = 1, verbose = TRUE) {
  # Fits a mixture of high-dimensional gaussian-gated mixture-of-experts by an EM-Lasso algorithm. The
  # algorithm maximizes an l1-regularized log-likelihood where the regularized parameters are the
  # regression coefficients for the Gaussian regressor expert network, and the means, for the
  # Gaussian-gating network

  n <- nrow(X)
  p <- ncol(X)

  if (dim(y)[1] != n) {
    y <- t(y)
  }

  max_iter_EM <- 2000
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
        # log_Phi_x <- mvgaussian_pdf(X, mu[k,], R[, , k], covtype = "diagonal")$log_fxi
        log_Phi_x <- mvgaussian_pdf(X, mu[k,], R[, , k], covtype = "full")$log_fxi

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
        tauk <- tau[, k] # weights

        # Gating Network

        # The mixing proportions
        alpha[k] <- sum(tauk) / n

        # the Gaussian's mean and covariance matrix
        out <- solve_weighted_LASSOgating(X, gamma, tauk, alpha[k], mu[k,], nu2[k,], verbose = FALSE)
        mu[k,] <- out$mu # Gaussian mean update
        nu2[k,] = out$nu2 # Covariance diagonal matrix elements update
        R[, , k] <- out$R # Covariance matrix update

        # Experts Network
        out <- solve_weighted_LASSOexpert(X, y, lambda, tauk, beta0[k], beta[, k], sigma2[k], verbose = FALSE)
        beta[, k] <- out$beta # Regression coefficients vector update
        beta0[k] <- out$beta0 # Intercept update
        sigma2[k] <- out$sigma2 # Regressor variance update
      }

      # Lasso-penalized log-likelihood
      pen <- gamma * sum(abs(mu)) + lambda * sum(abs(beta))
      loglik <- sum(log_sum_alpha_PhixPhiy) - pen

      if (verbose) {
        message("EM - HDMoE: Iteration: ", iter, " | log-likelihood: ", loglik)
      }

      # convergence test
      converge <- abs((loglik - prev_loglik)) <= threshold
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
    # Tau <- Tau(1:m,:);
    # solution.Piik <- Piik;
    solution$stats <- list()
    solution$stats$tau <- tau
    solution$stats$log_Phi_x <- log_Phi_x
    solution$stats$log_Phi_y <- log_Phi_y
    solution$stats$log_alpha_Phi_xy <- log_alpha_Phi_xy
    solution$stats$ml <- loglik
    solution$stats$stored_loglik <- stored_loglik

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
    #         Ey_k = X(1:m,:)*beta;
    #         solution.Ey_k = Ey_k;
    #         # E[yi]
    #         Ey = sum(Piik.*Ey_k,2);
    #         solution.Ey = Ey;
    #
    #         # Var[yi|zi=k]
    #         Vary_k = sigma2;
    #         solution.Vary_k = Vary_k;
    #
    #         # Var[yi]
    #         Vary = sum(Piik.*(Ey_k.^2 + ones(m,1)*Vary_k),2) - Ey.^2;
    #         solution.Vary = Vary;


    # BIC AIC and ICL
    if ((lambda == 0) && (gamma == 0)) {
      df <- length(psi)
    } else if (gamma == 0) {
      df <- length(c(as.vector(alpha), as.vector(mu), as.vector(nu2), as.vector(beta0), as.vector(beta[beta != 0]), as.vector(sigma2)))
    } else if (lambda == 0) {
      df <- length(c(as.vector(alpha), as.vector(mu[mu != 0]), as.vector(nu2), as.vector(beta0), as.vector(beta), as.vector(sigma2)))
    } else {
      df <- length(c(as.vector(alpha), as.vector(mu[mu != 0]), as.vector(nu2), as.vector(beta0), as.vector(beta[beta != 0]), as.vector(sigma2)))
    }
    solution$stats$df <- df

    solution$stats$BIC <- solution$stats$ml + pen - (df * log(n) / 2)
    solution$stats$AIC <- solution$stats$ml + pen - df
    ## CL: complete-data loglikelihood
    zik_log_alpha_Phi_xy <- Zik * log_alpha_Phi_xy
    sum_zik_log_fik <- apply(X = zik_log_alpha_Phi_xy, MARGIN = 1, sum)
    comp_loglik <- sum(sum_zik_log_fik)
    solution$stats$CL <- comp_loglik
    solution$stats$ICL <- solution$stats$CL - (df * log(n) / 2)
    ##
    if (nb_EM_runs > 1 && verbose) {
      message("Max value of the log-likelihood: ", solution$stats$ml, "\n\n")
    }
    if (loglik > best_loglik) {
      best_solution <- solution
      best_loglik <- loglik
    }
  }

  solution <- best_solution

  if (nb_EM_runs > 1 && verbose) {
    message("Max value of the log-likelihood: ", solution$stats$loglik, "\n")
  }

  solution$stats$cputime <- mean(stored_cputime)
  solution$stats$stored_cputime <- stored_cputime

  return(solution)

}
