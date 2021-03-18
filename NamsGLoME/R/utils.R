
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

loggausspdf = function(X, mu, Sigma){
  d = nrow(X) ; n = ncol(X) ;
  if (ncol(mu) == n)
  {X=X-mu;}#sweep(X,2,mu,"-")}
  else {X=sweep(X,1,mu,"-")}
  #X = #scale(X,center=mu,scale=FALSE) # d x n  ### X ou t(X)
  p = matrixcalc::is.positive.definite(Sigma,.Machine$double.eps)
  if (! p) {
    print("SNPD !");
    y= rep(-Inf,n)
  } else {
    U = chol(Sigma) # d x d
    Q = solve(t(U),X) # d x n
    q = colSums(Q^2) # 1xn quadratic term (M distance)
    c = d*log(2*pi)+2*sum(log(diag(U))) #1x1 normalization constant
    y = -(c+q)/2;} # 1xn
  return(y)
}

logsumexp = function(x,dimension=c(1,2)){
  # % Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
  # %   By default dim = 1 (row).
  if (is.null(dimension)) dimension=1 ;
  # % subtract the largest in each row
  y = apply(x,dimension,max)
  x = x-y#sweep(x,dimension,y,"-")
  #s = y+log(rowSums(exp(x)))
  s = y+ log(apply(exp(x), dimension, sum))
  i = is.infinite(y)
  if (sum(i) > 0) s[i]=y[i]
  return(s)
}

MAP = function(postProb) {

  ###########################################################################
  # function [klas, Z] = MAP(postProb)
  #
  # calculate a partition by applying the Maximum A Posteriori Bayes
  # allocation rule
  #
  #
  # Inputs :
  #   postProb, a matrix of dimensions [n x K] of the posterior
  #  probabilities of a given sample of n observations arizing from K groups
  #
  # Outputs:
  #   klas: a vector of n class labels (z_1, ...z_n) where z_i =k \in {1,...K}
  #       klas(i) = arg   max (postProb(i,k)) , for all i=1,...,n
  #                     1<=k<=K
  #               = arg   max  p(zi=k|xi;theta)
  #                     1<=k<=K
  #               = arg   max  p(zi=k;theta)p(xi|zi=k;theta)/sum{l=1}^{K}p(zi=l;theta) p(xi|zi=l;theta)
  #                     1<=k<=K
  #
  #
  #       Z : Hard partition data matrix [nxK] with binary elements Zik such
  #       that z_ik =1 iff z_i = k
  #
  ######################### Faicel Chamroukhi ########################################

  n <- nrow(postProb)
  K <- ncol(postProb)

  # Maximum a posteriori rule
  klas <- apply(X = postProb, MARGIN = 1, FUN = which.max)
  partition_MAP <- (klas %*% ones(1,K)) == (ones(n,1) %*% 1:K)
  Z <- partition_MAP

  return(list(klas = klas, Z = Z))
}

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


mvGMM_pdf <- function(X, gmm) {
  #
  # calcul d'une densite melange
  # Entrees :
  #       X         : tableau de n individus d dimentionnels [nxd]
  #       parametre : structure dont les chaps sont:
  #                  1. pik les K proportions
  #                  2. mk les K centres . matrice de dim dxK
  #                  3: sigmak = (sigma1,...,sigmak) : les
  #                     variances des K composants du melange. matrice de dim vecteur de dim [dxdxK]
  #
  #
  #
  # #########################################################################

  n <- nrow(X)

  K <- length(gmm$weights)

  log_pik_fik = matrix(data = 0, nrow = n, ncol = K)
  for (k in 1:K) {
    pik <- gmm$weights[k]
    muk <- gmm$means[k, ]
    sigmak <- gmm$covariances[, , k]

    log_fik <- mvgaussian_pdf(X, muk, sigmak)$log_fxi
    log_pik_fik[, k] <- ones(n, 1) * log(pik) + log_fik # nx1
  }

  # pik_fik = max(pik_fik,eps);
  #log_pik_fik = max(log_pik_fik,log(realmin));

  #pik_fik = exp(log_pik_fik);
  #fxi = sum(pik_fik,2);
  #log_fxi = log(fxi);

  return(list(log_pik_fik = log_pik_fik, log_fik = log_fik))

}

posterior_mvGMM = function(X, gmm) {
  #
  #
  #
  # calcul des probas a posteriori dans pour un melange de gaussiennes ? K
  # composants
  #
  # Entrees :
  #       X         : tableau de n individus d dimentionnels [nxd]
  #       para : structure dont les chaps sont:
  #                  1. pik les K proportions
  #                  2. mk les K centres . matrice de dim dxK
  #                  3: sigmak = (sigma1,...,sigmak) : les
  #                     variances des K composants du melange. matrice de dim vecteur de dim [dxdxK]
  #
  #
  # ###################################

  K <- length(gmm$weights)

  # [log_pik_fik, pik_fik, log_fxi, fxi] = mvGMM_pdf(X, gmm);
  log_pik_fik <- mvGMM_pdf(X, gmm)$log_pik_fik # nxK

  #log_pik_fik = min(log_pik_fik,log(realmax));
  #log_pik_fik = max(log_pik_fik,log(realmin));


  # log_sum_pik_fik = log(sum(exp(log_pik_fik),2));
  log_sum_pik_fik <- logsumexp(log_pik_fik, 1)  # nx1
  log_tauik <- log_pik_fik - log_sum_pik_fik %*% ones(1, K) # nxK
  tau = exp(log_tauik)  # nxK
  #tau <- exp(log_tauik) / (apply(exp(log_tauik), 1, sum) %*% ones(1, K)) # nxK
  loglik <- sum(log_sum_pik_fik) # 1x1

  return(list(tau = tau, loglik = loglik, log_pik_fik = log_pik_fik))
}

designmatrix = function(x, p, q = NULL, n = 1) {

  order_max <- p
  if (!is.null(q)) {
    order_max <- max(p, q)
  }

  X <- matrix(NA, length(x), order_max + 1)
  for (i in 0:(order_max)) {
    X[, i + 1] <- x ^ i
  }

  XBeta <- X[, 1:(p + 1)]
  # design matrix for Beta (the polynomial regressions)
  if (!is.null(q)) {
    Xw <- X[, 1:(q + 1)]
    Xw <- repmat(Xw, n, 1)
    # design matrix for w (the logistic regression)
  } else {
    Xw <- NULL
  }

  XBeta <- repmat(XBeta, n, 1)

  return(list(Xw = Xw, XBeta = XBeta))
}

ones <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(1, n, d))
  }
  else{
    return(array(1, dim = c(n, d, g)))
  }
}

zeros <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(0, n, d))
  }
  else{
    return(array(0, dim = c(n, d, g)))
  }
}

rand <- function(n, d, g = 1) {
  if (g == 1) {
    return(matrix(stats::runif(n * d), n, d))
  }
  else{
    return(array(stats::runif(n * d), dim = c(n, d, g)))
  }
}

repmat <- function(M, n, d) {
  return(kronecker(matrix(1, n, d), M))
}

drnorm <- function(n, d, mean, sd) {
  A <- matrix(nrow = n, ncol = d)
  for (i in 1:d) {
    A[, i] <- stats::rnorm(n, mean, sd)
  }
  return(A)
}

lognormalize <- function(M) {
  if (!is.matrix(M)) {
    M <- matrix(M)
  }
  n <- nrow(M)
  d <- ncol(M)
  a <- apply(M, 1, max)
  return(M - repmat(a + log(rowSums(exp(M - repmat(a, 1, d)))), 1, d))
}

normalize <- function(A, dim) {
  # Normalize makes the entries of a (multidimensional <= 2) array sum to 1.
  # Input
  # A = Array to be normalized
  # dim = dimension is specified to normalize.
  # Output
  # M = Array after normalize.
  # z is the normalize constant
  # Note:
  # If dim is specified, we normalize the specified dimension only,
  # Otherwise we normalize the whole array.
  # Dim = 1 normalize each column
  # Dim = 2 normalize each row

  if (nargs() < 2) {
    z <- sum(A)
    # Set any zeros to one before dividing
    # This is valid, since c = 0 ==> all i.A[i] = 0 ==> the anser should be 0/1 = 0.
    s <- z + (z == 0)
    M <- A / s
  } else if (dim == 1) {
    # normalize each column
    z <- colSums(A)
    s <- z + (z == 0)
    M <- A / matrix(s, nrow = dim(A)[1], ncol = length(s), byrow = TRUE)
  } else{
    z <- rowSums(A)
    s <- z + (z == 0)
    M <- A / matrix(s, ncol = dim(A)[2], nrow = length(s), byrow = FALSE)
  }
  output <- list(M = M, z = z)
  return(output)
}
