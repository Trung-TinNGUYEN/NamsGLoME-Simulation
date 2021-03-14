#' @export
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
    log_pik_fik[, k] <- ones(n, 1) * log(pik) + log_fik
  }

  # pik_fik = max(pik_fik,eps);
  #log_pik_fik = max(log_pik_fik,log(realmin));

  #pik_fik = exp(log_pik_fik);
  #fxi = sum(pik_fik,2);
  #log_fxi = log(fxi);

  return(list(log_pik_fik = log_pik_fik, log_fik = log_fik))

}
