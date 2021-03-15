#' @export
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
