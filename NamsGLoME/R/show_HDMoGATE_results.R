#' @export
show_HDMoGATE_results = function(estimated_model, mu = NULL, beta = NULL) {

  K <- length(estimated_model$param$alpha)

  if (!is.null(mu) && !is.null(beta)) {
    graphics::matplot(t(mu), type = "p", col = "black", pch = "o", ylab = "Coefficient value")

    for (k in 1:K) {
      graphics::points(estimated_model$param$mu[k,], col = "red", pch = "+")
    }
    graphics::title(bquote(mu ~ " coeffcients: True values and EM-LASSO estimates"))
    graphics::legend("topright", legend = c("True", "EM-LASSO estimates"), col = c("black", "red"), pch = c("o", "+"), cex = 0.8)

    graphics::matplot(beta, type = "p", col = "black", pch = "o", ylab = "Coefficient value")

    for (k in 1:K) {
      graphics::points(estimated_model$param$beta[, k], col = "red", pch = "+")
    }
    graphics::title(bquote(beta ~ " coeffcients: True values and EM-LASSO estimates"))
    graphics::legend("topright", legend = c("True", "EM-LASSO estimates"), col = c("black", "red"), pch = c("o", "+"), cex = 0.8)

  } else {
    graphics::matplot(t(estimated_model$param$mu), type = "p", col = "red", pch = "+", ylab = "Coefficient value")
    graphics::title(bquote(mu ~ " coeffcients: True values and EM-LASSO estimates"))

    graphics::matplot(estimated_model$param$beta, type = "p", col = "red", pch = "+", ylab = "Coefficient value")
    graphics::title(bquote(beta ~ " coeffcients: True values and EM-LASSO estimates"))
  }

  plot(estimated_model$stats$stored_loglik, type = "l", xlab = "EM-Lasso epoch", ylab = "l_1-penalized log-likelihood")
}
