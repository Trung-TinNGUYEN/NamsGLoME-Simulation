#' @export
# Permute a HDMO model if needed (K = 2 exclusively)
permute = function(model, sample, criterion = c("l2", "error.classif")) {
  criterion = match.arg(criterion)

  K <- length(model$param$alpha)
  n <- length(sample$y)

  classifRateModel <- sum(model$stats$klas == sample$stats$klasy) / n
  classifRateRevModel <- sum(abs(model$stats$klas - 2) + 1 == sample$stats$klasy) / n

  if (classifRateModel < classifRateRevModel) {
    index <- c(2, 1)

    model$param$alpha <- model$param$alpha[index]
    model$param$mu <- model$param$mu[index,]

    model$param$nu2 <- model$param$nu2[index,]

    for (k in 1:K) {
      model$param$R[, , k] <- diag(model$param$nu2[k,])
    }

    model$param$beta0 <- model$param$beta0[index]

    model$param$beta <- model$param$beta[, index]

    model$param$sigma2 <- model$param$sigma2[index]

    model$stats$psi <- c(as.vector(model$param$alpha), as.vector(model$param$mu), as.vector(model$param$nu2), as.vector(model$param$beta0), as.vector(model$param$beta), as.vector(model$param$sigma2))
    model$stats$klas <- abs(model$stats$klas - 2) + 1

  }

  return(model)

}
