rm(list = ls())

library(HDMoGATE)

# Sample specification
n <- 300
p <- 8

alpha <- c(1/2, 1/2)

# Gating network parameters
mu <- matrix(data = NA, nrow = 2, ncol = p)
mu[1, ] <- c(0, 1, -1, -1.5, 0, 0.5, 0, 0)
mu[2, ] <- c(2, 0, 1, -1.5, 0, -0.5, 0, 0)

nu2 <- matrix(data = 1, nrow = 2, ncol = p)
r <- array(NA, dim = c(p, p, 2))

A <- matrix(runif(p ^ 2) * 2 - 1, ncol = p)
tmp <- t(A) %*% A

mask <- matrix(sample(x = 0:1, size = p * p, prob = c(0.6, 0.4), replace = TRUE), nrow = p) + diag(p)
mask[mask >= 2] <- 1

r[, , 1] <- diag(nrow = p)
r[, , 2] <- diag(nrow = p)

# Experts' network parameters
beta0 <- c(1, 1)
beta <- matrix(data = NA, nrow = p, ncol = 2)
beta[, 1] <- c(0, 1.5, 0, 0, 0, 1, 0, -0.5)
beta[, 2] <- c(1, -1.5, 0, 0, 2, 0, 0, 0.5)

sigma <- matrix(data = 1, nrow = 1, ncol = 2)

groundpsi <- c(as.vector(alpha), as.vector(mu), as.vector(nu2), as.vector(beta0), as.vector(beta), as.vector(sigma))

###############################################################################
# Generate samples
###############################################################################

echan <- HDMoGATE::sample_MoGATE(alpha, mu, r, beta0, beta, sigma, n)

###############################################################################
# Estimation
###############################################################################

estimation <- HDMoGATE::learn_HDMoGATE(echan$X, echan$y, K = 2, lambda = 1, gamma = 2, nb_EM_runs = 1, verbose = TRUE)

estimation$param$R
estimation$param$R[, , 2] - r[, , 1]


###############################################################################
# Some plots
###############################################################################

HDMoGATE::show_HDMoGATE_results(estimated_model = estimation, mu = mu, beta = beta)
