
# Removes all objects from the current workspace (R memory).
#rm(list = ls())

###########################################################################
#                             Load packages
###########################################################################

# Install the NamsGLoME package for local machine.
# install.packages("devtools")
devtools::install("NamsGLoME")
library(NamsGLoME)

# Load CRAN Packages used in this paper.
#install.packages("capushe")
library(capushe)
#install.packages("xLLiM")
library(xLLiM)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("matrixcalc")
library(matrixcalc)

###########################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%% Simulated data set %%%%%%%%%%%%%%%%%%%%%%%%%%
###########################################################################

#######################################
#         Setting parameters
#######################################
# Different number of trials for several numerical experiments.

# Comparison histograms of selected K between WS and MS cases using jump and
# slope criteria over 100 trials.
#numTrial <- 100

# Box-plot of the tensorized Kullback-Leibler divergence according to the number of
# mixture components using the jump criterion over 100 trials.
#numTrial <- 100

# Tensorized Kullback-Leibler divergence between the true and selected densities based
# on the jump criterion, represented in a log-log scale, using 30 trials.
#numTrial <- 30

numTrial <- 2

# Sample size used for histograms of selected K
#n <- 2000
#n <- 10000

# Generate a vector of 11 elements in increasing order.
# Each element represents for different sample size used for studying the error decays.
n_size <- round(1000*2.^((0:10)/3))

# Resulted n_size: 1000  1260  1587  2000  2520  3175  4000  5040  6350  8000 10079

# Number data points for approximating KL using Monte Carlo method.
ny <- 30

# Collection of model based on the numer of components K of GLLiM.
# K = 1,...,Kmax
Kmax <- 20

#############################################################################################
# Saving data from experiments
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Create matrice containing the selected model (K-number of components of MoE)
selectedK_Djump_WS_JNLL <- matrix(0, nrow = numTrial, ncol = length(n_size))
selectedK_DDSE_WS_JNLL <- matrix(0, nrow = numTrial, ncol = length(n_size))

selectedK_Djump_MS_JNLL <- matrix(0, nrow = numTrial, ncol = length(n_size))
selectedK_DDSE_MS_JNLL <- matrix(0, nrow = numTrial, ncol = length(n_size))

dataCapushe_MoEList_WS_JNLL <- list()
dataCapushe_MoEList_MS_JNLL <- list()

# Using conditional log-likelihood

selectedK_Djump_WS_CNLL <- matrix(0, nrow = numTrial, ncol =  length(n_size))
selectedK_DDSE_WS_CNLL <- matrix(0, nrow = numTrial, ncol =  length(n_size))

selectedK_Djump_MS_CNLL <- matrix(0, nrow = numTrial, ncol =  length(n_size))
selectedK_DDSE_MS_CNLL <- matrix(0, nrow = numTrial, ncol =  length(n_size))

dataCapushe_MoEList_WS_CNLL <- list()
dataCapushe_MoEList_MS_CNLL <- list()


#--------------------------------------------------------------------------------------------
# Choosing the true parameter for simulated data sets
#------------------------------------------------------
p <- 1 # dimension of data

alpha <- c(1/2, 1/2)

# When X is sampled from Gaussian distribution
# Gating network parameters
mu <- matrix(data = NA, nrow = 2, ncol = p)
mu[1, ] <- c(0.2)
mu[2, ] <- c(0.8)

nu2 <- matrix(data = 1, nrow = 2, ncol = p)

r <- array(NA, dim = c(p, p, 2))
r[, , 1] <- 0.1*diag(nrow = p)
r[, , 2] <- 0.15*diag(nrow = p)

# Experts' network parameters:

# Well-Specified case: WS (GLLiM)

beta0_WS <- c(2, 0)
beta_WS <- matrix(data = NA, nrow = p, ncol = 2)
# beta[, 1] <- c(0, 1.5, 0, 0, 0, 1, 0, -0.5)
# beta[, 2] <- c(1, -1.5, 0, 0, 2, 0, 0, 0.5)

beta_WS[, 1] <- c(-5)
beta_WS[, 2] <- c(0.1)

sigma <- matrix(data = 0.3, nrow = 1, ncol = 2)
# Vector containing the standard deviation of Gaussian distribution.


# # When X is sampled from truncated Gaussian distribution. X shape.
# # Gating network parameters
# mu <- matrix(data = NA, nrow = 2, ncol = p)
# mu[1, ] <- c(0.3)
# mu[2, ] <- c(0.6)
#
# nu2 <- matrix(data = 1, nrow = 2, ncol = p)
# r <- array(NA, dim = c(p, p, 2))
#
# A <- matrix(runif(p ^ 2) * 2 - 1, ncol = p)
# tmp <- t(A) %*% A
#
# mask <- matrix(sample(x = 0:1, size = p * p, prob = c(0.6, 0.4), replace = TRUE), nrow = p) + diag(p)
# mask[mask >= 2] <- 1
#
# r[, , 1] <- 0.6*diag(nrow = p)
# r[, , 2] <- 0.3*diag(nrow = p)
#
# # Experts' network parameters
#
#
# # Well-Specified case: WS (GLLiM)
#
# beta0_WS <- c(2, -2)
# beta_WS <- matrix(data = NA, nrow = p, ncol = 2)
# beta_WS[, 1] <- c(-4)
# beta_WS[, 2] <- c(4)
#
# sigma <- matrix(data = 0.3, nrow = 1, ncol = 2)


# Misspecified case: MS (MoE with parabolic mean of Experts)

# Experts' network parameters
beta0_MS <- c(1, 0)
beta_MS <- matrix(data = NA, nrow = p, ncol = 2)

beta_MS[, 1] <- c(-6)
beta_MS[, 2] <- c(0)

beta2_MS <- matrix(data = NA, nrow = p, ncol = 2)

beta2_MS[, 1] <- c(1)
beta2_MS[, 2] <- c(-0.4)


# # # When X is sampled from truncated Gaussian distribution. X shape.
# # Experts' network parameters
# beta0_MS <- c(2, -2)
#
# beta_MS <- matrix(data = NA, nrow = p, ncol = 2)
# beta_MS[, 1] <- c(-1.5)
# beta_MS[, 2] <- c(1.5)
#
# beta2_MS <- matrix(data = NA, nrow = p, ncol = 2)
# beta2_MS[, 1] <- c(-3)
# beta2_MS[, 2] <- c(3)

############################################################################################
#--------------------------------------------------------------------------------------------
# Generating the simulated data set and estimating the parameters
#--------------------------------------------------------------------------------------------

# Saving the simulation data set in list() format
simData_WS <- list()
simData_MS <- list()

# Sample the fixed data set with size 10079.
echan_WS_test <- HDMoGATE::sample_GLLiM(alpha, mu, r, beta0_WS, beta_WS, sigma, n_size[length(n_size)], 1)
echan_MS_test <- HDMoGATE::sample_MoGATEParabola(alpha, mu, r, beta0_MS, beta_MS, beta2_MS, sigma, n_size[length(n_size)], 1)

# Collection of model based on the numer of components K of GLLiM.
# K = 1,...,Kmax

modelComplex_WS <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_WS_JNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_WS_CNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)

modelComplex_MS <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_MS_JNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_MS_CNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)


# Initalize matrices tensorized (Jensen)-Kullback-Leibler divergence
# between true and estimated density

# WS case
# KL over all trials, models and n_size samples
KL_WS <- Inf*array(1, dim = c(length(n_size), Kmax, numTrial))

# KL for the selected model over all trials on n_size samples  using Djump
KL_Sel_WS_Djump_JNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using DDSE
KL_Sel_WS_DDSE_JNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using Djump
KL_Sel_WS_Djump_CNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using Djump
KL_Sel_WS_DDSE_CNLL <- matrix(0, numTrial, length(n_size))

# For JKL
# JKL over all trials, models and n_size samples
JKL_WS <- Inf*array(1, dim = c(length(n_size), Kmax, numTrial))

# JKL for the selected model over all trials on n_size samples  using Djump
JKL_Sel_WS_Djump_JNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using DDSE
JKL_Sel_WS_DDSE_JNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using Djump
JKL_Sel_WS_Djump_CNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using Djump
JKL_Sel_WS_DDSE_CNLL <- matrix(0, numTrial, length(n_size))


# MS case
# KL over all trials, models and n_size samples
KL_MS <- Inf*array(1, dim = c(length(n_size), Kmax, numTrial))

# KL for the selected model over all trials on n_size samples  using Djump
KL_Sel_MS_Djump_JNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using DDSE
KL_Sel_MS_DDSE_JNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using Djump
KL_Sel_MS_Djump_CNLL <- matrix(0, numTrial, length(n_size))

# KL for the selected model over all trials on n_size samples using Djump
KL_Sel_MS_DDSE_CNLL <- matrix(0, numTrial, length(n_size))

# For JKL
# JKL over all trials, models and n_size samples
JKL_MS <- Inf*array(1, dim = c(length(n_size), Kmax, numTrial))

# JKL for the selected model over all trials on n_size samples  using Djump
JKL_Sel_MS_Djump_JNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using DDSE
JKL_Sel_MS_DDSE_JNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using Djump
JKL_Sel_MS_Djump_CNLL <- matrix(0, numTrial, length(n_size))

# JKL for the selected model over all trials on n_size samples using Djump
JKL_Sel_MS_DDSE_CNLL <- matrix(0, numTrial, length(n_size))


#########################################################################################
# Generate samples and plot samples
########################################################################################
# # Generate samples from GLLiM
#
# echan_WS <- HDMoGATE::sample_GLLiM(alpha, mu, r, beta0_WS, beta_WS, sigma, n, ny)
# simData_WS[[t]] <- echan_WS
#
#
# # Generate samples from parabolic mean of Experts
#
#
# echan_MS <- HDMoGATE::sample_MoGATEParabola(alpha, mu, r, beta0_MS, beta_MS, beta2_MS, sigma, n, ny)
# simData_MS[[t]] <- echan_MS

#
# # Show MS typical realization
#
#
# # WS case:
# X_WS <- echan_WS$X
# Y_WS <- echan_WS$y
# E_WS <- as.factor(echan_WS$stats$klasy)
#
# df_WS <- data.frame(X_WS, Y_WS, E_WS)
# names(df_WS) <- c('X', 'Y', 'Class')
# ggplot(df_WS, aes(x=X , y=Y , color = Class)) + geom_point() + scale_color_manual(values=c("blue", "red")) + ggtitle("Typical realization in WS case") + theme(plot.title = element_text(hjust = 0.5))
#
# # # 3D-Plot
# # library(plotly)
# # plot_ly(x=X[, 1], y=X[, 2], z=Y, type="scatter3d", mode="markers", color=E)
# # plot_ly(x=X[, 1], y=Y, type="scatter", mode="markers", color=E)
#
# # MS case:
# X_MS <- echan_MS$X
# Y_MS <- echan_MS$y
# E_MS <- as.factor(echan_MS$stats$klasy)
#
# df_MS <- data.frame(X_MS, Y_MS, E_MS)
# names(df_MS) <- c('X', 'Y', 'Class')
# ggplot(df_MS, aes(x = X, y=Y , color = Class)) + geom_point() + scale_color_manual(values=c("blue", "red")) + ggtitle("Typical realization in MS case") + theme(plot.title = element_text(hjust = 0.5))
#


############################################################################################
# Running the simulation over numTrial trials over 11 sample sizes n=1000,2000,...,10079.
#------------------------------------------------------------------------------------------
# Save time for calculation
start_time <- Sys.time()

for (t in 1:numTrial) {
  #t <- 1 # testing code
  for (n in 1:length(n_size)){

    # Test code
    #n <- 1
    # WS case
    # Initalizing the whole data set for estimators and Monte Carlo method
    echan_WS <- HDMoGATE::sample_GLLiM(alpha, mu, r, beta0_WS, beta_WS, sigma, n_size[n] , ny)
    # simData_WS[[t]] <- echan_WS

    # MS case
    echan_MS <- HDMoGATE::sample_MoGATEParabola(alpha, mu, r, beta0_MS, beta_MS, beta2_MS, sigma, n_size[n], ny)
    # simData_MS[[t]] <- echan_MS

    # For each number of components (k = 1,...,Kmax), we use gllim function to estimate the parameters.
    for (k in 1:Kmax) {
      # Using function gllim from xLLiM package when k > 1
      if (k>1) {

        #k <- 2 # Tesing

        # WS case:
        #while (KL_WS[t,k] == Inf){

        # # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_WS <- HDMoGATE::sample_GLLiM(alpha, mu, r, beta0_WS, beta_WS, sigma, n, ny)
        # simData_WS[[t]] <- echan_WS

        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        estimaMoEGLLiM_WS <- gllim(t(as.matrix(echan_WS_test$y[1:n_size[n]])),
                                   t(as.matrix(echan_WS_test$X[1:n_size[n]])), in_K = k)

        # Using n data points for estimators.
        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        gllim_inverse_WS <- gllim_inverse_dens(t(as.matrix(echan_WS_test$X[1:n_size[n]])),
                                               estimaMoEGLLiM_WS,t(as.matrix(echan_WS_test$y[1:n_size[n]])))

        modelComplex_WS[n, k] <- estimaMoEGLLiM_WS$nbpar
        contrast_WS_JNLL[n, k] <- -estimaMoEGLLiM_WS$LLf
        contrast_WS_CNLL[n, k] <- gllim_inverse_WS$CNLL

        # Calculate (Jensen)-Kullback-Leibler divergence JKL(dens1,dens2)
        # (or KL(dens1,dens2)) by Monte Carlo method using ny*n points

        # JKL case: Given rho such that 0 < rho < 1,
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho
        # KL(dens1(./x_i), (1-rho)*dens1(./x_i) + rho*dens2(./x_i))
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho 1/ny sum_{j=1..ny}
        # log(dens1(y_ij/x_i)/(1-rho)*dens1(./x_i) + rho*dens2(./x_i))
        # with y_ij is sampled from dens1(./x_i)

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        CLL_WS_vec <- gllim_inverse_dens(t(echan_WS$Xny),estimaMoEGLLiM_WS,t(echan_WS$yny))$CLL_vec

        JKL_WS[n,k,t] <- 1/(n_size[n]*ny*rho)*sum(log(echan_WS$gllim_pdf_vec/((1-rho)*echan_WS$gllim_pdf_vec+rho*CLL_WS_vec)))

        # KL case:
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} KL(dens1(./x_i),dens2(./x_i))
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/ny sum_{j=1..ny}
        # log(dens1(y_ij/x_i)/dens2(y_ij/x_i))
        # with y_ij is sampled from dens1(./x_i)

        # Using n*ny data points for Monte  Carlo.
        KL_WS[n,k,t] <- echan_WS$gllim_pdf+gllim_inverse_dens(t(echan_WS$Xny),estimaMoEGLLiM_WS,t(echan_WS$yny))$CNLL
        KL_WS[n,k,t] <- KL_WS[n,k,t]/(n_size[n]*ny)
        #}

        # MS case:

        #while (KL_MS[t,k] == Inf){
        # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_MS <- HDMoGATE::sample_MoGATEParabola(alpha, mu, r, beta0_MS, beta_MS, beta2_MS, sigma, n, ny)
        # simData_MS[[t]] <- echan_MS

        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        estimaMoEGLLiM_MS <- gllim(t(as.matrix(echan_MS_test$y[1:n_size[n]])),
                                   t(as.matrix(echan_MS_test$X[1:n_size[n]])), in_K = k)

        # Using n data points for estimators.
        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        gllim_inverse_MS <- gllim_inverse_dens(t(as.matrix(echan_MS_test$X[1:n_size[n]])),
                                               estimaMoEGLLiM_MS,t(as.matrix(echan_MS_test$y[1:n_size[n]])))

        modelComplex_MS[n, k] <- estimaMoEGLLiM_MS$nbpar
        contrast_MS_JNLL[n, k] <- -estimaMoEGLLiM_MS$LLf
        contrast_MS_CNLL[n, k] <- gllim_inverse_MS$CNLL

        # Calculate (Jensen)-Kullback-Leibler divergence JKL(dens1,dens2)
        # (or KL(dens1,dens2)) by Monte Carlo method using ny*n points

        # JKL case: Given rho such that 0 < rho < 1,
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho
        # KL(dens1(./x_i), (1-rho)*dens1(./x_i) + rho*dens2(./x_i))
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho 1/ny sum_{j=1..ny}
        # log(dens1(y_ij/x_i)/(1-rho)*dens1(./x_i) + rho*dens2(./x_i))
        # with y_ij is sampled from dens1(./x_i)

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        CLL_MS_vec <- gllim_inverse_dens(t(echan_MS$Xny),estimaMoEGLLiM_MS,t(echan_MS$yny))$CLL_vec

        JKL_MS[n,k,t] <- 1/(n_size[n]*ny*rho)*sum(log(echan_MS$moe2_pdf_vec/((1-rho)*echan_MS$moe2_pdf_vec+rho*CLL_MS_vec)))

        # KL case:
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} KL(dens1(./x_i),dens2(./x_i))
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/ny sum_{j=1..ny}
        # log(dens1(y_ij/x_i)/dens2(y_ij/x_i))
        # with y_ij is sampled from dens1(./x_i)

        # Using n*ny data points for Monte  Carlo.
        KL_MS[n,k,t] <- echan_MS$moe2_pdf+gllim_inverse_dens(t(echan_MS$Xny),estimaMoEGLLiM_MS,t(echan_MS$yny))$CNLL
        KL_MS[n,k,t] <- KL_MS[n,k,t]/(n_size[n]*ny)
        #}

        # k = 1 using lm function from R for Fitting Generalized Linear Models
        # Using MLE for linear regression in R
        # http://daviddalpiaz.github.io/appliedstats/simple-linear-regression.html#variance-estimation
      } else {
        # WS case:

        #while (KL_WS[t,k] == Inf){

        # # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_WS <- HDMoGATE::sample_GLLiM(alpha, mu, r, beta0_WS, beta_WS, sigma, n, ny)
        # simData_WS[[t]] <- echan_WS

        # Using lm function from R for Fitting Generalized Linear Models
        estimaMoE_WS_lm <- lm(as.matrix(echan_WS_test$y[1:n_size[n]]) ~ as.matrix(echan_WS_test$X[1:n_size[n]]))
        modelComplex_WS[n, k] <- attributes(logLik(estimaMoE_WS_lm))$df
        contrast_WS_JNLL[n, k] <- -logLik(estimaMoE_WS_lm)[1]
        contrast_WS_CNLL[n, k] <- -logLik(estimaMoE_WS_lm)[1]

        # # Testing function provides vectors of pdf over all samples from lm()
        # # mypdf and loglik() from R
        # num_obs = 30; sigma = 0.3
        # #set.seed(1)
        # epsilon = rnorm(n = num_obs, mean = 0, sd = sigma)
        # x_vals = seq(from = 0, to = 10, length.out = num_obs)
        #
        # sim_slr = function(x, beta_0 = 10, beta_1 = 5, sigma = 1) {
        #   n = length(x)
        #   epsilon = rnorm(n, mean = 0, sd = sigma)
        #   y = beta_0 + beta_1 * x + epsilon
        #   data.frame(predictor = x, response = y)
        # }
        #
        # sim_data = sim_slr(x = x_vals, beta_0 = 5, beta_1 = -2, sigma = 3)
        # x = sim_data$predictor; y = sim_data$response
        #
        # sim_fit = lm(response ~ predictor, data = sim_data)
        #
        # beta_0_hat = coef(sim_fit)[1]; beta_1_hat = coef(sim_fit)[2]
        #
        #
        # y_hat = beta_0_hat + beta_1_hat * x
        # e = y - y_hat
        # n = length(e)
        # s2_e = sum(e^2)/n
        # print(sum(sim_fit$residuals)/n)
        # print(sqrt(s2_e))
        #
        # mypdf <- sum(log(exp(loggausspdf(matrix(y,ncol = n), matrix(coef(sim_fit)[1] + coef(sim_fit)[2]*x,ncol = n),
        #                         matrix(s2_e,nrow = 1)))))
        # sprintf("mypdf = %s",mypdf)
        # logLik(sim_fit)


        # sum(log(exp(loggausspdf(t(echan_WS$y),
        #                         t(matrix(beta_0_hat + echan_WS$X%*% beta_1_hat,ncol=1)),
        #                         matrix(sum(y_hat^2)/n,ncol = 1)))))

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        y_hat_WS <- coef(estimaMoE_WS_lm)[1] + coef(estimaMoE_WS_lm)[2]*echan_WS$Xny
        s2_e_WS = sum((echan_WS$yny - y_hat_WS)^2)/(n_size[n]*ny)

        CLL_WS_vec <- exp(loggausspdf(t(echan_WS$yny), t(y_hat_WS), matrix(s2_e_WS,nrow = p)))
        estimaMoE_WS_lm_pdf <- sum(log(CLL_WS_vec))

        JKL_WS[n,k,t] <- 1/(n_size[n]*ny*rho)*sum(log(echan_WS$gllim_pdf_vec/((1-rho)*echan_WS$gllim_pdf_vec+rho*CLL_WS_vec)))

        KL_WS[n,k,t] <- 1/(n_size[n]*ny)*(echan_WS$gllim_pdf-estimaMoE_WS_lm_pdf)
        #}

        # MS case:
        #while (KL_MS[t,k] == Inf){

        # # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_MS <- HDMoGATE::sample_MoGATEParabola(alpha, mu, r, beta0_MS, beta_MS, beta2_MS, sigma, n, ny)
        # simData_MS[[t]] <- echan_MS

        # Using lm function from R for Fitting Generalized Linear Models
        estimaMoE_MS_lm <- lm(as.matrix(echan_MS_test$y[1:n_size[n]]) ~ as.matrix(echan_MS_test$X[1:n_size[n]]))

        modelComplex_MS[n, k] <- attributes(logLik(estimaMoE_MS_lm))$df
        contrast_MS_JNLL[n, k] <- -logLik(estimaMoE_MS_lm)[1]
        contrast_MS_CNLL[n, k] <- -logLik(estimaMoE_MS_lm)[1]

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        y_hat_MS <- coef(estimaMoE_MS_lm)[1] + coef(estimaMoE_MS_lm)[2]*echan_MS$Xny
        s2_e_MS = sum((echan_MS$yny - y_hat_MS)^2)/(n_size[n]*ny)

        CLL_MS_vec <- exp(loggausspdf(t(echan_MS$yny), t(y_hat_MS), matrix(s2_e_MS,nrow = p)))
        estimaMoE_MS_lm_pdf <- sum(log(CLL_MS_vec))

        JKL_MS[n, k, t] <- 1/(n_size[n]*ny*rho)*sum(log(echan_MS$moe2_pdf_vec/((1-rho)*echan_MS$moe2_pdf_vec+rho*CLL_MS_vec)))

        KL_MS[n, k, t] <- 1/(n_size[n]*ny)*(echan_MS$moe2_pdf-estimaMoE_MS_lm_pdf)

        #}
      }
    }
    #################################################
    # Using capushe package to select model model
    #################################################

    # Using joint log-likelihood
    # WS case:
    dataCapushe_MoE_WS_JNLL <- data.frame(c(1:Kmax), modelComplex_WS[n,], modelComplex_WS[n,], contrast_WS_JNLL[n,])
    names(dataCapushe_MoE_WS_JNLL) <- c("model", "pen", "complexity", "contrast")

    dataCapushe_MoEList_WS_JNLL[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_WS_JNLL)
    slopeHeur_MoE_WS <- capushe(dataCapushe_MoE_WS_JNLL)

    selectedK_Djump_WS_JNLL[t, n] <- slopeHeur_MoE_WS@Djump@model
    selectedK_DDSE_WS_JNLL[t, n] <- slopeHeur_MoE_WS@DDSE@model

    KL_Sel_WS_Djump_JNLL[t, n] <- KL_WS[n, as.numeric(selectedK_Djump_WS_JNLL[t, n]),t]
    KL_Sel_WS_DDSE_JNLL[t, n] <- KL_WS[n, as.numeric(selectedK_DDSE_WS_JNLL[t, n]),t]

    JKL_Sel_WS_Djump_JNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_Djump_WS_JNLL[t, n]),t]
    JKL_Sel_WS_DDSE_JNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_DDSE_WS_JNLL[t, n]),t]

    # MS case:
    dataCapushe_MoE_MS_JNLL <- data.frame(c(1:Kmax), modelComplex_MS[n,], modelComplex_MS[n,], contrast_MS_JNLL[n,])
    names(dataCapushe_MoE_MS_JNLL) <- c("model", "pen", "complexity", "contrast")

    dataCapushe_MoEList_MS_JNLL[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_MS_JNLL)
    slopeHeur_MoE_MS <- capushe(dataCapushe_MoE_MS_JNLL)

    selectedK_Djump_MS_JNLL[t, n] <- slopeHeur_MoE_MS@Djump@model
    selectedK_DDSE_MS_JNLL[t, n] <- slopeHeur_MoE_MS@DDSE@model

    KL_Sel_MS_Djump_JNLL[t, n] <- KL_MS[n, as.numeric(selectedK_Djump_MS_JNLL[t, n]),t]
    KL_Sel_MS_DDSE_JNLL[t, n] <- KL_MS[n, as.numeric(selectedK_DDSE_MS_JNLL[t, n]),t]

    JKL_Sel_MS_Djump_JNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_Djump_MS_JNLL[t, n]),t]
    JKL_Sel_MS_DDSE_JNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_DDSE_WS_JNLL[t, n]),t]

    ####
    #plot(Djump(dataCapushe_MoE_WS_JNLL),newwindow=FALSE)
    #plot(DDSE(dataCapushe_MoE_WS_JNLL),newwindow=FALSE)

    ####
    # Using conditional log-likelihood

    # WS case:
    dataCapushe_MoE_WS_CNLL <- data.frame(c(1:Kmax), modelComplex_WS[n,], modelComplex_WS[n,], contrast_WS_CNLL[n,])
    names(dataCapushe_MoE_WS_CNLL) <- c("model", "pen", "complexity", "contrast")

    dataCapushe_MoEList_WS_CNLL[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_WS_CNLL)
    slopeHeur_MoE_WS <- capushe(dataCapushe_MoE_WS_CNLL)

    selectedK_Djump_WS_CNLL[t, n] <- slopeHeur_MoE_WS@Djump@model
    selectedK_DDSE_WS_CNLL[t, n] <- slopeHeur_MoE_WS@DDSE@model

    KL_Sel_WS_Djump_CNLL[t, n] <- KL_WS[n, as.numeric(selectedK_Djump_WS_CNLL[t, n]),t]
    KL_Sel_WS_DDSE_CNLL[t, n] <- KL_WS[n, as.numeric(selectedK_DDSE_WS_CNLL[t, n]),t]

    JKL_Sel_WS_Djump_CNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_Djump_WS_CNLL[t, n]),t]
    JKL_Sel_WS_DDSE_CNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_DDSE_WS_CNLL[t, n]),t]

    # MS case:
    dataCapushe_MoE_MS_CNLL <- data.frame(c(1:Kmax), modelComplex_MS[n,], modelComplex_MS[n,], contrast_MS_CNLL[n,])
    names(dataCapushe_MoE_MS_CNLL) <- c("model", "pen", "complexity", "contrast")

    dataCapushe_MoEList_MS_CNLL[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_MS_CNLL)
    slopeHeur_MoE_MS <- capushe(dataCapushe_MoE_MS_CNLL)

    selectedK_Djump_MS_CNLL[t, n] <- slopeHeur_MoE_MS@Djump@model
    selectedK_DDSE_MS_CNLL[t, n] <- slopeHeur_MoE_MS@DDSE@model

    KL_Sel_MS_Djump_CNLL[t, n] <- KL_MS[n, as.numeric(selectedK_Djump_MS_CNLL[t, n]),t]
    KL_Sel_MS_DDSE_CNLL[t, n] <- KL_MS[n, as.numeric(selectedK_DDSE_MS_CNLL[t, n]),t]

    JKL_Sel_MS_Djump_CNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_Djump_MS_CNLL[t, n]),t]
    JKL_Sel_MS_DDSE_CNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_DDSE_WS_CNLL[t, n]),t]

  }
}

end_time <- Sys.time()
runningTime <- end_time - start_time
print(runningTime)
