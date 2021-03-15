##########################################################################################
#       Comparison histograms of selected GLoME models between well-specified (WS)
#       and misspecified (MS) cases using jump and slope criteria over 100 trials.
##########################################################################################

# Removes all objects from the current workspace (R memory).
#rm(list = ls())

##########################################################################################
#                                         Load packages
##########################################################################################
  
# Install the NamsGLoME package for local machine.
# install.packages("devtools")
# devtools::document("NamsGLoME")
devtools::install("NamsGLoME")
library(NamsGLoME)

# Load CRAN Packages used in this paper.

# CAlibrating Penalities Using Slope HEuristics (CAPUSHE):
# The capushe function proposes two algorithms based on the slope heuristics
# to calibrate penalties in the context of model selection via penalization.
#install.packages("capushe")
library(capushe)

# High Dimensional Locally-Linear Mapping:
# Provides a tool for non linear mapping (non linear regression)
# using a mixture of regression model and an inverse regression strategy.
#install.packages("xLLiM")
library(xLLiM)

# Create Elegant Data Visualisations Using the Grammar of Graphics:
# A system for 'declaratively' creating graphics, based on "The Grammar of Graphics".
# You provide the data, tell 'ggplot2' how to map variables to aesthetics,
# what graphical primitives to use, and it takes care of the details.
#install.packages("ggplot2")
library(ggplot2)

# A collection of functions to support matrix calculations
# for probability, econometric and numerical analysis.
#install.packages("matrixcalc")
library(matrixcalc)

##########################################################################################
#                           Setting parameters for experiments
##########################################################################################

# Number of trials for numerical experiments.
#numTrial <- 100

numTrial <- 2

# Sample sizes used for this experiment.
n_size <- c(2000,10000)

# Collection of model based on the numer of components K of GLLiM.
# K = 1,...,Kmax
Kmax <- 20

##########################################################################################
#                               Saving data of experiments
##########################################################################################

# Matrix contains the selected models based on djump and DDSE criteria from capushe package
# for different sample size over different trials on well-specified (WS) case.
selectedK_Djump_WS <- matrix(0, nrow = numTrial, ncol =  length(n_size))
selectedK_DDSE_WS <- matrix(0, nrow = numTrial, ncol =  length(n_size))

# Matrix contains the selected models based on djump and DDSE criteria from capushe package
# for different sample size over different trials on misspecified (MS) case.
selectedK_Djump_MS <- matrix(0, nrow = numTrial, ncol =  length(n_size))
selectedK_DDSE_MS <- matrix(0, nrow = numTrial, ncol =  length(n_size))

# List of data from capushe package.
dataCapushe_List_WS <- list()
dataCapushe_List_MS <- list()

##########################################################################################
#     Choosing the true parameters for simulated data sets where the true model 
#                         (mixture of components K) equals 2.
##########################################################################################

D <- 1 # dimension of data
pi_True <- c(1/2, 1/2)

####
# Gating network parameters:
# True means of gates
c_True <- matrix(data = NA, nrow = 2, ncol = D)
c_True[1, ] <- c(0.2)
c_True[2, ] <- c(0.8)

# True covariance matrices of Gates
Gamma_True <- array(NA, dim = c(D, D, 2))
Gamma_True[, , 1] <- 0.1*diag(nrow = D)
Gamma_True[, , 2] <- 0.15*diag(nrow = D)

####
# Experts network parameters:
# Well-Specified (WS) case. This model is identical
# with supervised Gaussian Locally Linear Mapping (GLLiM).

# Means of Gaussian Experts
b_True_WS <- c(2, 0)
A_True_WS <- matrix(data = NA, nrow = D, ncol = 2)

A_True_WS[, 1] <- c(-5)
A_True_WS[, 2] <- c(0.1)

# Misspecified (MS) case where we choose parabolic means for Gaussian Experts:

# Means of Gaussian Experts
beta0_MS <- c(1, 0)
beta_MS <- matrix(data = NA, nrow = D, ncol = 2)

beta_MS[, 1] <- c(-6)
beta_MS[, 2] <- c(0)

beta2_MS <- matrix(data = NA, nrow = D, ncol = 2)

beta2_MS[, 1] <- c(1)
beta2_MS[, 2] <- c(-0.4)

# Standard deviations of Gaussian Experts for both cases WS and MS.
sigma_True <- matrix(data = 0.3, nrow = 1, ncol = 2)

##########################################################################################
# Generate simulated data sets (WS and MS) and estimate the parameters of GLoME models
##########################################################################################

# Save the simulation data set in list() format
simData_WS <- list()
simData_MS <- list()

# Sample the data sets in WS and MS cases based on n_size sample sizes.
sample_Data_WS <- NamsGLoME::sample_GLLiM(pi_True, c_True, Gamma_True, b_True_WS, A_True_WS, sigma_True, n_size[length(n_size)], 1)
sample_Data_MS <- NamsGLoME::sample_MoGATEParabola(pi_True, c_True, Gamma_True, beta0_MS, beta_MS, beta2_MS, sigma_True, n_size[length(n_size)], 1)

# Collection of model based on the numer of components K of GLLiM.
# K = 1,...,Kmax

modelComplex_WS <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_WS_JNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_WS_CNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)

modelComplex_MS <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_MS_JNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)
contrast_MS_CNLL <- matrix(0, nrow = length(n_size), ncol = Kmax)


# Initalize matrices tensorized Kullback-Leibler divergence
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




#########################################################################################
# Generate samples and plot samples
########################################################################################
# # Generate samples from GLLiM
#
# echan_WS <- HDMoGATE::sample_GLLiM(pi_True, c_True, Gamma_True, b_True_WS, A_True_WS, sigma_True, n, ny)
# simData_WS[[t]] <- echan_WS
#
#
# # Generate samples from parabolic mean of Experts
#
#
# echan_MS <- HDMoGATE::sample_MoGATEParabola(pi_True, c_True, Gamma_True, beta0_MS, beta_MS, beta2_MS, sigma_True, n, ny)
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
    echan_WS <- HDMoGATE::sample_GLLiM(pi_True, c_True, Gamma_True, b_True_WS, A_True_WS, sigma_True, n_size[n] , ny)
    # simData_WS[[t]] <- echan_WS

    # MS case
    echan_MS <- HDMoGATE::sample_MoGATEParabola(pi_True, c_True, Gamma_True, beta0_MS, beta_MS, beta2_MS, sigma_True, n_size[n], ny)
    # simData_MS[[t]] <- echan_MS

    # For each number of components (k = 1,...,Kmax), we use gllim function to estimate the parameters.
    for (k in 1:Kmax) {
      # Using function gllim from xLLiM package when k > 1
      if (k>1) {

        #k <- 2 # Tesing

        # WS case:
        #while (KL_WS[t,k] == Inf){

        # # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_WS <- HDMoGATE::sample_GLLiM(pi_True, c_True, Gamma_True, b_True_WS, A_True_WS, sigma_True, n, ny)
        # simData_WS[[t]] <- echan_WS

        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        estimaMoEGLLiM_WS <- gllim(t(as.matrix(sample_Data_WS$y[1:n_size[n]])),
                                   t(as.matrix(sample_Data_WS$X[1:n_size[n]])), in_K = k)

        # Using n data points for estimators.
        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        gllim_inverse_WS <- gllim_inverse_dens(t(as.matrix(sample_Data_WS$X[1:n_size[n]])),
                                               estimaMoEGLLiM_WS,t(as.matrix(sample_Data_WS$y[1:n_size[n]])))

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
        # echan_MS <- HDMoGATE::sample_MoGATEParabola(pi_True, c_True, Gamma_True, beta0_MS, beta_MS, beta2_MS, sigma_True, n, ny)
        # simData_MS[[t]] <- echan_MS

        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        estimaMoEGLLiM_MS <- gllim(t(as.matrix(sample_Data_MS$y[1:n_size[n]])),
                                   t(as.matrix(sample_Data_MS$X[1:n_size[n]])), in_K = k)

        # Using n data points for estimators.
        # Estimate the parameters \widehat{s}_{\widehat{m}} using inverse regression from GLLiM
        gllim_inverse_MS <- gllim_inverse_dens(t(as.matrix(sample_Data_MS$X[1:n_size[n]])),
                                               estimaMoEGLLiM_MS,t(as.matrix(sample_Data_MS$y[1:n_size[n]])))

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
        # echan_WS <- HDMoGATE::sample_GLLiM(pi_True, c_True, Gamma_True, b_True_WS, A_True_WS, sigma_True, n, ny)
        # simData_WS[[t]] <- echan_WS

        # Using lm function from R for Fitting Generalized Linear Models
        estimaMoE_WS_lm <- lm(as.matrix(sample_Data_WS$y[1:n_size[n]]) ~ as.matrix(sample_Data_WS$X[1:n_size[n]]))
        modelComplex_WS[n, k] <- attributes(logLik(estimaMoE_WS_lm))$df
        contrast_WS_JNLL[n, k] <- -logLik(estimaMoE_WS_lm)[1]
        contrast_WS_CNLL[n, k] <- -logLik(estimaMoE_WS_lm)[1]

        # # Testing function provides vectors of pdf over all samples from lm()
        # # mypdf and loglik() from R
        # num_obs = 30; sigma_True = 0.3
        # #set.seed(1)
        # epsilon = rnorm(n = num_obs, mean = 0, sd = sigma_True)
        # x_vals = seq(from = 0, to = 10, length.out = num_obs)
        #
        # sim_slr = function(x, beta_0 = 10, beta_1 = 5, sigma_True = 1) {
        #   n = length(x)
        #   epsilon = rnorm(n, mean = 0, sd = sigma_True)
        #   y = beta_0 + beta_1 * x + epsilon
        #   data.frame(predictor = x, response = y)
        # }
        #
        # sim_data = sim_slr(x = x_vals, beta_0 = 5, beta_1 = -2, sigma_True = 3)
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

        CLL_WS_vec <- exp(loggausspdf(t(echan_WS$yny), t(y_hat_WS), matrix(s2_e_WS,nrow = D)))
        estimaMoE_WS_lm_pdf <- sum(log(CLL_WS_vec))

        JKL_WS[n,k,t] <- 1/(n_size[n]*ny*rho)*sum(log(echan_WS$gllim_pdf_vec/((1-rho)*echan_WS$gllim_pdf_vec+rho*CLL_WS_vec)))

        KL_WS[n,k,t] <- 1/(n_size[n]*ny)*(echan_WS$gllim_pdf-estimaMoE_WS_lm_pdf)
        #}

        # MS case:
        #while (KL_MS[t,k] == Inf){

        # # Initalizing the whole data set for estimators and Monte Carlo method
        # echan_MS <- HDMoGATE::sample_MoGATEParabola(pi_True, c_True, Gamma_True, beta0_MS, beta_MS, beta2_MS, sigma_True, n, ny)
        # simData_MS[[t]] <- echan_MS

        # Using lm function from R for Fitting Generalized Linear Models
        estimaMoE_MS_lm <- lm(as.matrix(sample_Data_MS$y[1:n_size[n]]) ~ as.matrix(sample_Data_MS$X[1:n_size[n]]))

        modelComplex_MS[n, k] <- attributes(logLik(estimaMoE_MS_lm))$df
        contrast_MS_JNLL[n, k] <- -logLik(estimaMoE_MS_lm)[1]
        contrast_MS_CNLL[n, k] <- -logLik(estimaMoE_MS_lm)[1]

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        y_hat_MS <- coef(estimaMoE_MS_lm)[1] + coef(estimaMoE_MS_lm)[2]*echan_MS$Xny
        s2_e_MS = sum((echan_MS$yny - y_hat_MS)^2)/(n_size[n]*ny)

        CLL_MS_vec <- exp(loggausspdf(t(echan_MS$yny), t(y_hat_MS), matrix(s2_e_MS,nrow = D)))
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

    dataCapushe_List_WS[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_WS_CNLL)
    slopeHeur_MoE_WS <- capushe(dataCapushe_MoE_WS_CNLL)

    selectedK_Djump_WS[t, n] <- slopeHeur_MoE_WS@Djump@model
    selectedK_DDSE_WS[t, n] <- slopeHeur_MoE_WS@DDSE@model

    KL_Sel_WS_Djump_CNLL[t, n] <- KL_WS[n, as.numeric(selectedK_Djump_WS[t, n]),t]
    KL_Sel_WS_DDSE_CNLL[t, n] <- KL_WS[n, as.numeric(selectedK_DDSE_WS[t, n]),t]

    JKL_Sel_WS_Djump_CNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_Djump_WS[t, n]),t]
    JKL_Sel_WS_DDSE_CNLL[t, n] <- JKL_WS[n, as.numeric(selectedK_DDSE_WS[t, n]),t]

    # MS case:
    dataCapushe_MoE_MS_CNLL <- data.frame(c(1:Kmax), modelComplex_MS[n,], modelComplex_MS[n,], contrast_MS_CNLL[n,])
    names(dataCapushe_MoE_MS_CNLL) <- c("model", "pen", "complexity", "contrast")

    dataCapushe_List_MS[[t]] <- list(n = n_size[n], dataCapushe = dataCapushe_MoE_MS_CNLL)
    slopeHeur_MoE_MS <- capushe(dataCapushe_MoE_MS_CNLL)

    selectedK_Djump_MS[t, n] <- slopeHeur_MoE_MS@Djump@model
    selectedK_DDSE_MS[t, n] <- slopeHeur_MoE_MS@DDSE@model

    KL_Sel_MS_Djump_CNLL[t, n] <- KL_MS[n, as.numeric(selectedK_Djump_MS[t, n]),t]
    KL_Sel_MS_DDSE_CNLL[t, n] <- KL_MS[n, as.numeric(selectedK_DDSE_MS[t, n]),t]

    JKL_Sel_MS_Djump_CNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_Djump_MS[t, n]),t]
    JKL_Sel_MS_DDSE_CNLL[t, n] <- JKL_MS[n, as.numeric(selectedK_DDSE_WS[t, n]),t]

  }
}

end_time <- Sys.time()
runningTime <- end_time - start_time
print(runningTime)
