##########################################################################################################
#       1. Comparison HISTOGRAMS of selected GLoME models between well-specified (WS)
#           and misspecified (MS) cases using jump and slope criteria over 100 trials.
#       2. BOX-PLOTS of the tensorized Kullback-Leibler divergence according to the number
#             of mixture components using the jump criterion over 100 trials.
##########################################################################################################
# RUNNING TIME on Dell Lattitude 5490 (Intel(R) Core(TM) i5-8250U CPU @ 1.6GHz, 8GB RAM).
# We run our experiments on 2000 and 10000 data points over 100 trials,
#  using 30 samples for Monte Carlo method to approximate tensorized Kullback-Leibler divergence:
# 1. Histograms: around 10 hours.
# 2. Box-plot: around 18 hours.
# Total: 28 hours.
##########################################################################################################

# Removes all objects from the current workspace (R memory).
#rm(list = ls())

##########################################################################################################
#                                         Load packages
##########################################################################################################

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

##########################################################################################################
#                           Setting parameters for experiments
##########################################################################################################

# Number of trials for numerical experiments.
num_trials <- 100

# Sample sizes used for this experiment.
num_obs <- c(2000,10000)

# Collection of model based on the numer of components K of GLLiM.
# K = 1,...,Kmax
Kmax <- 20

# Number data points for approximating KL using Monte Carlo method.
ny <- 30

# Hyperparameter rho for (Jensen)-Kullback-Leibler divergence.
rho = 1/2

##########################################################################################################
#                               Saving data of experiments
##########################################################################################################

# Matrix contains the selected models based on djump and DDSE criteria from capushe package
# for different sample size over different trials on well-specified (WS) case.
model_Djump_WS <- model_DDSE_WS <- matrix(0, nrow = num_trials, ncol =  length(num_obs))

# Matrix contains the selected models based on djump and DDSE criteria from capushe package
# for different sample size over different trials on misspecified (MS) case.
model_Djump_MS <- model_DDSE_MS <- matrix(0, nrow = num_trials, ncol =  length(num_obs))

# List of data used for capushe package.
data_capushe_WS <- data_capushe_MS <- list()

####
# Collection of model complexity values and minimum contrast value for each model
# based on the numer of mixture components K of GLoME: K = 1,...,Kmax
complexity_WS <- contrast_WS <- matrix(0, nrow = length(num_obs), ncol = Kmax)
complexity_MS <- contrast_MS <- matrix(0, nrow = length(num_obs), ncol = Kmax)

####
# Initalize tensorized (Jensen)-Kullback-Leibler divergence between true and estimated densities:

# WS case:
# KL over all trials, models and num_obs samples
KL_WS <- Inf*array(1, dim = c(length(num_obs), Kmax, num_trials))

# KL for a selected model over all trials on num_obs samples using Djump and DDSE.
KL_model_hat_Djump_WS <- KL_model_hat_DDSE_WS <- matrix(0, num_trials, length(num_obs))

# JKL over all trials, models and num_obs samples
JKL_WS <- Inf*array(1, dim = c(length(num_obs), Kmax, num_trials))

# JKL for a selected model over all trials on num_obs samples using Djump and DDSE.
JKL_model_hat_Djump_WS <- JKL_model_hat_DDSE_WS <- matrix(0, num_trials, length(num_obs))

# MS case:
# KL over all trials, models and num_obs samples
KL_MS <- Inf*array(1, dim = c(length(num_obs), Kmax, num_trials))

# KL for the selected model over all trials on num_obs samples using Djump and DDSE.
KL_model_hat_Djump_MS <- KL_model_hat_DDSE_MS <- matrix(0, num_trials, length(num_obs))

# JKL over all trials, models and num_obs samples
JKL_MS <- Inf*array(1, dim = c(length(num_obs), Kmax, num_trials))

# KL for the selected model over all trials on num_obs samples using Djump and DDSE.
JKL_model_hat_Djump_MS <- JKL_model_hat_DDSE_MS <- matrix(0, num_trials, length(num_obs))

##########################################################################################################
#     Choosing the true parameters for simulated data sets where the true model
#                         (mixture of components K) equals 2.
##########################################################################################################

D <- 1 # dimension of data
pi_True <- c(1/2, 1/2)

####
# Gating network parameters:
# True means of gates
c_true <- matrix(data = NA, nrow = 2, ncol = D)
c_true[1, ] <- c(0.2)
c_true[2, ] <- c(0.8)

# True covariance matrices of Gates
Gamma_true <- array(NA, dim = c(D, D, 2))
Gamma_true[, , 1] <- 0.1*diag(nrow = D)
Gamma_true[, , 2] <- 0.15*diag(nrow = D)

####
# Experts network parameters:
# Well-Specified (WS) case. This model is identical
# with supervised Gaussian Locally Linear Mapping (GLLiM).

# Means of Gaussian Experts
b_true_WS <- c(2, 0)
A_true_WS <- matrix(data = NA, nrow = D, ncol = 2)

A_true_WS[, 1] <- c(-5)
A_true_WS[, 2] <- c(0.1)

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
sigma_true <- matrix(data = 0.3, nrow = 1, ncol = 2)

##########################################################################################################
# Generate simulated data sets (WS and MS) and estimate the parameters of GLoME models
##########################################################################################################

# Save running time for this experiments
start_time <- Sys.time()

####
# For each trial (t), each sample size (n), each GLoME model (K), we first use GLLiM model
# to estimate the parameters of inverse regression.
# This leads to the parameters of forward regression.
# Next, we make use of capushe package to calibrate penalties in the context of model
# selection via penalization based on the slope heuristics.

for (t in 1:num_trials) {
  for (n in 1:length(num_obs)){

    # Initalize the whole data sets used for penalized maximum likelihood estimators (PMLE)
    # and Monte Carlo method:
    # WS case
    sample_data_WS <- NamsGLoME::sample_GLLiM(pi_True, c_true, Gamma_true, b_true_WS,
                                              A_true_WS, sigma_true, num_obs[n] , ny)

    # MS case
    sample_data_MS <- NamsGLoME::sample_GLoME_parabola(pi_True, c_true, Gamma_true, beta0_MS,
                                                       beta_MS, beta2_MS, sigma_true, num_obs[n], ny)

    for (K in 1:Kmax) {
      ####
      # K > 1: Using gllim function from xLLiM package:
      if (K>1) {

        # WS case:
        # Using n data points for estimators.
        # Estimate the inverese parameters \widehat{s}_{\widehat{m}} using GLLiM.
        estimate_GLoME_WS <- gllim(t(as.matrix(sample_data_WS$y[1:num_obs[n]])),
                                   t(as.matrix(sample_data_WS$X[1:num_obs[n]])), in_K = K)

        # Estimate the parameters \widehat{s}^*_{\widehat{m}} using inverse regression trick.
        gllim_inverse_WS <- gllim_inverse_dens(t(as.matrix(sample_data_WS$X[1:num_obs[n]])),
                                               estimate_GLoME_WS,t(as.matrix(sample_data_WS$y[1:num_obs[n]])))

        complexity_WS[n, K] <- estimate_GLoME_WS$nbpar
        contrast_WS[n, K] <- gllim_inverse_WS$CNLL

        ####
        # Approximate (Jensen)-Kullback-Leibler divergence KL(dens1,dens2)
        # (or JKL(dens1,dens2)) by Monte Carlo method using ny*n points

        # KL case:
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} KL(dens1(./X_i),dens2(./X_i))
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/ny sum_{j=1..ny}
        # log(dens1(Y_ij/x_i)/dens2(Y_ij/X_i))
        # with Y_ij is sampled from dens1(./X_i).

        # Using n*ny data points for Monte  Carlo.
        KL_WS[n,K,t] <- sample_data_WS$gllim_pdf+gllim_inverse_dens(t(sample_data_WS$Xny),estimate_GLoME_WS,t(sample_data_WS$yny))$CNLL
        KL_WS[n,K,t] <- KL_WS[n,K,t]/(num_obs[n]*ny)
        #}

        # JKL case: Given rho such that 0 < rho < 1,
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho
        # KL(dens1(./X_i), (1-rho)*dens1(./X_i) + rho*dens2(./X_i))
        # JKL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/rho 1/ny sum_{j=1..ny}
        # log(dens1(Y_ij/X_i)/(1-rho)*dens1(./X_i) + rho*dens2(./X_i))
        # with Y_ij is sampled from dens1(./Y_i)

        # Vector of pdf values over all samples \hat{s}^*_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        CLL_WS_vec <- gllim_inverse_dens(t(sample_data_WS$Xny),estimate_GLoME_WS,t(sample_data_WS$yny))$CLL_vec

        JKL_WS[n,K,t] <- 1/(num_obs[n]*ny*rho)*sum(log(sample_data_WS$gllim_pdf_vec/
                                                        ((1-rho)*sample_data_WS$gllim_pdf_vec+rho*CLL_WS_vec)))

        # MS case:
        # Using n data points for estimators.
        # Estimate the inverese parameters \widehat{s}_{\widehat{m}} using GLLiM.
        estimate_GLoME_MS <- gllim(t(as.matrix(sample_data_MS$y[1:num_obs[n]])),
                                   t(as.matrix(sample_data_MS$X[1:num_obs[n]])), in_K = K)

        # Estimate the parameters \widehat{s}^*_{\widehat{m}} using inverse regression trick.
        gllim_inverse_MS <- gllim_inverse_dens(t(as.matrix(sample_data_MS$X[1:num_obs[n]])),
                                               estimate_GLoME_MS,t(as.matrix(sample_data_MS$y[1:num_obs[n]])))

        complexity_MS[n, K] <- estimate_GLoME_MS$nbpar
        contrast_MS[n, K] <- gllim_inverse_MS$CNLL

        # Approximate (Jensen)-Kullback-Leibler divergence KL(dens1,dens2)
        # (or JKL(dens1,dens2)) by Monte Carlo method using ny*n points

        # KL case:
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} KL(dens1(./x_i),dens2(./x_i))
        # KL(dens1,dens2) ~= 1/n sum_{i=1..n} 1/ny sum_{j=1..ny}
        # log(dens1(y_ij/x_i)/dens2(y_ij/x_i))
        # with y_ij is sampled from dens1(./x_i)

        # Using n*ny data points for Monte  Carlo.
        KL_MS[n,K,t] <- sample_data_MS$moe2_pdf+gllim_inverse_dens(t(sample_data_MS$Xny),estimate_GLoME_MS,
                                                                   t(sample_data_MS$yny))$CNLL
        KL_MS[n,K,t] <- KL_MS[n,K,t]/(num_obs[n]*ny)
        #}

        # Vector of pdf values over all sample \hat{s}^*_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        CLL_MS_vec <- gllim_inverse_dens(t(sample_data_MS$Xny),estimate_GLoME_MS,t(sample_data_MS$yny))$CLL_vec

        JKL_MS[n,K,t] <- 1/(num_obs[n]*ny*rho)*sum(log(sample_data_MS$moe2_pdf_vec/
                                                        ((1-rho)*sample_data_MS$moe2_pdf_vec+rho*CLL_MS_vec)))
      ####
      # K = 1: Using lm function from stats package in R to estimate MLE for linear regression:
      } else {
        ####
        # WS case:
        # Using lm function from stats package in R to estimate MLE for linear regression
        estimate_lm_WS <- lm(as.matrix(sample_data_WS$y[1:num_obs[n]]) ~ as.matrix(sample_data_WS$X[1:num_obs[n]]))
        complexity_WS[n, K] <- attributes(logLik(estimate_lm_WS))$df
        contrast_WS[n, K] <- -logLik(estimate_lm_WS)[1]

        # Approximate (Jensen)-Kullback-Leibler divergence KL(dens1,dens2)
        # (or JKL(dens1,dens2)) by Monte Carlo method using ny*n points:

        # Vector of pdf values over all sample \hat{s}^*_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        y_hat_WS <- coef(estimate_lm_WS)[1] + coef(estimate_lm_WS)[2]*sample_data_WS$Xny
        s2_e_WS = sum((sample_data_WS$yny - y_hat_WS)^2)/(num_obs[n]*ny)

        CLL_WS_vec <- exp(loggausspdf(t(sample_data_WS$yny), t(y_hat_WS), matrix(s2_e_WS,nrow = D)))
        estimaMoE_WS_lm_pdf <- sum(log(CLL_WS_vec))

        KL_WS[n,K,t] <- 1/(num_obs[n]*ny)*(sample_data_WS$gllim_pdf-estimaMoE_WS_lm_pdf)
        JKL_WS[n,K,t] <- 1/(num_obs[n]*ny*rho)*sum(log(sample_data_WS$gllim_pdf_vec/
                                                         ((1-rho)*sample_data_WS$gllim_pdf_vec+rho*CLL_WS_vec)))

        ####
        # MS case:
        # Using lm function from stats package in R to estimate MLE for linear regression.
        estimaMoE_MS_lm <- lm(as.matrix(sample_data_MS$y[1:num_obs[n]]) ~ as.matrix(sample_data_MS$X[1:num_obs[n]]))

        complexity_MS[n, K] <- attributes(logLik(estimaMoE_MS_lm))$df
        contrast_MS[n, K] <- -logLik(estimaMoE_MS_lm)[1]

        # Approximate (Jensen)-Kullback-Leibler divergence KL(dens1,dens2)
        # (or JKL(dens1,dens2)) by Monte Carlo method using ny*n points:

        # Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
        y_hat_MS <- coef(estimaMoE_MS_lm)[1] + coef(estimaMoE_MS_lm)[2]*sample_data_MS$Xny
        s2_e_MS = sum((sample_data_MS$yny - y_hat_MS)^2)/(num_obs[n]*ny)

        CLL_MS_vec <- exp(loggausspdf(t(sample_data_MS$yny), t(y_hat_MS), matrix(s2_e_MS,nrow = D)))
        estimaMoE_MS_lm_pdf <- sum(log(CLL_MS_vec))

        KL_MS[n, K, t] <- 1/(num_obs[n]*ny)*(sample_data_MS$moe2_pdf-estimaMoE_MS_lm_pdf)
        JKL_MS[n, K, t] <- 1/(num_obs[n]*ny*rho)*sum(log(sample_data_MS$moe2_pdf_vec/
                                                           ((1-rho)*sample_data_MS$moe2_pdf_vec+rho*CLL_MS_vec)))
      }
    }
    ####
    # Using capushe package to select model:
    ####

    # WS case:
    data_capushe_WS_df <- data.frame(c(1:Kmax), complexity_WS[n,], complexity_WS[n,], contrast_WS[n,])
    names(data_capushe_WS_df) <- c("model", "pen", "complexity", "contrast")

    data_capushe_WS[[t]] <- list(n = num_obs[n], dataCapushe = data_capushe_WS_df)
    slope_heuristics_WS <- capushe(data_capushe_WS_df)

    model_Djump_WS[t, n] <- slope_heuristics_WS@Djump@model
    model_DDSE_WS[t, n] <- slope_heuristics_WS@DDSE@model

    KL_model_hat_Djump_WS[t, n] <- KL_WS[n, as.numeric(model_Djump_WS[t, n]),t]
    KL_model_hat_DDSE_WS[t, n] <- KL_WS[n, as.numeric(model_DDSE_WS[t, n]),t]

    JKL_model_hat_Djump_WS[t, n] <- JKL_WS[n, as.numeric(model_Djump_WS[t, n]),t]
    JKL_model_hat_DDSE_WS[t, n] <- JKL_WS[n, as.numeric(model_DDSE_WS[t, n]),t]

    # MS case:
    data_capushe_MS_df <- data.frame(c(1:Kmax), complexity_MS[n,], complexity_MS[n,], contrast_MS[n,])
    names(data_capushe_MS_df) <- c("model", "pen", "complexity", "contrast")

    data_capushe_MS[[t]] <- list(n = num_obs[n], dataCapushe = data_capushe_MS_df)
    slope_heuristics_MoE_MS <- capushe(data_capushe_MS_df)

    model_Djump_MS[t, n] <- slope_heuristics_MoE_MS@Djump@model
    model_DDSE_MS[t, n] <- slope_heuristics_MoE_MS@DDSE@model

    KL_model_hat_Djump_MS[t, n] <- KL_MS[n, as.numeric(model_Djump_MS[t, n]),t]
    KL_model_hat_DDSE_MS[t, n] <- KL_MS[n, as.numeric(model_DDSE_MS[t, n]),t]

    JKL_model_hat_Djump_MS[t, n] <- JKL_MS[n, as.numeric(model_Djump_MS[t, n]),t]
    JKL_model_hat_DDSE_MS[t, n] <- JKL_MS[n, as.numeric(model_DDSE_MS[t, n]),t]

  }
}

end_time <- Sys.time()
running_time <- end_time - start_time
print(running_time)

##########################################################################################################
# Save simulated WS and MS data sets and its related terms to files in local machine.
##########################################################################################################

# If save_data = TRUE, the simulated WS and MS data sets and its related terms
# is exported to files in local machine, otherwise, skip this step.
# Default value: save_data <-  FALSE

save_data <-  FALSE

if (save_data == TRUE){

  saveRDS(KL_WS, file = "KL_WS.rds")
  saveRDS(KL_model_hat_Djump_WS, file = "KL_model_hat_Djump_WS.rds")
  saveRDS(KL_model_hat_DDSE_WS, file = "KL_model_hat_DDSE_WS.rds")

  saveRDS(KL_MS, file = "KL_MS.rds")
  saveRDS(KL_model_hat_Djump_MS, file = "KL_model_hat_Djump_MS.rds")
  saveRDS(KL_model_hat_DDSE_MS, file = "KL_model_hat_DDSE_MS.rds")

  saveRDS(JKL_WS, file = "JKL_WS.rds")
  saveRDS(JKL_model_hat_Djump_WS, file = "JKL_model_hat_Djump_WS.rds")
  saveRDS(JKL_model_hat_DDSE_WS, file = "JKL_model_hat_DDSE_WS.rds")

  saveRDS(JKL_MS, file = "JKL_MS.rds")
  saveRDS(JKL_model_hat_Djump_MS, file = "JKL_model_hat_Djump_MS.rds")
  saveRDS(JKL_model_hat_DDSE_MS, file = "JKL_model_hat_DDSE_MS.rds")

  saveRDS(data_capushe_WS, file = "data_capushe_WS.rds")
  saveRDS(model_Djump_WS, file = "model_Djump_WS.rds")
  saveRDS(model_DDSE_WS, file = "model_DDSE_WS.rds")

  saveRDS(data_capushe_MS, file = "data_capushe_MS.rds")
  saveRDS(model_Djump_MS, file = "model_Djump_MS.rds")
  saveRDS(model_DDSE_MS, file = "model_DDSE_MS.rds")
}

##########################################################################################################
#             Plot histograms of selected models for WS and MS cases
#                 using jump andslope criteria over 100 trials.
##########################################################################################################

  pdf("Histograms_Selected_K_100_Trials_Djump_DDSE_WS_MS.pdf",  width = 8.27, height = 11.69)
  op <- par(mfrow = c(4, 2))
  hist(as.numeric(model_Djump_WS[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(a) 2000 WS data points using jump criterion")

  hist(as.numeric(model_Djump_WS[,2]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(b) 10000 WS data points using jump criterion")

  hist(as.numeric(model_Djump_MS[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(c) 2000 MS data points using jump criterion")

  hist(as.numeric(model_Djump_MS[,2]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(d) 10000 MS data points using jump criterion")

  hist(as.numeric(model_DDSE_WS[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(e) 2000 WS data points using slope criterion")

  hist(as.numeric(model_DDSE_WS[,2]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(f) 10000 WS data points using slope criterion")

  hist(as.numeric(model_DDSE_MS[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(g) 2000 MS data points using slope criterion")

  hist(as.numeric(model_DDSE_MS[,2]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(h) 10000 MS data points using slope criterion")
  op <- par(mfrow = c(1, 1))
  dev.off()

##########################################################################################################
# Box-plot of the tensorized (Jensen)-Kullback-Leibler divergence according to the number of
#           mixture components using the jump criterion over 100 trials.
##########################################################################################################

###############################
# Kullback-Leibler divergence.
###############################
pdf("Boxplot_KL_100_Trials_Djump_WS_MS.pdf",  width = 11.69, height = 8.27)
op <- par(mfrow = c(2, 2))

####
# Example WS with 2000 data points.
####

# Create a data frame that combines all interested boxplots for WS case using jump criterion
KL_WS_df_2000 <- data.frame(t(KL_WS[1,,]), KL_model_hat_Djump_WS[,1])
names(KL_WS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(KL_WS_df_2000, border = "blue", ylim = c(min(KL_WS[1,-1,]),max(KL_WS[1,-1,])))

# Add some lines
lines(1:Kmax, complexity_WS[1,1:Kmax]/(2*num_obs[1]),
      lty = "dashed", col = "black") # asymptotic E[KL]

lines(1:Kmax, colMeans(t(KL_WS[1,,])), lty = "solid", col = "blue") # asymptotic E[KL]

lines(1:21, mean(KL_model_hat_Djump_WS[,1])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[KL]

# Make a legend for lines
legend(8, 0.98*(min(KL_WS[1,-1,])+max(KL_WS[1,-1,])),
       legend = c("asymptotic E[KL]", "empirical E[KL]", "E[KL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(a) Example WS with 2000 data points")

####
# Example WS with 10000 data points.
####
KL_WS_df_2000 <- data.frame(t(KL_WS[2,,]), KL_model_hat_Djump_WS[,2])
names(KL_WS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(KL_WS_df_2000, border = "blue", ylim = c(min(KL_WS[2,-1,]),max(KL_WS[2,-1,])))

# Add some lines
lines(1:Kmax, complexity_WS[2,1:Kmax]/(2*num_obs[2]),
      lty = "dashed", col = "black") # asymptotic E[KL]

lines(1:Kmax, colMeans(t(KL_WS[2,,])), lty = "solid", col = "blue") # asymptotic E[KL]

lines(1:21, mean(KL_model_hat_Djump_WS[,2])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[KL]

# Make a legend for lines
legend(8, 0.96*(min(KL_WS[2,-1,])+max(KL_WS[2,-1,])),
       legend = c("asymptotic E[KL]", "empirical E[KL]", "E[KL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(b) Example WS with 10000 data points")

####
# Example MS with 2000 data points.
####
KL_MS_df_2000 <- data.frame(t(KL_MS[1,,]), KL_model_hat_Djump_MS[,1])
names(KL_MS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(KL_MS_df_2000, border = "blue", ylim = c(min(KL_MS[1,-1,]),max(KL_MS[1,-1,])))

# Add some lines
lines(1:Kmax, complexity_MS[1,1:Kmax]/(2*num_obs[1]),
      lty = "dashed", col = "black") # asymptotic E[KL]

lines(1:Kmax, colMeans(t(KL_MS[1,,])), lty = "solid", col = "blue") # asymptotic E[KL]

lines(1:21, mean(KL_model_hat_Djump_MS[,1])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[KL]

# Make a legend for lines
legend(8, 0.98*(min(KL_MS[1,-1,])+max(KL_MS[1,-1,])),
       legend = c("asymptotic E[KL]", "empirical E[KL]", "E[KL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(c) Example MS with 2000 data points")

####
# Example MS with 10000 data points.
####
KL_MS_df_2000 <- data.frame(t(KL_MS[2,,]), KL_model_hat_Djump_MS[,2])
names(KL_MS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(KL_MS_df_2000, border = "blue", ylim = c(min(KL_MS[2,-1,]),max(KL_MS[2,-1,])))

# Add some lines
lines(1:Kmax, complexity_MS[2,1:Kmax]/(2*num_obs[2]),
      lty = "dashed", col = "black") # asymptotic E[KL]

lines(1:Kmax, colMeans(t(KL_MS[2,,])), lty = "solid", col = "blue") # asymptotic E[KL]

lines(1:21, mean(KL_model_hat_Djump_MS[,2])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[KL]

# Make a legend for lines
legend(8, 0.96*(min(KL_MS[2,-1,])+max(KL_MS[2,-1,])),
       legend = c("asymptotic E[KL]", "empirical E[KL]", "E[KL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(d) Example MS with 10000 data points")

op <- par(mfrow = c(1, 1))
dev.off()

#####################################
# Jensen-Kullback-Leibler divergence.
#####################################
pdf("Boxplot_JKL_100_Trials_Djump_WS_MS.pdf",  width = 11.69, height = 8.27)
op <- par(mfrow = c(2, 2))

####
# Example WS with 2000 data points.
####

# Create a data frame that combines all interested boxplots for WS case using jump criterion
JKL_WS_df_2000 <- data.frame(t(JKL_WS[1,,]), JKL_model_hat_Djump_WS[,1])
names(JKL_WS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(JKL_WS_df_2000, border = "blue", ylim = c(min(JKL_WS[1,-1,]),max(JKL_WS[1,-1,])))

# Add some lines
lines(1:Kmax, complexity_WS[1,1:Kmax]/(2*num_obs[1]),
      lty = "dashed", col = "black") # asymptotic E[JKL]

lines(1:Kmax, colMeans(t(JKL_WS[1,,])), lty = "solid", col = "blue") # asymptotic E[JKL]

lines(1:21, mean(JKL_model_hat_Djump_WS[,1])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[JKL]

# Make a legend for lines
legend(8, 0.98*(min(JKL_WS[1,-1,])+max(JKL_WS[1,-1,])),
       legend = c("asymptotic E[JKL]", "empirical E[JKL]", "E[JKL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(a) Example WS with 2000 data points")

####
# Example WS with 10000 data points.
####
JKL_WS_df_2000 <- data.frame(t(JKL_WS[2,,]), JKL_model_hat_Djump_WS[,2])
names(JKL_WS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(JKL_WS_df_2000, border = "blue", ylim = c(min(JKL_WS[2,-1,]),max(JKL_WS[2,-1,])))

# Add some lines
lines(1:Kmax, complexity_WS[2,1:Kmax]/(2*num_obs[2]),
      lty = "dashed", col = "black") # asymptotic E[JKL]

lines(1:Kmax, colMeans(t(JKL_WS[2,,])), lty = "solid", col = "blue") # asymptotic E[JKL]

lines(1:21, mean(JKL_model_hat_Djump_WS[,2])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[JKL]

# Make a legend for lines
legend(8, 0.96*(min(JKL_WS[2,-1,])+max(JKL_WS[2,-1,])),
       legend = c("asymptotic E[JKL]", "empirical E[JKL]", "E[JKL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(b) Example WS with 10000 data points")

####
# Example MS with 2000 data points.
####
JKL_MS_df_2000 <- data.frame(t(JKL_MS[1,,]), JKL_model_hat_Djump_MS[,1])
names(JKL_MS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(JKL_MS_df_2000, border = "blue", ylim = c(min(JKL_MS[1,-1,]),max(JKL_MS[1,-1,])))

# Add some lines
lines(1:Kmax, complexity_MS[1,1:Kmax]/(2*num_obs[1]),
      lty = "dashed", col = "black") # asymptotic E[JKL]

lines(1:Kmax, colMeans(t(JKL_MS[1,,])), lty = "solid", col = "blue") # asymptotic E[JKL]

lines(1:21, mean(JKL_model_hat_Djump_MS[,1])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[JKL]

# Make a legend for lines
legend(8, 0.98*(min(JKL_MS[1,-1,])+max(JKL_MS[1,-1,])),
       legend = c("asymptotic E[JKL]", "empirical E[JKL]", "E[JKL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(c) Example MS with 2000 data points")

####
# Example MS with 10000 data points.
####
JKL_MS_df_2000 <- data.frame(t(JKL_MS[2,,]), JKL_model_hat_Djump_MS[,2])
names(JKL_MS_df_2000) <- c(as.character(c(1:Kmax)),"SK")

# Plot a boxplot: set the y range in boxplot graph without the first column.
boxplot(JKL_MS_df_2000, border = "blue", ylim = c(min(JKL_MS[2,-1,]),max(JKL_MS[2,-1,])))

# Add some lines
lines(1:Kmax, complexity_MS[2,1:Kmax]/(2*num_obs[2]),
      lty = "dashed", col = "black") # asymptotic E[JKL]

lines(1:Kmax, colMeans(t(JKL_MS[2,,])), lty = "solid", col = "blue") # asymptotic E[JKL]

lines(1:21, mean(JKL_model_hat_Djump_MS[,2])*matrix(1,1,21),
      col = "green", pch = 4, type = "o") # asymptotic E[JKL]

# Make a legend for lines
legend(8, 0.96*(min(JKL_MS[2,-1,])+max(JKL_MS[2,-1,])),
       legend = c("asymptotic E[JKL]", "empirical E[JKL]", "E[JKL] of the selected K"),
       col = c("black", "blue","green"),
       lty = c("dashed","solid","solid"),
       pch = c(NA, NA, 4))
title("(d) Example MS with 10000 data points")

op <- par(mfrow = c(1, 1))
dev.off()


##########################################################################################################
# Plot error decay of Tensorized Kullback-Leibler divergence between the true and selected
#  densities based on the jump criterion, represented in a log-log scale, using 30 trials.
##########################################################################################################

