#' @export
simulation = function(num_trials, num_obs, GLoME_true, Kmax = 20, ny = 30, rho = 1/2, save_data = FALSE,
                      plot_histogram = FALSE, plot_boxplot_KL = FALSE, plot_boxplot_JKL = FALSE,
                      plot_error_decay_KL = FALSE, plot_error_decay_JKL = FALSE,
                      t_constant_WS = 3, t_constant_MS = 20, plot_clustering_2000_samples == FALSE){

  # %%%%%%%%%%%%%%%%% Non-asymptotic Model Selection in Mixture of Experts Models %%%%%%%%%%%%%%%%%%%%%%
  # %% Author: TrungTin Nguyen (14-03-2021) - tinnguyen0495@gmail.com

  # % Description: For each trial (t), each sample (n), each GLoME model (K):

  # 1. We first initalize the whole data sets in two case: well-specified (WS) and misspecified (MS),
  #   which are used for the penalized maximum likelihood estimators (PMLE) and Monte Carlo method.

  # 2. Then, we make use of GLLiM model to estimate the parameters of inverse regression.

  # 3. Using gllim_inverse_dens() leads to the parameters of forward regression.

  # 4. Finally, we utilize capushe package to calibrate penalties in the context of model
  #   selection via penalization based on the slope heuristics.

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # %%%% Input %%%%

  # % - num_trials (1x1): Number of trials for numerical experiments.
  # % - num_obs (1xlength(num_obs)): Sample sizes used for this experiment.
  # % - Kmax (1x1): Collection of model based on the numer of components K of GLoME.
  # % - ny (1x1): Number data points for approximating KL using Monte Carlo method.
  # % - rho in (0,1) (1x1): Hyperparameter rho for (Jensen)-Kullback-Leibler ((J)KL) divergence.

  # % - GLoME_true list(): structure contains true parameters for simulation.
  # %   By default, it comprises the following components:
  # %   - GLoME_true$K_true (1x1): True model (number of mixture components).
  # %   - GLoME_true$D (1x1): Dimensional of covariates X(Dx1).
  # %   - GLoME_true$L (1x1): Dimensional of response variables Y(Lx1).

  # %   Gaussian gating functions:
  # %     - GLoME_true$pi_true (1xK_true): Prior on mixture.
  # %     - GLoME_true$c_true (K_truexD): True means.
  # %     - GLoME_true$Gamma_true (DxDxK_true): True covariance matrices.
  # %   Means of Gaussian Experts: WS case with First degree polynomials (linear polynomials).
  # %     - GLoME_true$b_true_WS (LxK_true): True intercept vectors.
  # %     - GLoME_true$A_true_WS (LxDxK_true): True slope coefficient matrices.

  # %   Means of Gaussian Experts:
  # %   MS case with second degree polynomials (quadratic polynomials, parabola).
  # %     - GLoME_true$beta0_MS (LxK_true): True coefficient vectors of the first term.
  # %     - GLoME_true$beta_MS (LxDxK_true): True coefficient matrices of the second term.
  # %     - GLoME_true$beta2_MS (LxDxK_true): True coefficient matrices of the third term.
  # %   Covariance matrices of Gaussian Experts:
  # %     - GLoME_true$sigma_true (LxL): True standard deviations of Gaussian Experts.

  # % - save_data == TRUE: Save the outputs of this numerical experiments
  # %    to local machine on Working directory.
  # % - plot_histogram == TRUE: Comparison histogram of selected GLoME models between well-specified (WS)
  # %   and misspecified (MS) cases using jump and slope criteria over num_trials.
  # % - plot_boxplot_KL == TRUE, plot_boxplot_JKL == TRUE: Box-plots of the tensorized
  # %   (J)KL divergence according to the number of mixture components
  # %   using the jump criterion over num_trials.
  # % - plot_error_decay_KL = TRUE, plot_error_decay_JKL = TRUE: Error decay of tensorized
  # %   (J)KL divergence between the true and selected  densities based on the jump criterion,
  # %   represented in a log-log scale, using 30 trials.
  # % - t_constant_WS = 3, t_constant_MS = 20: Default values for contanst in error decays.

  # %%%% Output %%%%

  # % WS case:
  # %   - model_Djump_WS, model_DDSE_WS (num_trials x num_obs): Matrix contains the selected models
  # %     based on djump and DDSE criteria from capushe package for different sample size over
  # %     different trials on well-specified (WS) case.
  # %   - data_capushe_WS, data_capushe_MS list(): List of data used for capushe package.
  # %   - complexity_WS, contrast_WS (num_obs x Kmax): Collection of model complexity values and
  # %     minimum contrast value for each model based on the numer of mixture components K of GLoME:
  # %     K = 1,...,Kmax.

  # %   - KL_WS (1x num_obs x num_trials): KL over all trials, models and num_obs samples.
  # %   - KL_model_hat_Djump_WS, KL_model_hat_DDSE_WS (num_obs x num_trials): KL for a selected model.
  # %   - JKL_WS (1x num_obs x num_trials): KL over all trials, models and num_obs samples.
  # %   - JKL_model_hat_Djump_WS, JKL_model_hat_DDSE_WS (num_obs x num_trials): JKL for a selected model.

  # % MS case:
  # %   - model_Djump_MS, model_DDSE_MS (num_trials x num_obs): Matrix contains the selected models
  # %     based on djump and DDSE criteria from capushe package for different sample size over
  # %     different trials on well-specified (MS) case.
  # %   - data_capushe_MS, data_capushe_MS list(): List of data used for capushe package.
  # %   - complexity_MS, contrast_MS (num_obs x Kmax): Collection of model complexity values and
  # %     minimum contrast value for each model based on the numer of mixture components K of GLoME:
  # %     K = 1,...,Kmax.

  # %   - KL_MS (1x num_obs x num_trials): KL over all trials, models and num_obs samples.
  # %   - KL_model_hat_Djump_MS, KL_model_hat_DDSE_MS (num_obs x num_trials): KL for a selected model.
  # %   - JKL_MS (1x num_obs x num_trials): KL over all trials, models and num_obs samples.
  # %   - JKL_model_hat_Djump_MS, JKL_model_hat_DDSE_MS (num_obs x num_trials): JKL for a selected model.
  # % - running_time (s/mins/hours): total running time for the function NamsGLoME_simulation().
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##########################################################################################################
#                           Import some functions used for this package.
##########################################################################################################
#

# High Dimensional Locally-Linear Mapping:
# Provides a tool for non linear mapping (non linear regression)
# using a mixture of regression model and an inverse regression strategy.
#install.packages("xLLiM")
import::from(xLLiM, gllim)

# CAlibrating Penalities Using Slope HEuristics (CAPUSHE):
# The capushe function proposes two algorithms based on the slope heuristics
# to calibrate penalties in the context of model selection via penalization.
#install.packages("capushe")
import::from(capushe, capushe)

# Create Elegant Data Visualisations Using the Grammar of Graphics:
# A system for 'declaratively' creating graphics, based on "The Grammar of Graphics".
# You provide the data, tell 'ggplot2' how to map variables to aesthetics,
# what graphical primitives to use, and it takes care of the details.
#install.packages("ggplot2")
import::from(ggplot2, ggplot, geom_line, geom_smooth, ggtitle, aes, scale_x_log10, scale_y_log10,
             scale_linetype_manual, scale_colour_manual, xlab, ylab, theme, margin)

# Provides a number of user-level functions to work with "grid" graphics,
# notably to arrange multiple grid-based plots on a page, and draw tables.
import::from(gridExtra, grid.arrange) # Side-by-side plots with ggplot2.


##########################################################################################################
# Create the true parameters for simulated data sets based on the GLoME_true structure.
##########################################################################################################

K_true <- GLoME_true$K_true # (1x1)
D <- GLoME_true$D # (1x1)
L <- GLoME_true$L # (1x1)

if ((K_true != 2) || (L != 1) || (D != 1)){
  print('Warning: The default parameters: K_true = 2, D = 1, L = 1.
        Please modify the following parts so that it is consistent with your customized parameters!')
  break
}

# %   Gaussian gating functions:
pi_true <- GLoME_true$pi_true # (1x1)
c_true <- GLoME_true$c_true # (K_truexD)
Gamma_true <- GLoME_true$Gamma_true # (DxDxK_true)

# Means of Gaussian Experts: WS case with First degree polynomials (linear polynomials).
b_true_WS <- GLoME_true$b_true_WS # (LxK_true)
A_true_WS <- GLoME_true$A_true_WS # (LxDxK_true)

# Means of Gaussian Experts: MS case with second degree polynomials (quadratic polynomials, parabola).
beta0_MS <- GLoME_true$beta0_MS # (LxK_true)
beta_MS <- GLoME_true$beta_MS # (LxDxK_true)
beta2_MS <- GLoME_true$beta2_MS # (LxDxK_true)

# Covariance matrices of Gaussian Experts:
sigma_true <- GLoME_true$sigma_true # (LxL)

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
    sample_data_WS <- NamsGLoME::sample_GLLiM(pi_true, c_true, Gamma_true, b_true_WS,
                                              A_true_WS, sigma_true, num_obs[n] , ny)

    # MS case
    sample_data_MS <- NamsGLoME::sample_GLoME_parabola(pi_true, c_true, Gamma_true, beta0_MS,
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


##########################################################################################################
#  Plot histograms of selected models for WS and MS cases using jump andslope criteria over num_trials.
##########################################################################################################

if (plot_histogram == TRUE){
  if (length(num_obs) < 2){
    print('Default length of num_obs is at least 2. Please modify the hyperparameters!')
    break
  }
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
}

##########################################################################################################
# Box-plot of the tensorized (Jensen)-Kullback-Leibler divergence according to the number of
#           mixture components using the jump criterion over 100 trials.
##########################################################################################################

#####################################
# Kullback-Leibler divergence.
#####################################

if (plot_boxplot_KL == TRUE){
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

}

#####################################
# Jensen-Kullback-Leibler divergence.
#####################################

if (plot_boxplot_JKL == TRUE){
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

}

##########################################################################################################
# Error decay of tensorized (J)KL divergence between the true and selected
#   densities based on the jump criterion, represented in a log-log scale, using 30 trials.
#   A free least-square regression with standard error and a regression with slope −1 were
#   added to stress the two different behavior for each graph.
##########################################################################################################

#####################################
# Kullback-Leibler divergence.
#####################################

if (plot_error_decay_KL == TRUE){

  pdf("KL_model_hat_Djump_WS_MS_error_decay.pdf",  width = 11.69, height = 8.27)

  # WS case:
  KL_model_hat_Djump_WS_df <- data.frame(num_obs = num_obs, E_KL = colMeans(KL_model_hat_Djump_WS))
  KL_model_hat_Djump_WS_asym_df <- data.frame(num_obs = num_obs, E_KL_asym = t_constant_WS/num_obs)

  KL_model_hat_Djump_WS_ggplot <- ggplot() +
    geom_line(data = KL_model_hat_Djump_WS_df,
              aes(x = num_obs, y = E_KL, col = "E[KL]", linetype = "E[KL]"))+
    scale_x_log10() + scale_y_log10() +
    geom_smooth(data = KL_model_hat_Djump_WS_df, method='lm', formula = y ~ x,
                aes(x = num_obs, y = E_KL,col = "linear regression of E[KL]",linetype = "linear regression of E[KL]")) +
    geom_line(data = KL_model_hat_Djump_WS_asym_df,
              aes(x = num_obs, y = E_KL_asym, col = "n -> t/n", linetype = "n -> t/n"))+
    xlab("Sample size")+ylab("Mean of Kullback Leibler distance")+
    scale_linetype_manual(name = "Error decay",
                          values =c("E[KL]" = "solid","linear regression of E[KL]"= "dashed","n -> t/n" = "dotdash"))+
    scale_colour_manual(name = "Error decay",
                        values =c("E[KL]" = "blue","linear regression of E[KL]" = "red", "n -> t/n"="black"))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(10, 10, 10, 10)) +
      ggtitle("(a) WS Example.")

  # MS case:
  KL_model_hat_Djump_MS_df <- data.frame(num_obs = num_obs, E_KL = colMeans(KL_model_hat_Djump_MS))
  KL_model_hat_Djump_MS_asym_df <- data.frame(num_obs = num_obs, E_KL_asym = t_constant_MS/num_obs)

  KL_model_hat_Djump_MS_ggplot <- ggplot() +
    geom_line(data = KL_model_hat_Djump_MS_df,
              aes(x = num_obs, y = E_KL, col = "E[KL]", linetype = "E[KL]"))+
    scale_x_log10() + scale_y_log10() +
    geom_smooth(data = KL_model_hat_Djump_MS_df, method='lm', formula = y ~ x,
                aes(x = num_obs, y = E_KL,col = "linear regression of E[KL]",linetype = "linear regression of E[KL]")) +
    geom_line(data = KL_model_hat_Djump_MS_asym_df,
              aes(x = num_obs, y = E_KL_asym, col = "n -> t/n", linetype = "n -> t/n"))+
    xlab("Sample size")+ylab("Mean of Kullback Leibler distance")+
    scale_linetype_manual(name = "Error decay",
                          values =c("E[KL]" = "solid","linear regression of E[KL]"= "dashed","n -> t/n" = "dotdash"))+
    scale_colour_manual(name = "Error decay",
                        values =c("E[KL]" = "blue","linear regression of E[KL]" = "red", "n -> t/n"="black"))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(10, 10, 10, 10)) +
      ggtitle("(b) MS Example.")

  grid.arrange(KL_model_hat_Djump_WS_ggplot, KL_model_hat_Djump_MS_ggplot, nrow = 1)
  dev.off()
}

#####################################
# Jensen-Kullback-Leibler divergence.
#####################################

if (plot_error_decay_JKL == TRUE){

  pdf("JKL_model_hat_Djump_WS_MS_error_decay.pdf",  width = 11.69, height = 8.27)

  # WS case:
  JKL_model_hat_Djump_WS_df <- data.frame(num_obs = num_obs, E_JKL = colMeans(JKL_model_hat_Djump_WS))
  JKL_model_hat_Djump_WS_asym_df <- data.frame(num_obs = num_obs, E_JKL_asym = t_constant_WS/num_obs)

  JKL_model_hat_Djump_WS_ggplot <- ggplot() +
    geom_line(data = JKL_model_hat_Djump_WS_df,
              aes(x = num_obs, y = E_JKL, col = "E[JKL]", linetype = "E[JKL]"))+
    scale_x_log10() + scale_y_log10() +
    geom_smooth(data = JKL_model_hat_Djump_WS_df, method='lm', formula = y ~ x,
                aes(x = num_obs, y = E_JKL,col = "linear regression of E[JKL]",linetype = "linear regression of E[JKL]")) +
    geom_line(data = JKL_model_hat_Djump_WS_asym_df,
              aes(x = num_obs, y = E_JKL_asym, col = "n -> t/n", linetype = "n -> t/n"))+
    xlab("Sample size")+ylab("Mean of Kullback Leibler distance")+
    scale_linetype_manual(name = "Error decay",
                          values =c("E[JKL]" = "solid","linear regression of E[JKL]"= "dashed","n -> t/n" = "dotdash"))+
    scale_colour_manual(name = "Error decay",
                        values =c("E[JKL]" = "blue","linear regression of E[JKL]" = "red", "n -> t/n"="black"))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(10, 10, 10, 10)) +
    ggtitle("(a) WS Example.")

  # MS case:
  JKL_model_hat_Djump_MS_df <- data.frame(num_obs = num_obs, E_JKL = colMeans(JKL_model_hat_Djump_MS))
  JKL_model_hat_Djump_MS_asym_df <- data.frame(num_obs = num_obs, E_JKL_asym = t_constant_MS/num_obs)

  JKL_model_hat_Djump_MS_ggplot <- ggplot() +
    geom_line(data = JKL_model_hat_Djump_MS_df,
              aes(x = num_obs, y = E_JKL, col = "E[JKL]", linetype = "E[JKL]"))+
    scale_x_log10() + scale_y_log10() +
    geom_smooth(data = JKL_model_hat_Djump_MS_df, method='lm', formula = y ~ x,
                aes(x = num_obs, y = E_JKL,col = "linear regression of E[JKL]",linetype = "linear regression of E[JKL]")) +
    geom_line(data = JKL_model_hat_Djump_MS_asym_df,
              aes(x = num_obs, y = E_JKL_asym, col = "n -> t/n", linetype = "n -> t/n"))+
    xlab("Sample size")+ylab("Mean of Kullback Leibler distance")+
    scale_linetype_manual(name = "Error decay",
                          values =c("E[JKL]" = "solid","linear regression of E[JKL]"= "dashed","n -> t/n" = "dotdash"))+
    scale_colour_manual(name = "Error decay",
                        values =c("E[JKL]" = "blue","linear regression of E[JKL]" = "red", "n -> t/n"="black"))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(10, 10, 10, 10)) +
    ggtitle("(b) MS Example.")

  grid.arrange(JKL_model_hat_Djump_WS_ggplot, JKL_model_hat_Djump_MS_ggplot, nrow = 1)
  dev.off()
}

##########################################################################################################
# Save simulated MS and MS data sets and its related terms to files in local machine.
##########################################################################################################

# If save_data = TRUE, the simulated WS and MS data sets and its related terms
# is exported to files in local machine, otherwise, skip this step.
# Default value: save_data  = FALSE

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

  KL_model_hat_Djump_WS_slope_lm <- lm(log(colMeans(KL_model_hat_Djump_WS)) ~ log(num_obs))
  saveRDS(KL_model_hat_Djump_WS_slope_lm$coefficients[2], file = "KL_model_hat_Djump_WS_slope_lm.rds")
  KL_model_hat_Djump_MS_slope_lm <- lm(log(colMeans(KL_model_hat_Djump_MS)) ~ log(num_obs))
  saveRDS(KL_model_hat_Djump_MS_slope_lm$coefficients[2], file = "KL_model_hat_Djump_MS_slope_lm.rds")

}

end_time <- Sys.time()
running_time <- end_time - start_time

# Output list
output <- list(running_time = running_time, model_Djump_WS = model_Djump_WS, model_DDSE_WS = model_DDSE_WS,
             model_Djump_MS = model_Djump_MS, model_DDSE_MS = model_DDSE_MS, data_capushe_WS = data_capushe_WS,
             data_capushe_MS = data_capushe_MS, complexity_WS = complexity_WS, contrast_WS = contrast_WS,
             complexity_MS = complexity_MS, contrast_MS = contrast_MS, KL_WS = KL_WS,
             KL_model_hat_Djump_WS = KL_model_hat_Djump_WS, KL_model_hat_DDSE_WS = KL_model_hat_DDSE_WS,
             JKL_WS = JKL_WS, JKL_model_hat_Djump_WS = JKL_model_hat_Djump_WS, JKL_model_hat_DDSE_WS = JKL_model_hat_DDSE_WS,
             KL_MS = KL_MS, KL_model_hat_Djump_MS = KL_model_hat_Djump_MS, KL_model_hat_DDSE_MS = KL_model_hat_DDSE_MS,
             JKL_MS = JKL_MS, JKL_model_hat_Djump_MS = JKL_model_hat_Djump_MS, JKL_model_hat_DDSE_MS = JKL_model_hat_DDSE_MS,
             KL_model_hat_Djump_WS_slope_lm =  lm(log(colMeans(KL_model_hat_Djump_WS)) ~ log(num_obs))$coefficients[2],
             KL_model_hat_Djump_MS_slope_lm = lm(log(colMeans(KL_model_hat_Djump_MS)) ~ log(num_obs))$coefficients[2])

return(output)
}
