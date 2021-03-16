#' @export
model_true = function(K_true = 2, D = 1, L = 1){
  # %%%%%%%%%%%%%%%%% Generate true parameters for GLoME models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %% Author: TrungTin Nguyen (14-03-2021) - tinnguyen0495@gmail.com

  # % Description:
  #     We initalize true parameters model for validating "Non-asymptotic penalization criteria for
  #     model selection in Gaussian localized mixture of experts models (GLoME)" in two cases:
  # 1. Well-specified (WS): Means of Gaussian Experts are first degree polynomials (linear polynomials).
  # 2. Misspecified (MS): Means of Gaussian Experts are second degree polynomials
  # %                                   (quadratic polynomials, parabola).
  # %   By defaut, we choose the following true parameters.
  # %  However, users can change them to arbitrary values.
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # %%%% Input %%%%
  # % - K_true (1x1): True model (number of mixture components).
  # % - D (1x1): Dimensional of covariates X(Dx1).
  # % - L (1x1): Dimensional of response variables Y(Lx1).
  # %%%% Output %%%%
  # % - GLoME_true list(): structure contains true parameters for simulation. By default,
  # %   it comprises the following components:
  # %   - GLoME_true$K_true (1x1): True model (number of mixture components).
  # %   - GLoME_true$D (1x1): Dimensional of covariates X(Dx1).
  # %   - GLoME_true$L (1x1): Dimensional of response variables Y(Lx1).
  # %   Gaussian gating functions:
  # %     - GLoME_true$pi_true (1xK_true): Prior on mixture proportion to generate covariates.
  # %     - GLoME_true$c_true (K_truexD): True means.
  # %     - GLoME_true$Gamma_true (DxDxK_true): True covariance matrices.
  # %   Means of Gaussian Experts: WS case with First degree polynomials (linear polynomials).
  # %     - GLoME_true$b_true_WS (LxK_true): True intercept vectors.
  # %     - GLoME_true$A_true_WS (LxDxK_true): True slope coefficient matrices.
  # %   Means of Gaussian Experts: MS case with second degree polynomials
  # %                                   (quadratic polynomials, parabola).
  # %     - GLoME_true$beta0_MS (LxK_true): True coefficient vectors of the first term.
  # %     - GLoME_true$beta_MS (LxDxK_true): True coefficient matrices of the second term.
  # %     - GLoME_true$beta2_MS (LxDxK_true): True coefficient matrices of the third term.
  # %   Covariance matrices of Gaussian Experts:
  # %     - GLoME_true$sigma_true (LxL): True standard deviations of Gaussian Experts.
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##########################################################################################################
# Create the true parameters for simulated data sets based on the GLoME_true structure.
##########################################################################################################
# Structure contains true parameters for simulation.
GLoME_true <- list()

GLoME_true$K_true <- K_true  # (1x1)
GLoME_true$D <- D  # (1x1)
GLoME_true$L <-L # (1x1)

#### Prior on mixture proportion to generate covariates from K_true mixture components.
GLoME_true$pi_true <- matrix(1/K_true, nrow = 1, ncol = K_true)

if ((K_true != 2) || (L != 1) || (D != 1)){
  print('Warning: The default parameters: K_true = 2, D = 1, L = 1.
        Please modify the following parts so that it is consistent with your customized parameters!')
  break
}

####
# Gating network parameters:
# True means of gates
c_true <- matrix(data = NA, nrow = K_true, ncol = D)
c_true[1, ] <- c(0.2)
c_true[2, ] <- c(0.8)
GLoME_true$c_true <- c_true

# True covariance matrices of Gates
Gamma_true <- array(NA, dim = c(D, D, K_true))
Gamma_true[, , 1] <- 0.1*diag(nrow = D)
Gamma_true[, , 2] <- 0.15*diag(nrow = D)
GLoME_true$Gamma_true <- Gamma_true

####
# Experts network parameters:
# Well-Specified (WS) case. This model is identical
# with supervised Gaussian Locally Linear Mapping (GLLiM).

# Means of Gaussian Experts
b_true_WS <- matrix(data = c(2, 0), nrow = K_true, ncol = D) # (LxK_true)
GLoME_true$b_true_WS <- b_true_WS

A_true_WS <- matrix(data = NA, nrow = D, ncol = K_true) # (LxDxK_true)
A_true_WS[, 1] <- c(-5)
A_true_WS[, 2] <- c(0.1)
GLoME_true$A_true_WS <- A_true_WS
# Misspecified (MS) case where we choose parabolic means for Gaussian Experts:

# Means of Gaussian Experts
beta0_MS <- matrix(data = c(1, 0), nrow = K_true, ncol = D)
GLoME_true$beta0_MS <- beta0_MS

beta_MS <- matrix(data = NA, nrow = D, ncol = K_true)
beta_MS[, 1] <- c(-6)
beta_MS[, 2] <- c(0)
GLoME_true$beta_MS <- beta_MS

beta2_MS <- matrix(data = NA, nrow = D, ncol = K_true)
beta2_MS[, 1] <- c(1)
beta2_MS[, 2] <- c(-0.4)
GLoME_true$beta2_MS <- beta2_MS

# Standard deviations of Gaussian Experts for both cases WS and MS.
sigma_true <- matrix(data = 0.3, nrow = 1, ncol = K_true)
GLoME_true$sigma_true <- sigma_true

return(GLoME_true)
}
