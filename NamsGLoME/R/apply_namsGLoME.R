#' @export
apply_namsGLoME = function(X, Y, num_trials = 1, Kmax = 12, plot_histogram = FALSE, save_data = FALSE,
                      plot_clustering_samples = FALSE, model_hat = 2, plot_slope_heuristic = FALSE){

  # %%%%%%%%%%%%%%%%% Non-asymptotic Model Selection in Mixture of Experts Models %%%%%%%%%%%%%%%%%%%%%%
  # %% Author: TrungTin Nguyen (14-03-2021) - tinnguyen0495@gmail.com

  # % Description: For each trial (t), each GLoME model (K):
  #
  # 1. We make use of GLLiM model to estimate the parameters of inverse regression.
  # 2. Using gllim_inverse_dens() leads to the parameters of forward regression.
  # 3. Finally, we utilize capushe package to calibrate penalties in the context of model
  #   selection via penalization based on the slope heuristics.

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # %%%% Input %%%%
  # % - X (DxN): covariate variables.
  # % - Y (LxN): response variables.
  # % - num_trials (1x1): Number of trials for numerical experiments.
  # % - Kmax (1x1): Collection of model based on the numer of components K of GLoME.
  # % - save_data == TRUE: Save the outputs of this numerical experiments
  # %    to local machine on Working directory.
  # % - plot_histogram == TRUE: Histogram of selected GLoME models using jump and slope criteria over num_trials.
  # % - plot_slope_heuristic = TRUE: Plot of the selected model dimension using the jump and slope criteria.
  # % - plot_clustering_samples = TRUE: Perform clustering and regression tasks on input data sets.

  # %%%% Output %%%%

  # % Data from capushe package:
  # %   - model_Djump, model_DDSE (num_trials x 1): Matrix contains the selected models
  # %     based on Djump and DDSE criteria from capushe package over num_trials.
  # %   - data_capushe list(): List of data used for capushe package.
  # %   - complexity, contrast (num_obs x Kmax): Collection of model complexity values and
  # %     minimum contrast value for each model based on the numer of mixture components K of GLoME:
  # %     K = 1,...,Kmax.
  # % - model_hat: after running num_trials on the input data sets,
  # %   the model with the highest frequency on the histogram of selected model is used.
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
import::from(capushe, capushe, Djump, DDSE)

# Create Elegant Data Visualisations Using the Grammar of Graphics:
# A system for 'declaratively' creating graphics, based on "The Grammar of Graphics".
# You provide the data, tell 'ggplot2' how to map variables to aesthetics,
# what graphical primitives to use, and it takes care of the details.
#install.packages("ggplot2")
import::from(ggplot2, ggplot, geom_line, geom_smooth, ggtitle, aes, scale_x_log10, scale_y_log10,
             scale_linetype_manual, scale_colour_manual, xlab, ylab, theme, margin, geom_point,
             scale_color_manual, scale_shape_manual, labs)

# Provides a number of user-level functions to work with "grid" graphics,
# notably to arrange multiple grid-based plots on a page, and draw tables.
import::from(gridExtra, grid.arrange) # Side-by-side plots with ggplot2.

# R package plot3D (Soetaert 2013b) contains functions for plotting multi-dimensional
# data. Many functions are derived from the persp function, other functions start from the
# image or contour function.
import::from(plot3D, persp3D)

# fields: Tools for Spatial Data
# For curve, surface and function fitting with an emphasis on splines,
# spatial data, geostatistics, and spatial statistics.
import::from(fields, image.plot)

# gridBase: Integration of base and grid graphics
# Integration of base and grid graphics
import::from(grid, viewport, unit)

import::from(pracma, Mode)

##########################################################################################################
#                             Extract information from the input data.
##########################################################################################################

num_obs <- c(ncol(X))
D <- nrow(X)
L <- nrow(Y)

##########################################################################################################
#                               Saving data of experiments
##########################################################################################################

# Matrix contains the selected models based on Djump and DDSE criteria from
# capushe package over different trials
model_Djump <- model_DDSE <- matrix(0, nrow = num_trials, ncol =  length(num_obs))

# List of data used for capushe package.
data_capushe <- list()

# Collection of model complexity values and minimum contrast value for each model
# based on the numer of mixture components K of GLoME: K = 1,...,Kmax
complexity <- contrast <- matrix(0, nrow = length(num_obs), ncol = Kmax)


##########################################################################################################
#         Perform model selection and estimation the parameters of GLoME models over num_Trials
##########################################################################################################

# Save running time for this experiments
start_time <- Sys.time()

####
# For each trial (t), each GLoME model (K), we first use GLLiM model
# to estimate the parameters of inverse regression.
# This leads to the parameters of forward regression.
# Next, we make use of capushe package to calibrate penalties in the context of model
# selection via penalization based on the slope heuristics.

for (t in 1:num_trials) {
  for (n in 1:length(num_obs)){
    for (K in 1:Kmax) {
      ####
      # K > 1: Using gllim function from xLLiM package:
      if (K>1) {

        # Estimate the inverese parameters \widehat{s}_{\widehat{m}} using GLLiM.
        estimate_GLoME <- gllim(Y, X, in_K = K)

        # Estimate the parameters \widehat{s}^*_{\widehat{m}} using inverse regression trick.
        forward_model_GLoME <- gllim_inverse_dens(X, estimate_GLoME, Y)

        complexity[n, K] <- estimate_GLoME$nbpar
        contrast[n, K] <- forward_model_GLoME$CNLL
        ####
        # K = 1: Using lm function from stats package in R to estimate MLE for linear regression:
      } else {
        ####
        # Using lm function from stats package in R to estimate MLE for linear regression
        estimate_lm <- lm(t(Y) ~ t(X))
        complexity[n, K] <- attributes(logLik(estimate_lm))$df
        contrast[n, K] <- -logLik(estimate_lm)[1]

        }
    }
    ####
    # Using capushe package to select model:
    ####

    # WS case:
    data_capushe_df <- data.frame(c(1:Kmax), complexity[n,], complexity[n,], contrast[n,])
    names(data_capushe_df) <- c("model", "pen", "complexity", "contrast")

    data_capushe[[t]] <- list(n = num_obs[n], dataCapushe = data_capushe_df)
    slope_heuristics <- capushe(data_capushe_df)

    model_Djump[t, n] <- slope_heuristics@Djump@model
    model_DDSE[t, n] <- slope_heuristics@DDSE@model
  }
}

##########################################################################################################
# The selected model based on the highest frequency over num_trials.
##########################################################################################################
model_hat_Djump_mode <- Mode(matrix(as.numeric(model_Djump)))
model_hat_Djump_index <- match(as.character(model_hat_Djump_mode), model_Djump)

model_hat_DDSE_mode <- Mode(matrix(as.numeric(model_Djump)))
model_hat_DDSE_index <- match(as.character(model_hat_DDSE_mode), model_Djump)

##########################################################################################################
# Slope heuristics: plot of the selected model dimension based on the highest frequency.
##########################################################################################################

if (plot_slope_heuristic == TRUE){

  pdf("Plot_Slope_Heuristics.pdf",  width = 15, height = 10)
  op <- par(mfrow = c(2, 1))
  # Jump criterion:
  plot(Djump(as.data.frame(data_capushe[[model_hat_Djump_index]]$dataCapushe)), newwindow=FALSE)
  # Slope criterion:
  plot(DDSE(as.data.frame(data_capushe[[model_hat_DDSE_index]]$dataCapushe)), newwindow=FALSE)
  op <- par(mfrow = c(1, 1))
  dev.off()
}

##########################################################################################################
#  Plot histograms of selected models for WS and MS cases using jump andslope criteria over num_trials.
##########################################################################################################

if (plot_histogram == TRUE){

  pdf("Histograms_Selected_K_100_Trials_Djump_DDSE_MS.pdf",  width = 8.27, height = 11.69)
  op <- par(mfrow = c(2, 1))
  hist(as.numeric(model_Djump[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(a) Selected model using jump criterion")

  hist(as.numeric(model_DDSE[,1]), breaks = c(1:Kmax), col = "blue", xlab = "Selected number of classes",
       ylab = "Empirical Probability", freq = 0,main = "(e) Selected model using slope criterion")

  op <- par(mfrow = c(1, 1))
  dev.off()
}

##########################################################################################################
#                       Save capushe data sets to files in local machine.
##########################################################################################################

# If save_data = TRUE, the capushe data sets are exported to files in local machine, otherwise, skip this step.
# Default value: save_data  = FALSE

if (save_data == TRUE){
  saveRDS(data_capushe, file = "data_capushe.rds")
  saveRDS(model_Djump, file = "model_Djump.rds")
  saveRDS(model_DDSE, file = "model_DDSE.rds")
}

##########################################################################################################
#                   Clustering deduced from the estimated conditional density of GLoME:
# Using by a MAP principle with 2000 data points of example WS and MS. The dash and solid black curves
# present the true and estimated mean functions. With our default simulated data sets, in WS case, we
# select model with K = 2. For MS case, model with K = 4 is utilized to perform clustering and regression.
##########################################################################################################
if ((plot_clustering_samples == TRUE)&&((D > 1)|| (L > 1))){
  print("The current version does not support to plot the high-dimension data! We update it in the furture work!")
}

if ((plot_clustering_samples == TRUE)&&(D == 1) && (L == 1)){
  #########
  #WS case
  #########
  pdf("Clustering_2000_realization.pdf", width = 11.69, height = 8.27)
  # ggplot first
  ggplot(ethanol, aes(x= NOx , y= E)) +
    geom_point(size = 3, colour = "blue", shape = 1) + labs(x = "NO", y = "Equivalence Ratio")
  ####
  # Figures of clustering deduced from the estimated conditional density by a MAP principle
  # Maximum a posteriori probability for the laten variable Z, p(Z_i=k|X_i,Y_i, \widehat{m})
  # Find the maximum position for each row of a matrix, breaking ties at random.
  ####

  # Choose the best initializtion based on its likelihood values
  # using emgm over n_init_emgm iteration.
  n_init_emgm <- 100
  init_emgm_llh <- -Inf
  init_emgm <- list()
  inital_df <- rbind(EquivRatio, NO)

  for (i in 1:n_init_emgm){
    init_emgm_temp <- emgm(inital_df, init=selectedK_Djump_Ethanol_CNLL_mode)

    if (init_emgm_temp$llh > init_emgm_llh){
      init_emgm <- init_emgm_temp
      init_emgm_llh <- init_emgm_temp$llh
    }
  }

  # Estimate the inverese parameters \widehat{s}_{\widehat{m}} using GLLiM.
  inverse_model_hat <- gllim(Y, X, in_K = model_hat_Djump_mode, in_r = init_emgm, maxiter = 1000)

  # Estimate the parameters \widehat{s}^*_{\widehat{m}} using inverse regression trick.
  forward_model_GLoME <- gllim_inverse_dens(X, inverse_model_hat, Y)


  estimate_class <- as.factor(max.col(inverse_model_hat$r))
  # Visualize the resulted clustering on WS data set.
  estimate_class_df <- data.frame(X, Y, estimate_class, t(forward_model_GLoME$x_exp))
  names(estimate_class_df) <- c('X', 'Y', 'Class','PostMeans')

  # Calculate the whole estimated mean function

  # Creating vector of colors for each class.
  color_class_clustering <- rainbow(model_hat)
  shape_class_clustering <- c(1:model_hat)

  clustering_2000_realization <- ggplot() +
    geom_point(data = estimate_class_df,
               aes(x=  X, y= Y, shape=Class, color=Class)) +
    scale_color_manual(values = color_class_clustering) + scale_shape_manual(values=shape_class_clustering)+
    labs(x = " X" , y =  "Y")

  for (k in 1:model_hat){
    gllim_As<-
      forward_model_GLoME$As[,,k]

    gllim_bs<-
      forward_model_GLoME$bs[,k]

    data_k <- estimate_class_df[which(estimate_class==k),]

    data_k_x <- seq(from = min(data_k$X), to = max(data_k$X), length.out = length_plot)
    data_k_y <- gllim_As*data_k_x + gllim_bs

    data_subClass <- data.frame(data_k_x = data_k_x, data_k_y = data_k_y)
    clustering_2000_realization <- clustering_2000_realization +
      geom_line(data = data_subClass, aes(x= data_k_x , y= data_k_y), color = color_class_clustering[k])
  }
  # Calculate the whole estimated mean function
  data_x_estimate <- seq(from = min(estimate_class_df$X),
                            to = max(estimate_class_df$X), length.out = length_plot)
  data_y_estimate <- seq(from = min(estimate_class_df$Y),
                            to = max(estimate_class_df$Y), length.out = length_plot)

  data_y_estimate_gllim_dens <- gllim_inverse_dens(matrix(data_x_estimate,ncol = length_plot),
                                                      inverse_model_hat,
                                                      matrix(data_y_estimate, ncol = length_plot))

  data_y_estimate_estiMeans <- data_y_estimate_gllim_dens$x_exp

  data_y_estimate_estiMeans <- data.frame(data_x_estimate = data_x_estimate,
                                             data_y_estimate_estiMeans = t(data_y_estimate_estiMeans))

  clustering_2000_realization <- clustering_2000_realization +
    geom_line(data = data_y_estimate_estiMeans, aes(x= data_x_estimate , y= data_y_estimate_estiMeans), color = "black")+
    theme(legend.position = "none") + ggtitle("(b) Clustering by GLoME in WS case.")

  # ####
  # # Plot the estimated posterior of mixture proportions.
  # ####
  #
  # # Creating vector of colors for each class.
  # color_class_posterior <- rainbow(model_hat)
  # #color_class_posterior <- c("#00FFFFFF","#FF0000FF" )
  # shape_class_posterior <- c(1:model_hat)
  #
  # clustering_2000_posterior <- ggplot()
  # for (k in 1:model_hat){
  #
  #   data_Posterior <- estimate_class_df
  #   data_k_x <- seq(from = min(data_Posterior$X), to = max(data_Posterior$X), length.out = length_plot)
  #   data_k_y <- seq(from = min(data_Posterior$Y), to = max(data_Posterior$Y), length.out = length_plot)
  #
  #   # Calculate the estimated posterior for the mixing proportion for each mixture components
  #   data_k_y_gllim_dens <- gllim_inverse_dens(matrix(data_k_x,ncol = length_plot),
  #                                             forward_model_GLoME,  matrix(data_k_y, ncol = length_plot))
  #   data_k_y_posterior <- data_k_y_gllim_dens$alpha[,k]
  #
  #   data_estimPosterior <- data.frame(data_k_x = data_k_x, data_k_y_posterior = data_k_y_posterior)
  #   clustering_2000_posterior <- clustering_2000_posterior +
  #     geom_line(data = data_estimPosterior, aes(x= data_k_x , y= data_k_y_posterior),
  #               color = color_class_posterior[k]) +
  #     labs(x = "X", y = " Mixing probabilities")
  # }
  # clustering_2000_posterior <- clustering_2000_posterior +
  #   theme(legend.position = "none") + ggtitle("(d) Gating network probabilities.")

  # create an apporpriate viewport.  Modify the dimensions and coordinates as needed
  vp.11 <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"),  just=c("left","top"),
                    y = 1, x = 0)

  vp.12<- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"),  just=c("left","top"),
                   y = 1, x = 0.5)

  vp.22 <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"),  just=c("left","top"),
                    y = 0.5, x = 0.5)

  # Plot your base graphics
  par(mfrow=c(2,2))

  length_plot_3D  <- 200

  data_x_estimate_3D <- seq(from = min(X), to = max(X), length.out = length_plot_3D)
  data_y_estimate_3D <- seq(from = min(Y), to = max(Y), length.out = length_plot_3D)

  forward_model_GLoME_3D <- matrix(NA, nrow = length_plot_3D, ncol = length_plot_3D)
  for (i in 1:length_plot_3D){
    for (j in 1:length_plot_3D){
      forward_model_GLoME_3D[i,j] <-
        gllim_inverse_dens(matrix(data_x_estimate_3D[i]),
                           inverse_model_hat, matrix(data_y_estimate_3D[j]))$CLL_vec
    }
  }

  plot(0,type='n',axes=FALSE,ann=FALSE)
  plot(0,type='n',axes=FALSE,ann=FALSE)
  forward_estimate_2D <- image.plot(data_x_estimate_3D, data_y_estimate_3D, forward_model_GLoME_3D,
                                       xlab = "X", ylab = "Y",col = topo.colors(length_plot_3D^2),
                                       main = "(c) 2D view of the resulting conditional density
                                            with the 2 regression components")
  # Plot the ggplot using the print command
  print(forward_model_true_2000_realization, vp = vp.11)
  #print(clustering_2000_realization, vp = vp.12)
  #print(clustering_2000_posterior, vp = vp.22)

  dev.off()


}

end_time <- Sys.time()
running_time <- end_time - start_time

###########################################################################################################
# Output list
###########################################################################################################

output <- list(running_time = running_time, model_Djump = model_Djump, model_DDSE = model_DDSE,
                 data_capushe = data_capushe, complexity = complexity, contrast = contrast)

return(output)

}
