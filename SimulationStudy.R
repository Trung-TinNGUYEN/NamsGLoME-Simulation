##########################################################################################################
#   Non-asymptotic Penalization Criteria for Model Selection in Mixture of Experts (MoE) Models:
#
# 1. Comparison histogram of selected Gaussian localized MoE (GLoME) models between well-specified (WS)
#     and misspecified (MS) cases using jump and slope criteria over 100 trials.
# 2. Box-plots of the tensorized (Jensen)-Kullback-Leibler ((J)KL) divergence according to the number
#     of mixture components using the jump criterion over 100 trials.
# 3. Error decay of tensorized (J)KL divergence between the true and selected
#   densities based on the jump criterion, represented in a log-log scale, using 30 trials.
#   A free least-square regression with standard error and a regression with slope −1 were
#   added to stress the two different behavior for each graph.
#
# A personal access token (This personal access token is used for install package NamsGLoME) with 
# repo:invite scope was recently added to your account (17 March, 2021).
# 4669acbc1009b7bd4d491d3b43a21f0c9343d6d1

##############################################################################################################

##############################################################################################################
#                                       Install and load NamsGLoME package
##############################################################################################################
  
# Removes all objects from the current workspace (R memory).
rm(list = ls())

# Install the NamsGLoME package to your local machine.
# install.packages("devtools")
# devtools::document("NamsGLoME")
devtools::install("NamsGLoME")
library(NamsGLoME)

##############################################################################################################
#                         Customize hyperparameters for the simulated data sets
##############################################################################################################

# Collection of model based on the maximum numer of components K of GLLiM.
# K = 1,...,Kmax
Kmax <- 20

# Number data points for approximating KL using Monte Carlo method.
ny <- 30
  
# Hyperparameter rho for (Jensen)-Kullback-Leibler divergence.
rho = 1/2

##############################################################################################################
#                     Initalize the true parameters for simulated data sets:
# 1. The default parameters: K_true = 2, D = 1, L = 1.
# 2. Customized parameters: Please modify the function model_true() and simulation()!
##############################################################################################################

GLoME_true <- NamsGLoME::model_true(K_true = 2, D = 1, L = 1)


##############################################################################################################
#                     Perform a specified task based on your seletion:
# If input_task == :
# "1". Clustering deduced from the estimated conditional density of GLoME.
# "2". Histogram and boxplot (HB) numerical experiments.
# "3". Error Decay (HB) numerical experiments.
# "4". Stop the this experiment.
##############################################################################################################
input_task = 0 # Defaut task

while (input_task != 4){
  
  input_task <- readline(prompt = "1 = Clustering, 2 = Histogram, 3 = Box-plot, 4 = Stop: Type input_task = ")
  
  if (input_task == 1){
    ##########################################################################################################
    #                   Clustering deduced from the estimated conditional density of GLoME: 
    # Using by a MAP principle with 2000 data points of example WS and MS. The dash and solid black curves
    # present the true and estimated mean functions.
    ##########################################################################################################
    
    clustering_WS <- NamsGLoME::simulation(num_obs = 2000, GLoME_true = GLoME_true, 
                                           plot_clustering_samples = TRUE)
  } else if (input_task == 2){
    ##########################################################################################################
    #                         Histogram and boxplot (HB) numerical experiments:
    # WARNING: RUNNING TIME = 24 hours.
    # Machine: Dell Lattitude 5490 (Intel(R) Core(TM) i5-8250U CPU @ 1.6GHz, 8GB RAM).
    # Hyperparameters: We run our experiments on 2000 and 10000 data points over 100 trials,
    #  using ny = 30 samples for Monte Carlo method to approximate tensorized Kullback-Leibler divergence.
    ##########################################################################################################
    
    ####
    # Hyperparameters:
    ####
    # Number of trials
    #num_trials_HB <- 100
    num_trials_HB <- 2
    
    # Sample sizes used for HB numerical experiments: 2000 10000
    num_obs_HB <- c(2000, 10000)
    
    # If plot_histogram = TRUE, NamsGLoME::simulation() exports histogram of selected GLoME models 
    # between well-specified (WS) and misspecified (MS) cases using jump and slope criteria over 100 trials.
    
    # If plot_boxplot_KL = TRUE, plot_boxplot_JKL == TRUE, NamsGLoME::simulation() exports box-plots of 
    # the tensorized (Jensen)-Kullback-Leibler divergence according to the number of mixture components
    # using the jump criterion over num_trials.
    
    ####
    # Perform the HB simulation
    ####
    
    simulation_HB <- NamsGLoME::simulation(num_trials_HB, num_obs_HB, GLoME_true, Kmax, ny, rho,
                                        plot_histogram = TRUE, plot_boxplot_KL = TRUE)
    
  } else if (input_task == 3){
    ##########################################################################################################
    #                         Error Decay (HB) numerical experiments:
    # WARNING: RUNNING TIME = 20 hours.
    # Machine: Dell Lattitude 5490 (Intel(R) Core(TM) i5-8250U CPU @ 1.6GHz, 8GB RAM).
    # Hyperparameters: We run our experiments from 1000 to 10079 (11 data points) over 30 trials,
    #  using ny = 30 samples for Monte Carlo method to approximate tensorized Kullback-Leibler divergence.
    ##########################################################################################################
    
    ####
    # Hyperparameters:
    ####
    
    # Number of trials.
    # num_trials_ED <- 30
    num_trials_ED <- 1 # Testing code with 1 trial.
    
    # Sample sizes: 1000  1260  1587  2000  2520  3175  4000  5040  6350  8000 10079.
    # num_obs_ED <- round(seq(1000, 2000, length.out = 11)) # Testing with small samples
    num_obs_ED <- round(1000*2.^((0:10)/3))
    
    # If plot_error_decay_KL = TRUE, NamsGLoME::simulation() exports error decay (ED) of 
    # Tensorized Kullback-Leibler divergence between the true and selected densities 
    # based on the jump criterion, represented in a log-log scale, using num_trials_ED:
    
    ####  
    # Perform the ED simulation   
    ####
    
    # plot_error_decay_KL = TRUE
    simulation_ED <- NamsGLoME::simulation(num_trials_ED, num_obs_ED, GLoME_true, Kmax, ny,
                                           rho, plot_error_decay_KL = TRUE)
    
  }
  if ((input_task != 1)&&(input_task != 2)&&(input_task != 3)&&(input_task != 4)){
    print("Please run the program again such that input_task belongs to {1,2,3,4}!")
  }
}


