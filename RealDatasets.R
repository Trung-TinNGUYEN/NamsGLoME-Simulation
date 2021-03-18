  ##########################################################################################################
  #   Non-asymptotic Penalization Criteria for Model Selection in Mixture of Experts (MoE) Models:
  #
  # 1. Performing clustering and regression tasks on a real data set.
  # 2. For instance, we illustrate how to use NamsGLoME on ethanol data set of Brinkman (1981):
  # The data comprises of 88 observations, which represent the relationship between the engine’s concentration of nitrogen
  # oxide (NO) emissions and the equivalence ratio (ER), a measure of the air-ethanol mix, used as a
  # spark-ignition engine fuel in a single-cylinder automobile test (Figures 8a and 8e). Our goal is then
  # to estimate the parameters of a GLoME model, as well as the number of mixture components.
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
  #                                         Load data sets
  ##############################################################################################################
  
  library(SemiPar)
  data(ethanol)
  
  EquivRatio<- matrix(ethanol$E, ncol = length(ethanol$E))
  NO <- matrix(ethanol$NOx, ncol = length(ethanol$NOx))
  
  
  ##############################################################################################################
  #                         Customize hyperparameters for the simulated data sets
  ##############################################################################################################
  
  # Collection of model based on the maximum numer of components K of GLLiM.
  # K = 1,...,Kmax
  Kmax <- 12 # This Kmax is recommended for the Ethanol real data set.
  
  # Number of trials to perform model selection.
  num_trials = 100
  
  ##############################################################################################################
  #                     Perform a specified task based on your seletion:
  # If input_task == 
  # "1". Consider NO as input variable and ER as response.
  # "2". Consider ER as input variable and NO as response.
  # "3". Stop the this experiment.
  # WARNING: RUNNING TIME = 1 hour.
  # Machine: Dell Lattitude 5490 (Intel(R) Core(TM) i5-8250U CPU @ 1.6GHz, 8GB RAM).
  # Hyperparameters: We run two experiments with input_task = 1 and input_task 2 over 100 num_trials = 100, Kmax = 12.
  ##############################################################################################################
  input_task = 0 # Defaut task
  
  while (input_task != 3){
    
    input_task <- readline(prompt = "1: X = NO, 2: X = ER, 3 = Stop: Type input_task = ")
        
    if (input_task == 1){
      ##########################################################################################################
      # Perform the clustering and regression task with NO as input variable and ER as responses
      ##########################################################################################################
      clustering_NO <- NamsGLoME::apply_namsGLoME(X = NO, Y = EquivRatio, num_trials = num_trials, Kmax = Kmax, 
                                                  plot_histogram = TRUE, input_task = 1,
                                                  plot_clustering_ethanol = TRUE, plot_slope_heuristic = TRUE
                                                  )
    } else if (input_task == 2){
      ##########################################################################################################
      # Here, instead of considering the variable NO as the covariate, we use it as the response variable.
      # Then, the resulting clustering, the estimated mean function (black curve)
      #  and mixing probabilities are more easily interpretable.
      ##########################################################################################################
      clustering_EquivRatio <- NamsGLoME::apply_namsGLoME(X = EquivRatio, Y = NO, num_trials = num_trials, Kmax = Kmax,
                                                          plot_histogram = TRUE, plot_clustering_ethanol = TRUE,
                                                          plot_slope_heuristic = TRUE, input_task = 2)
  
    } 
    if ((input_task != 1)&&(input_task != 2)&&(input_task != 3)){
      print("Please run the program again such that the input_task belongs to {1,2,3}!")
    }
    if (input_task == 3){
      print("Thank you for selecting model with GLoME models. Hope to see you again soon!")
    }
  }
  
  
