# NamsGLoME-Simulation
This repository contains all numerical experiments (R code) for "A non-asymptotic approach for model selection via penalization in high-dimensional mixture of experts models", which will appear in the Electronic Journal of Statistics.

# Article
[arXiv: ](https://arxiv.org/abs/2104.02640) A non-asymptotic approach for model selection via penalization in high-dimensional mixture of experts models. Electronic Journal of Statistics, 2022.

## Description of the codes:

    SimulationStudy.R : Run simulations and draw results on simulated data sets.
    
    RealDatasets.R : Use the NamsGLoME package and draw results on the Ethanol data set.
    
    Results: Expected results when running completely the SimulationStudy.R and RealDatasets.R.
    
    NamsGLoME: package used for the SimulationStudy.R and RealDatasets.R, comprising the following functions:
    
        sample_GLLiM.R : Draws an n-sample from a supervised Gaussian locally-linear mapping (GLLiM).
    
        gllim_inverse_dens.R : Evaluate conditional log-likelihood of forward regression model from GLLiM model.
        
        sample_GLoME_parabola.R : Draws an n-sample from a Gaussian-gated localized Mixture of Experts (GLoME) 
                                with parabolic mean experts.
        
        model_true.R : Initalize true parameters model for validating NamsGLoME on simulated data sets:
                            1. Well-specified (WS): Means of Gaussian Experts are first degree polynomials 
                            (linear polynomials).
                            2. Misspecified (MS): Means of Gaussian Experts are second degree polynomials 
                            (quadratic polynomials, parabola).
                            
        simulation.R :      Create numerical schemes and then perform model selection by NamsGLoME:
                            1. We first initalize the whole data sets in two case: well-specified (WS) and 
                            misspecified (MS), which are used for the penalized maximum likelihood estimators 
                            (PMLE) and Monte Carlo method.
                            2. Then, we make use of GLLiM model to estimate the parameters of inverse 
                            regression.
                            3. Using gllim_inverse_dens() leads to the parameters of forward regression.
                            4. Finally, we utilize capushe package to calibrate penalties in the context of 
                            model selection via penalization based on the slope heuristics.
                            
        apply_namsGLoME.R : Apply NamsGLoME to some real applications, for instance, Ethanol data set:
                            1. We make use of GLLiM model to estimate the parameters of inverse regression. 
                            2. Using gllim_inverse_dens() leads to the parameters of forward regression.
                            3. Finally, we utilize capushe package to calibrate penalties in the context of
                            model selection via penalization based on the slope heuristics.
                            
        utils.R : Some complementary functions needed for the packaged NamsGLoME.
    
    
