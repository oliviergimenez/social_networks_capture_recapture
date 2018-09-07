# Inferring animal social networks from capture-recapture data 

This repository contains `R` codes to reproduce the analyses in a paper submitted that aims to infer animal social networks from capture-recapture data. There are two directories. The directory `dolphin_case_study` contains codes for the real case study on Commerson's dolphin data and is used, among others, to produce Table 1 and Figures 1 and 2. The directory `simulations` is used to asses the bias in parameter estimates.

We used a Bayesian approach with MCMC simulations to implement our model. We used non-informative uniform U[0,1] prior distributions for all parameters. Convergence was assessed using the Brooks-Gelman-Rubin statistic. Mixing was checked visually by inspecting the chains. 
