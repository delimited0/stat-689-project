----
title: STAT 689 Project
author: Patrick Ding, James Dole, Naveed Merchant
----

# Proposal

Our group wants to look further in to MCMC methods other than just metropolis hastings and gibbs sampling. We want to look into Hamiltonian MC and find out how this method compares to others on different models. By this we mean if we choose different models, do some methods converge more quickly then others or do some methods simply fail to successfully approximate the posterior distribution at all.

## Goals

1. We need to find a way to actually evaluate how well a MCMC method performs. This is one of the things we will research and apply. Afterwards we will explain how these methods compare with each other using this evaluation, and explain exactly what this evaluation is.

2. We will implement HMC and metropolis hastings in R for a variety of models and compare their performance. 

3. The models we want to fit are: 
  * The conjugate normal-normal model
  * A multimodal normal mixture model
  * A logistic regression model
  * If we have more time we intend to examine how HMC works with other models.
  
4. We will evaluate these models and inference algorithms on simulated data. 

5. We will also try to apply these methods to some real data.

## Aspirational Goals

1. We can explore extentions to HMC through other implementations like Stan and Edward. 

2. If we have more time we will look into Langevian dynamics which we can implement through Edward.

