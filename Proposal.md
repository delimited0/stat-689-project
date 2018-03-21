# Proposal

1. Our group wants to look further in to MCMC methods other than just metropolis hastings and gibbs sampling. We want to look into Hamiltonian MC, find out how this method compares to others on different models. By this we mean if we choose different models, do some methods converge more quickly then others or do some methods simply fail to successfully approximate the posterior distribution at all.

2. We need to find a way to actually evaluate how well a MCMC method performs. This is one of the things we will research and try to apply.

3. We will compare HMC with metropolis hastings. We can explore extentions to HMC through other implementations like Stan and Edward. If we still have more time remaining we wanted to look into Langevian dynamics which we can implement through Edward (this is in python).

4. The models we want to fit are: 
  a. The conjugate normal-normal model
  b. A multimodal normal mixture model
  c. A logistic regression model
  d. If we have more time we intend to examine how HMC works with other models.
  
5. We can easily simulate data to check how or methods work. We will try to apply these methods to some real data
