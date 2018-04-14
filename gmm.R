#Gaussian Mixture functions

#Metropolis Hastings     ----

library(mvtnorm)
library(MCMCpack)

gmm_likelihood <- function(x, k)
{
  function(params)
  {
    mus <- params[1:k]
    sigmas <- params[(k+1):(2*k)]
    piprop <- params[[2*k+1]]
    sum(log(
      mapply(function(mu, covar) {
        dmvnorm(x, mu, covar)
      }, 
      mus, sigmas) %*% piprop
    ))
  }
}

gmm_prior <- function(k, alpha, v, S, mu0)
{  
  function(params)
  {
    mus <- params[1:k]
    sigmas <- params[(k+1):(2*k)]
    piprop <- params[[2*k+1]]
    sum(mapply(function(mu,covar) {
      log(diwish(covar, v, S)) + 
        dmvnorm(x = mu, mean = mu0, sigma = covar, log = TRUE) 
    }, mus, sigmas)) + log(ddirichlet(piprop, alpha))
    
  
  }
}

gmm_proposal<-function(k, mu_sigma, v)
{
    function(params)
    {
      mus <- params[1:k]
      sigmas <- params[(k+1):(2*k)]
      piprop <- params[[2*k+1]]
      mus <- lapply(mus, function(mu) rmvnorm(1, mu, mu_sigma))
      sigmas <- lapply(sigmas, function(sigma) 
        rWishart(1, v, (v - nrow(sigma) - 1) * sigma))
      piprop <- rdirichlet(1, piprop)
      list(mus, sigmas, piprop)
    }
}


gmm_prop_density <- function(k, mu_sigma, v)
{
  function(params_p, params_q)
  {
    mus_p <- params_p[1:k]
    sigmas_p <- params_p[(k+1):(2*k)]
    piprop_p <- params_p[[2*k+1]]
    mus_q <- params_q[1:k]
    sigmas_q <- params_q[(k+1):(2*k)]
    piprop_q <- params_q[[2*k+1]]
    
    sigmas_density <- mapply(function(sigma_p, sigma_q) {
      log(dwish(sigma_p, v, (v - nrow(sigma_q) - 1) * sigma_q)) 
      }, sigmas_p, sigmas_q)
    
    piprop_density <- log(ddirichlet(piprop_p, piprop_q))
    sum(sigmas_density) + piprop_density
  }
}



## unit tests ----
# 
# u <- runif(10) < 1/2
# 
# draws <- u * rmvnorm(10, mean = c(1,2), sigma = diag(2)) + 
#   (1-u) * rmvnorm(10, mean = c(8, 10), sigma = 2 * diag(2))
# 
# 
# L1 <- gmm_likelihood(draws, 2)
# L1(list(c(0, 0), c(0, 0), 9 * diag(2), 5 * diag(2), c(1/3, 2/3)))
# L1(list(c(1, 2), c(8, 10), diag(2), 2 * diag(2), c(1/2, 1/2)))


#P1 <- gmm_prior(2, c(1/2, 1/2), 2, diag(2), c(0, 0))
#P1(list(c(0, 0), c(0, 0), 9 * diag(2), 5 * diag(2), c(1/3, 2/3)))
#P1(list(c(0, 0), c(0, 0), diag(2), diag(2), c(1/2, 1/2)))


#Prop1 <- gmm_props(2, .5*diag(2), 4)
#propdraw <- Prop1(list(c(1, 2), c(8, 10), diag(2), 2 * diag(2), c(1/2, 1/2)))


#propdens1 <- gmm_prop_density(2, .5*diag(2), 4)

# propdens1(list(c(0, 0), c(0, 0), 9 * diag(2), 5 * diag(2), c(1/3, 2/3)),
#          list(c(0, 0), c(0, 0), diag(2), diag(2), c(1/2, 1/2)))
