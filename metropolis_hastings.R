metropolis_hastings <- function(init_param, proposal, iters, 
                                likelihood, prior, prop_density) {
  # init_param: initial parameters
  # proposal: function that samples from proposal distribution
  # iters: number of iterations
  # likelihood: log likelihood function, function of parameters only
  # prior: log prior function
  # prop_density: proposal distribution function
  # return iter x D matrix of samples
  
  param_history <- matrix(data = NA, nrow = iters+1, ncol = length(init_param))
  param_history[1, ] <- init_param
  param_curr <- param_history[1, ]
  
  for (k in 2:iters) {
    propose_param <- proposal(param_curr)
    
    a <- exp(likelihood(propose_param) + prior(propose_param) + 
               prop_density(propose_param, param_curr) -
             likelihood(param_curr) - prior(param_curr) -
               prop_density(param_curr, propose_param))
    a <- min(1, a) 
    u <- runif(1)
    
    param_curr <- (u < a) * propose_param + (u >= a) * param_curr
    
    param_history[k, ] <- param_curr
  }

  param_history
}