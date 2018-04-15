metropolis_hastings <- function(param_init, likelihood, prior, proposal, 
                                prop_density, iters) {
# metropolis_hastings
#
### arguments ###
# param_init: list of initial parameters
#
# likelihood: log likelihood function, function of parameters only
#   arguments: params, vector of parameters
#   value: log likelihood evaluated at params
#
# prior: log prior function
#   arguments: params, vector of parameters
#   value: prior evaluated at params
#
# proposal: function that samples from proposal distribution, OR list of such 
#           functions
#   arguments: params, vector of parameters
#   value: proposed parameter vector
#
# prop_density: proposal density function, OR list of such functions
#   arguments: curr_param, vector of current parameters
#              prop_param, vector of proposed parameters
#   value: proposal density evaluated at prop_param, given curr_param
#
# iters: number of iterations
#  
### return value ###
#   samples - iter x D matrix of samples
#   accept_rate - acceptance rate
  
  param_history <- make_history(param_init, iters)
  param_curr <- param_init
  accepts <- 0
  
  for (k in 2:iters) {
    
    if (is.function(proposal)) {
      # block sampling
      propose_param <- proposal(param_curr)
      a <- exp(likelihood(propose_param) + prior(propose_param) + 
                 prop_density(propose_param, param_curr) -
               likelihood(param_curr) - prior(param_curr) -
                 prop_density(param_curr, propose_param))
      a <- min(1, a) 
      u <- runif(1)
      
      if (is.na(a) | is.infinite(a) | is.nan(a))
        break
      
      if (u < a) param_curr <- propose_param
      
      accepts <- accepts + (u < a)
    }
    else {
      # component sampling
      
      for (j in 1:length(proposal)) {
        comp_proposal <- proposal[[j]]
        comp_density <- prop_density[[j]]
        
        propose_param <- param_curr
        propose_param[[j]] <- comp_proposal(propose_param)
        
        a <- exp(likelihood(propose_param) + prior(propose_param) +
                  comp_density(propose_param, param_curr) -
                 likelihood(param_curr) - prior(param_curr) -
                  comp_density(param_curr, propose_param))
        a <- min(1, a) 
        u <- runif(1)
        
        if (u < a) param_curr <- propose_param
        
        accepts <- accepts + (u < a)
      }
    }
    
    param_history <- save_sample(param_history, param_curr, k)
  }

  list(samples = param_history, accept_rate = accepts / iters)
}