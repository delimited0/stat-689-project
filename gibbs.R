gibbs <- function(param_init, iters, ...) {
  #' gibbs
  #' gibbs sampling
  #' 
  #' Arguments
  #' param_init: initial parameters
  #' iters: number of samples
  #' ...: functions for sampling from full conditionals
  #'
  #' Return Value
  #' param_history: iters x D posterior sample matrix

  conditionals <- list(...)
  
  param_history <- matrix(data = NA, nrow = iters+1, ncol = length(param_init))
  param_history[1, ] <- param_init
  param_curr <- param_history[1, ]
  
  for (k in 2:iters) {
    for (j in 1:length(conditionals)) {
      cond <- conditionals[[j]]
      param_curr[j] <- cond(param_curr[-j])
    }
    param_history[k, ] <- param_curr
  }
  
  param_history
}