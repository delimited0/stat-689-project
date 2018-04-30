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
  
  param_history <- make_history(param_init, iters)
  param_curr <- param_init
  
  for (k in 2:iters) {
    for (j in 1:length(conditionals)) {
      cond <- conditionals[[j]]
      param_curr[[j]] <- cond(param_curr[-j])
    }
    param_history <- save_sample(param_history, param_curr, k)
  }
  
  list(samples = param_history)
}