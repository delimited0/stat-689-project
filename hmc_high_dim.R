hmc <- function(param_init, U, dU, e, L, iters) {
  # hamiltonian monte carlo
  # 
  # arguments #
  # param_init: initial parameters
  #
  # U: potential energy function
  #   arguments: params
  # 
  # dU: potential energy gradient function
  #  arguments: params
  # 
  # e: leapfrog step size
  # 
  # L: number of leapfrog steps
  # 
  # iters: number of iterations
  
  mom_history <- matrix(data = NA, nrow = iters+1, ncol = length(param_init))
  param_history <- make_history(param_init, iters)
  param_curr <- param_init
  
  for (k in 2:iters) {
    q <- param_curr
    p <- rnorm(length(param_curr), 0, 1)  
    p_curr <- p
    mom_history[k, ] <- p_curr
    
    # leapfrog steps
    p <- hmc_step(p, dU(q), -e / 2)
    for (i in 1:L) {
      q <- hmc_step(q, p, e) 
      if (i != L) p <- hmc_step(p, dU(q), -e)
    } 
    p <- - hmc_step(p, dU(q), -e / 2)
    
    # start end energies
    U_curr <- U(param_curr)
    K_curr <- sum(p_curr ^ 2) / 2
    U_prop <- U(q)
    K_prop <- sum(p^2) / 2
    
    # accept-reject
    if (runif(1) < exp(U_curr + K_curr - U_prop - K_prop))
      param_curr <- q
    
    param_history <- save_sample(param_history, param_curr, k)
  }
  
  return(list(samples = param_history, momentum = mom_history))
}

hmc_step <- function(a, b, e) {
  # move a by e * b
  a + e * b
}