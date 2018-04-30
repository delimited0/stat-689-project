# utilities

make_history <- function(param_init, iters) {
# make sample history list and intialize
# arguments
# param_init: list of initial parameters
# iters: number of sampling iterations

  lapply(param_init, function(p) {
    if (is.vector(p) || nrow(p) == 1) {
      pmat <- matrix(NA, nrow = iters, ncol = length(p))
      pmat[1, ] <- p
    }
    else if (is.matrix(p)) {
      pmat <- array(NA, dim = c(iters, nrow(p), ncol(p)))
      pmat[1, , ] <- p
    }
    return(pmat)
  })
}

save_sample <- function(param_history, param, n) {
# add sample to history
# arguments
# param_history: list of param histories
# param: list of new parameters
# n: iteration number
  
  mapply(function(ph, p) {
    if (is.vector(p) || nrow(p) == 1) 
      ph[n, ] <- p
    else if (is.matrix(p))
      ph[n, , ] <- p
    ph
  }, param_history, param, SIMPLIFY = FALSE)
}

make_mom_history <- function(param_init, iters) {
  lapply(param_init, function(p) {
    if (is.vector(p)) {
      pmat <- matrix(NA, nrow = iters, ncol = length(p))
    }
    else if (is.matrix(p)) {
      pmat <- array(NA, dim = c(iters, nrow(p), ncol(p)))
    }
    return(pmat)
  })
}