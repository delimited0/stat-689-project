# normal normal functions

nn_likelihood <- function(x) {
  function(params) {
    sum(dnorm(x, params[1], 1, log = TRUE))
  }
}

nn_prior <- function() {
  function(params) 0
}

nn_proposal <- function(sigma_mean) {
  function(params) {
    rnorm(1, mean = params[1], sd = sigma_mean)
  }
}

nn_prop_density <- function(sigma_mean, lambda_rate) {
  function(p, q) 0
}

