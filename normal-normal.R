# normal normal functions

nn_likelihood <- function(x) {
  function(params) {
    sum(dnorm(x, params[1], params[2], log = TRUE))
  }
}

nn_prior <- function() {
  function(params) 0
}

nn_proposal <- function(sigma_mean, sigma_sd) {
  function(params) {
    mu <- rnorm(1, mean = params[1], sd = sigma_mean)
    sigma <- rnorm(1, mean = params[2], sd = sigma_sd)
    c(mu, sigma * sign(sigma))
  }
}

nn_prop_density <- function(sigma_mean, lambda_rate) {
  function(p, q) 0
}

