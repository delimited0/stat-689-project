# normal inverse gamma

library(truncdist)

# Metropolis Hastings ----
# component wise proposal

nig_likelihood <- function(x) {
  function(params) {
    sum(dnorm(x, params[[1]], params[[2]], log = TRUE))
  }
}

nig_prior <- function() {
  function(params) {
    - log(params[[2]]) 
  }
}

nig_proposal <- function(mu_s, sigma_s) {
  list(
    function(params) rnorm(1, params[[1]], mu_s),
    function(params) rtrunc(1, spec = "norm", a = 0, b = Inf, 
                            mean = params[[2]], sd = sigma_s, log = TRUE)
  )
}

nig_prop_density <- function(mu_s, sigma_s) {
  list(
    function(p, q) dnorm(p[[1]], q[[1]], mu_s, log = TRUE),
    function(p, q) dtrunc(p[[2]], spec = "norm", a = 0, b = Inf, 
                          mean = q[[2]], sd = sigma_s, log = TRUE)
  )
}


# Gibbs sampler ----

nig_cond_mu <- function(x) {
  function(params) {
    rnorm(1, mean(x), params[1] / length(x))
  }
} 

nig_cond_sigma <- function(x) {
  function(params) {
    n <- length(x)
    1 / rgamma(1, shape = (n - 1) / 2, .5 * sum((x - params[1]) ^ 2))
  }
}