# normal inverse gamma

library(truncdist)

# Metropolis Hastings ----
# component wise proposal

nig_likelihood <- function(x) {
  n <- length(x)
  function(params) {
    v <- params[[2]]
    mu <- params[[1]]
    
    - ((n + 2) / 2) * log(v) - sum((x - mu) ^ 2) / (2 * v)
  }
}

nig_prior <- function() {
  function(params) {
    0#* log(params[[2]])
  }
}

nig_proposal <- function(mu_s, sigma_s) {
  list(
    function(params) rnorm(1, params[[1]], mu_s),
    function(params) {
      r <- rnorm(1, params[[2]], sigma_s)
      r * sign(r)
    }
    # rtrunc(1, spec = "norm", a = 0, b = Inf, 
    #                       mean = params[[2]], sd = sigma_s, log = TRUE)
  )
}

nig_prop_density <- function(mu_s, sigma_s) {
  list(
    function(p, q) dnorm(p[[1]], q[[1]], mu_s, log = TRUE),
    function(p, q) dnorm(p[[2]], q[[2]], sigma_s, log = TRUE)
      # dtrunc(p[[2]], spec = "norm", a = 0, b = Inf, 
      # mean = q[[2]], sd = sigma_s, log = TRUE)
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

# HMC ----

nig_U <- function(x) {
  n <- length(x)
  function(params) {
    v <- params[2]
    mu <- params[1]
    ((n + 2) / 2) * log(v) + sum((x - mu) ^ 2) / (2 * v)
  }
}

nig_dU <- function(x) {
  function(params) {
    n <- length(x)
    v <- params[2]
    mu <- params[1]
    c(
      - sum(x - mu) / v,
      (n + 2) / (2 * v) - sum((x - mu) ^ 2) / (2 * v ^ 2)
    )
  }
}
