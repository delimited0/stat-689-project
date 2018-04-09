# normal inverse gamma

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