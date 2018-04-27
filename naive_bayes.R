# naive bayes

library(mvnfast)
library(MCMCpack)

# Metropolis hastings ----

# Gibbs ----

nb_cond_pi <- function(y, alpha) {
  ns <- table(y) + alpha
  function(params) {
    rdirichlet(1, ns)  
  }
}

nb_cond_mu <- function(y, X) {
  ns <- table(y)
  d <- ncol(X)
  X_k_sum <- lapply(sort(unique(y)) + 1, function(k) colSums(X[y == k, ]))
  function(params) {
    sigmas <- params[["sigmas"]]
    mapply(FUN = function(n, sig, X_sum) {
      Q <- diag((1 + n / sig) ^ (-1), nrow = d, ncol = d)
      b <- X_sum / sig
      rmvn(1, Q %*% b, Q)  
    }, ns, sigmas, X_k_sum)
  }
}

nb_cond_sigmas <- function(y, X, a, b) {
  ns <- table(y)
  d <- ncol(X)
  X_k_sum <- lapply(sort(unique(y)) + 1, function(k) colSums(X[y == k, ]))
  X_k_norm <- lapply(sort(unique(y)) + 1, function(k) sum(X[y == k, ] ^ 2))
  K <- length(unique(y))
  function(params) {
    mus <- params[["mus"]]
    mapply(function(n, k, X_sum, X_norm) {
      shape <- (n * d + a) / 2
      rate <- 
        .5 * (X_norm - 2 * mus[, k] %*% X_sum + n * mus[, k] %*% mus[, k] + b)
      rinvgamma(1, shape, rate)
    }, ns, 1:K, X_k_sum, X_k_norm)
  }
}

# HMC ----
