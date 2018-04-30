# naive bayes

library(mvnfast)
library(MCMCpack)
library(greta)

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

# hmc ----

nb_hmc_model <- function(y, X, alpha, a, b) {
  K <- length(unique(y))
  d <- ncol(X)
  N <- nrow(X)
  ns <- table(y)
  X_k <- lapply(1:K, function(k) as_data(X[y == (k - 1), ]))
  # y <- as_data(y)
  
  pi = dirichlet(alpha)
  sigmas = inverse_gamma(a/2, b/2, dim = K)
  mu = multivariate_normal(rep(0, d), diag(1, nrow = d), dim = K)
  
  identity <- diag(1, d)
  
  for (k in 1:K) {
    distribution(X_k[[k]]) = 
      multivariate_normal(t(mu[k,]), sigmas[k] * identity, dim = nrow(X_k[[k]]))
  }
  # y = categorical(t(pi), dim = N)
  
  model(pi, sigmas, mu)
}

# prediction ----

predict_nb <- function(mus, sigmas, X_new) {
# mus: matrix of class means
# sigmas: vector class variances
# X_new: matrix of observations
  K <- length(sigmas)
  d <- ncol(X_new)
  probs <- matrix(NA, nrow = nrow(X_new), ncol = K)
  
  for (k in 1:K)
    probs[, k] <- dmvn(X_new, mus[, k], diag(sigmas[k], nrow = d), log = TRUE) 
  
  apply(probs, 1, which.max)
}
