library(mvnfast)

nb_metropolis_hastings <- function(param_init, iters, X, y, alpha, a, b, 
                                   mu_s = 1 / ncol(X), sigma_s = 1, pi_s = 1) {
  # metropolis hastings for Gaussian Naive Bayes 
  #
  # param_init: list with arguments
  #     1. pi: vector of class probabilities
  #     2. mus: matrix of class means
  #     3. sigmas: vector of class variances                                  
  #
  # iters: number of samples
  # X:  
  # y:
  #
  # mu_s should be VERY SMALL
  
  ns <- table(y)
  K <- length(unique(y))
  d <- ncol(X)
  
  pi_posterior <- function(pi, ns, alpha)
    sum((ns + alpha - 1) * log(pi))
  
  mu_posterior <- function(mu, sigma, X_sum, n) {
    .5 * ((2 * mu %*% X_sum - n * mu %*% mu) / sigma - mu %*% mu)
  }
  
  sigma_posterior <- function(sigma, mu, X_sum, n) {
    - (n * d / 2) * log(sigma) - (a + 1) * log(sigma) - (b / sigma) -
      .5 * (2 * mu %*% X_sum - n * mu %*% mu) / sigma
  }

  param_history <- make_history(param_init, iters)
  param_curr <- param_init
  accept_pi <- 0; accept_mu <- 0; accept_sigma <- 0
  X_sums <- lapply(sort(unique(y)) + 1, function(k) colSums(X[y == k, ]))
  
  for (k in 2:iters) {
    
    # sample pi -- 
    pi_curr <- param_curr[["pi"]]        
    pi_prop <- rdirichlet(1, rep(pi_s, 10))[1, ]
    pi_ratio <- 
      pi_posterior(pi_prop, ns, alpha) + log(ddirichlet(pi_prop, pi_curr)) -
      pi_posterior(pi_curr, ns, alpha) - log(ddirichlet(pi_curr, pi_prop))
    a <- min(exp(pi_ratio), 1) 
    if (a > runif(1)) {
      param_curr[["pi"]] <- pi_prop                                    
      accept_pi <- accept_pi + 1
    }
    
    # sample mu --
    sigmas <- param_curr[["sigmas"]]
    for (j in 1:K) {    
      mu_curr <- param_curr[["mus"]]
      mu_curr_j <- mu_curr[, j]
      mu_prop_j <- rmvn(n = 1, mu = mu_curr_j, sigma = mu_s * diag(d))[1, ]
      mu_ratio <- 
        mu_posterior(mu_prop_j, sigmas[j], X_sums[[j]], ns[j]) - 
        mu_posterior(mu_curr_j, sigmas[j], X_sums[[j]], ns[j])
      a <- min(exp(mu_ratio), 1) 
      if (a > runif(1)) {
        mu_curr[, j] <- mu_prop_j
        param_curr[["mus"]] <- mu_curr
        accept_mu <- accept_mu + 1
      }
    }
    
    # sample sigma --  
    mus <- param_curr[["mus"]]
    for (j in 1:K) {
      sigma_curr <- param_curr[["sigmas"]]
      sigma_curr_j <- sigma_curr[j]
      sigma_prop_j <- rnorm(n = 1, mean = sigma_curr_j, sd = sigma_s)
      if (sigma_prop_j < 0)
        sigma_prop_j <- - sigma_prop_j        
      sigma_ratio <- 
        sigma_posterior(sigma_prop_j, mus[, j], X_sums[[j]], ns[j]) -
        sigma_posterior(sigma_curr_j, mus[, j], X_sums[[j]], ns[j])
      a <- min(exp(sigma_ratio), 1) 
      if (a > runif(1)) {
        sigma_curr[j] <- sigma_prop_j    
        param_curr[["sigmas"]] <- sigma_curr
        accept_sigma <- accept_sigma + 1
      }
    }
    
    param_history <- save_sample(param_history, param_curr, k)
  }
  
  list(samples = param_history, 
       accept_rate = c("accept_pi" = accept_pi / K, 
                       "accept_mu" = accept_mu / K,
                       "accept_sigma" = accept_sigma / K
       ) / iters
  )
}