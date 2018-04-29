// naive bayes

data {
  // training data
  int<lower = 1> K;               // num classes
  int<lower = 0> N;               // num docs
  int<lower = 1> D;               // image dimension
  vector[D] x[N];                 // image n
  int<lower = 0, upper = 10> y[N]; // class n
  
  // hyperparameters
  vector<lower = 0>[K] alpha;     // class prob prior
  real a;                         // sigma shape
  real b;                         // sigma rate
}

parameters {
  simplex[K] theta;   // topic prevalence
  vector[D] mu[K];    // class k mean
  real<lower = 0> sigsq[K];  // class k variance
}

model {
  // priors
  theta ~ dirichlet(alpha);
  for (k in 1:K) {
    mu[k] ~ multi_normal(rep_vector(0.0 ,D), diag_matrix(rep_vector(1.0, D)));
    sigsq[k] ~ inv_gamma(a / 2, b / 2);
  }
  
  // likelihood, including latent category
  for (n in 1:N) {
    x[n] ~ multi_normal(mu[y[n]], diag_matrix(rep_vector(sigsq[y[n]], D)));
  }
}