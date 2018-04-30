u = function(q){
  a=0
  for (i in 1:21){
    a[i]=(q[i]-11 + i)^2
  }
  return(sum(a))
}
du = function(q){
  b = 0
  for (i in 1:21){
    b[i] = 2*(q[i] - 11 + i)
  }
  return(b)
}
hmc <- function(param_init, U, dU, e, L, iters) {
  # hamiltonian monte carlo
  #
  # arguments #
  # param_init:
  #
  
  mom_history <- matrix(data = NA, nrow = iters+1, ncol = length(param_init))
  param_history <- matrix(data = NA, nrow = iters+1, ncol = length(param_init))
  param_history[1, ] <- param_init
  param_curr <- param_history[1, ]
  
  for (k in 2:iters) {
    q <- param_curr
    p <- rnorm(length(param_curr), 0, 1) 
    p_curr <- p
    mom_history[k, ] <- p_curr
    
    # leapfrog steps
    p <- hmc_step(p, dU(q), -e / 2)
    for (i in 1:L) {
      q <- hmc_step(q, p, e)
      if (i != L) p <- hmc_step(p, dU(q), -e)
    }
    p <- -hmc_step(p, dU(q), -e / 2) #the negation is unnecessary
    # start end energies
    U_curr <- U(param_curr)
    K_curr <- sum(p_curr ^ 2) / 2
    U_prop <- U(q)
    K_prop <- sum(p^2) / 2
    
    # accept-reject
    if (runif(1) < exp(U_curr + K_curr - U_prop - K_prop)) param_curr <- q
    
    param_history[k, ] <- param_curr
  }
  
  # return(list(samples = param_history, momentum = mom_history))
  return(param_history)
}

hmc_step <- function(a, b, e) {
  # move a by e * b
  a + e * b
}
set.seed(3)
qtracker = hmc(rep(0,21), u, du, 0.18, 10, 10000)
qtracker = qtracker[1:10000,]
length(unique(qtracker))/(21*10000)




xtracker = matrix(NA, nrow = 10000, ncol = 21)
xtracker[1,] = rep(0,21)
for (i in 2:10000){
  xold = xtracker[i-1,]
  xprop = xold + rnorm(length(xold), sd=0.3)
  reject = runif(1)
  if  (-u(xprop) +u(xold) > log(reject)) xold = xprop
  xtracker[i,] = xold
}
length(unique(xtracker))/(21*10000)


par(mfrow=c(4,6))
for(i in 1:21){
  plot(qtracker[,i], ylab = c("q", i))
}

par(mfrow=c(4,6))
for(i in 1:21){
  plot(xtracker[,i], xlab = c("x", i))
}

par(mfrow=c(4,6))
for(i in 1:21){
  acf(xtracker[2500:10000,i], main = c("x",i))
}

par(mfrow=c(4,6))
for(i in 1:21){
  acf(qtracker[2500:10000,i], main = c("q",i))
}


sum(colMeans(qtracker[2500:10000,]) - 10:-10)^2
sum(colMeans(xtracker[2500:10000,]) - 10:-10)^2


effectiveSize(mcmc(xtracker))
effectiveSize(mcmc(qtracker))



microbenchmark(qtracker = hmc(rep(0,21), u, du, 0.2, 10, 10000))

microbenchmark(
  {xtracker[1,] = rep(0,21)
for (i in 2:10000){
  xold = xtracker[i-1,]
  xprop = xold + rnorm(length(xold), sd=0.3)
  reject = runif(1)
  if  (-u(xprop) +u(xold) > log(reject)) xold = xprop
  xtracker[i,] = xold
}
})
