set.seed(3)
qtracker = hmc(1, u, du, 0.005, 10, 100000)
qtracker = qtracker[1:100000]

xtracker = 1
for (i in 2:100000){
  xold = xtracker[i-1]
  xprop = xold + rnorm(1, sd = 0.7)
  reject = runif(1)
  if (exp(-u(xprop))/exp(-u(xold)) > reject) xold = xprop
  xtracker[i] = xold
}
par(mfrow =c(2,2))
plot(qtracker)
plot(xtracker)
hist(qtracker, xlim=c(1.9,2.1), breaks = 10000)
hist(xtracker, xlim=c(1.9,2.1), breaks = 10000)
