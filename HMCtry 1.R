require(rootSolve)

dens<-function(q)
{
 exp(-0.5*(q+10)^2)+exp(-0.5*(q-10)^2) + exp(-0.5*(q+20)^2)+exp(-0.5*(q-20)^2)
}

U = function(q){
  return(-log(dens(q)))
  }



gradU<-function(q,U)
{
    h<-.000000001
    for(i in 1)
    gradg<-(U(q+h)-U(q-h))/(2*h)
}

HMC=function(U2,eps,L,currq,gradU2)
{
  q=currq
  #p=rnorm(length(q),0,1)
  p = rt(length(q), df = 2)
  currp=p
  
  p=p-eps*gradU2(q)
    #gradient(U2,q,centered = TRUE,pert=1*10^-10)/2
  if (is.nan(p)==T){
    p = 0
  }
  for(i in 1:L)
  {
    q=q+eps*p
   # print(q)
    if(i!=L)
    {
      p=p-eps*gradient(U2,q,centered = TRUE,pert = 1*10^-10)
   #   print(p)
      if(is.nan(p) == T){
        p = 0
      }
    }
    
  }
  p=p-eps*gradient(U2,q,centered = TRUE,pert=1*10^-10)/2
  if(is.nan(p) == T){
    p = 0
  }
  #p=-p # this does nothing but satisfies some theory
  
  currU=U2(currq)
  a <<- q
  currK=sum(currp^2)/2
  propU=U2(q)
  propK=sum(p^2)/2
  #print(c(currU, propU, currK, propK))
  if(runif(1)<exp(currU-propU+currK-propK))
  {
    return(q)
  }
  else
  {
    return(currq)
  }

}
#qtracker<-c()
#qtracker[1]<-1
qtracker = 1
for(i in 2:100000)
{
  qtracker[i]<-HMC(U,2,5,qtracker[i-1])
}
xtracker = 1
for (i in 2:100000)
{
xtrackernew = xtracker[i-1] + rnorm(1, sd = 7)
if ((dens(xtrackernew) / dens(xtracker[i-1]))> runif(1))
    {
      xtracker[i] = xtrackernew
}
else xtracker[i] = xtracker[i-1]
}

par(mfrow=c(2,1))
plot(qtracker[10000:100000])
plot(xtracker[10000:100000])
