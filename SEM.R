#Methode SEM

#C on simule les donn?es qu'on a pas pu observ?es (qu'on avait mises ? c)
beta=1.5
eta=100
c=40
n=50
t <- rweibull(n,beta,eta)
t[t>c]=c
m=length(t[t<c])


simu=function(c,eta,beta,j)
{
  u=runif(j)
  return(eta*((-log(u)+(c/eta)^beta)^(1/beta)))
}


f2=function(b,t,n)
{
  1/b + (sum(log(t))/n) - ((sum((t^b)*log(t)))/(sum(t^b)))
}

EMV2=function(t)
{
  b.est = uniroot(f2,lower=0.001,upper=10,t=t,n=n)$root
  eta.est = (1/n * sum(t^b.est))^(1/b.est)
  return(list(beta_EMV=b.est,eta_EMV=eta.est))
}


eta0=50
beta0=3
K=1000
res=matrix(nrow = K+1,ncol = 2)
res[1,1]=beta0
res[1,2]=eta0
for (k in 1:K)
{
  t[t>=c]=simu(c,eta0,beta0,n-m)
  beta0=EMV2(t)$beta_EMV
  eta0=EMV2(t)$eta_EMV
  res[k+1,1]=beta0
  res[k+1,2]=eta0
}

x=seq(from=1,to=K+1,by=1)
par(mfrow = c(1,2))
plot(x,res[,1],xlab="nb it?rations",ylab="beta estim?")
abline(h=beta,col="red")
plot(x,res[,2],xlab="nb it?rations",ylab="eta estim?")
abline(h=eta,col="red")



