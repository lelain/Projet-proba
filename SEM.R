#Methode SEM

#fonction permettant de simuler N realisations de loi de Weibull de parametre eta0 et beta0, conditionnee par X>c
simu=function(c,eta0,beta0,N)
{
  u=runif(N)
  return(eta0*((-log(u)+(c/eta0)^beta0)^(1/beta0)))
}

#la log-vraisemblance pour l'?chantillon de d?part

logVrai=function(ech,c,beta,eta)
{
  m=length(ech[ech<c])
  m*log(beta) - m*log(eta) + (beta-1)*(sum(ech[ech<c])) 
           -m*(beta-1)*log(eta) - (sum(ech^beta))/(eta^beta)                               
}


f_beta_SEM=function(b,ech,N)
{
  1/b + (sum(log(ech))/N) - ((sum((ech^b)*log(ech)))/(sum(ech^b)))
}

EMV_SEM=function(x)
{
  n=length(x)
  beta.est = uniroot(f_beta_SEM,lower=0.001,
                         upper=20,ech=x,N=n)$root
  eta.est = (1/n * sum(x^beta.est))^(1/beta.est)
  return(list(beta_EMV_SEM=beta.est,eta_EMV_SEM=eta.est))
}

SEM=function(t,c,beta0=1,eta0=50,epsilon=0.001,iterMax=30)
{
  res=c(beta0,eta0) 
  n=length(t)
  m=length(t[t<c])  
  V=logVrai(t,c,beta0,eta0)   
  beta.est=beta0
  eta.est=eta0
  
  i=1   
  repeat
  {
    tSimu=t
    tSimu[tSimu>=c]=simu(c,eta.est,beta.est,n-m)
    beta.est=EMV_SEM(tSimu)$beta_EMV_SEM
    eta.est=EMV_SEM(tSimu)$eta_EMV_SEM
    res=rbind(res,c(beta.est,eta.est))
    V=c(V,logVrai(t,c,beta.est,eta.est))
    if (abs(V[i+1]-V[i])<epsilon) {break}
    if (i==iterMax) {break}
    i=i+1
  }
  return(list(beta_SEM=beta.est,eta_SEM=eta.est,iter=i,logVrai=V))
}


#test convergence
betaVrai=0.5
etaVrai=100
c=40
n=seq(from=50,to=1500,by=10)
SEM.est=c()

for (i in n)
{
  t <- rweibull(i,betaVrai,etaVrai)
  t[t>c]=c
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
}

par(mfrow = c(1,2))
plot(n,SEM.est[,1],ylab=expression(beta* " estime"))
abline(h=betaVrai,col=2)
plot(n,SEM.est[,2],ylab=expression(eta*" estime"))
abline(h=etaVrai,col=2)

  
#figure 2 - loi des estimateurs
betaVrai=1.5
etaVrai=100
n=500
c=40
k=1000   
SEM.est=matrix(nrow=k,ncol=2) 
for (i in 1:k)
{
  t <- rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  SEM.est[i,1]=SEM(t,c)$beta_SEM
  SEM.est[i,2]=SEM(t,c)$eta_SEM
}
par(mfrow = c(1,2))
hist(SEM.est[,1],breaks=20,
     xlab="",main=expression("histogramme pour "*hat(beta)))
hist(SEM.est[,2],breaks=20,
     xlab="",main=expression("histogramme pour "*hat(eta)))


