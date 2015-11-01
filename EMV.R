

f_beta=function(b,ech,nbdef,censure)
{
  1/b + (sum(log(ech[ech<censure])))/nbdef - (sum((ech^b)*log(ech)))/(sum(ech^b))
}

EMV=function(t,c)
{
  m=length(t[t<c])
  beta.est = uniroot(f_beta,lower=0.001,upper=20,ech=t,nbdef=m,censure=c)$root
  eta.est = (1/m * sum(t^beta.est))^(1/beta.est)   #si m=0 on a des problèmes!! et cela peut arriver...
  return(list(beta_EMV=beta.est,eta_EMV=eta.est))
}

#figure 1 - convergence
betaVrai=0.5
etaVrai=100
c=40
n=seq(from=50,to=1500,by=10)
EMV.est=c()

for (i in n)
{
  repeat
  {
    t <- rweibull(i,betaVrai,etaVrai)
    t[t>c]=c
    if (length(t[t<c])!=0) {break}  #pour ne pas avoir m=0
  }
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
}

par(mfrow = c(1,2))
plot(n,EMV.est[,1],ylab=expression(beta* " estime"))
abline(h=betaVrai,col=2)
plot(n,EMV.est[,2],ylab=expression(eta*" estime"))
abline(h=etaVrai,col=2)


#figure 2 - normalité asymptotique
betaVrai=1.5
etaVrai=100
n=500
c=40
k=1000  #on fait k realisations 
EMV.est=matrix(nrow=k,ncol=2) 
for (i in 1:k)
{
  repeat
  {
    t <- rweibull(n,betaVrai,etaVrai)
    t[t>c]=c
    if (length(t[t<c])!=0) {break}  
  }
  EMV.est[i,1]=EMV(t,c)$beta_EMV
  EMV.est[i,2]=EMV(t,c)$eta_EMV
}
par(mfrow = c(1,2))
hist(EMV.est[,1],breaks=20,
   xlab="",main=expression("histogramme pour "*hat(beta)))
hist(EMV.est[,2],breaks=20,
   xlab="",main=expression("histogramme pour "*hat(eta)))

