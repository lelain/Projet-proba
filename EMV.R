

f_beta=function(b,ech,nbdef,censure)
{
  1/b + (sum(log(ech[ech<censure])))/nbdef - (sum((ech^b)*log(ech)))/(sum(ech^b))
}

EMV=function(t,c)
{
  m=length(t[t<c])
  beta.est = uniroot(f_beta,lower=0.001,upper=10,ech=t,nbdef=m,censure=c)$root
  eta.est = (1/m * sum(t^beta.est))^(1/beta.est)   #si m=0 on a des problèmes!! et cela peut arriver...
  return(list(beta_EMV=beta.est,eta_EMV=eta.est))
}

#figure 1 - convergence
betaVrai=1.5
etaVrai=100
c=40
N=300
n=rep(0,N)
EMV.est=matrix(nrow=N,ncol=2)

for (i in 1:N)
{
  n[i]=20+5*i
  t <- rweibull(n[i],betaVrai,etaVrai)
  t[t>c]=c
  
  EMV.est[i,1]=EMV(t,c)$beta_EMV
  EMV.est[i,2]=EMV(t,c)$eta_EMV
}

par(mfrow = c(1,2))
plot(n,EMV.est[,1],ylab=expression(beta* " estimé"))
abline(h=betaVrai,col=2)
plot(n,EMV.est[,2],ylab=expression(eta*" estimé"))
abline(h=etaVrai,col=2)


#figure 2 - normalité asymptotique
betaVrai=1.5
etaVrai=100
n=500
c=40
k=1000  #on fait k réalisations 
EMV.est=matrix(nrow=k,ncol=2) 
for (i in 1:k)
{
  t <- rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  EMV.est[i,1]=EMV(t,c)$beta_EMV
  EMV.est[i,2]=EMV(t,c)$eta_EMV
}
par(mfrow = c(1,2))
hist(EMV.est[,1],breaks=20,xlab="",main=expression("histogramme pour "*hat(beta)))
hist(EMV.est[,2],breaks=20,xlab="",main=expression("histogramme pour "*hat(eta)))

#3. test , évaluation de la moyenne et variance
#estimation de la variance et du biais
betaVrai=2
etaVrai=100
n=25
c=40
k=50  #on fait k réalisations 
EMV.est=matrix(nrow=k,ncol=2) 
for (i in 1:k)
{
  t <- rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  EMV.est[i,1]=EMV(t,c)$beta_EMV
  EMV.est[i,2]=EMV(t,c)$eta_EMV
}

var(EMV.est[,1])
var(EMV.est[,2])

mean(EMV.est[,1])
mean(EMV.est[,2])



n=50
c=40
k=1000  #on fait k r?alisations 
EMV.est=matrix(nrow=k,ncol=2) 
for (i in 1:k)
{
  t <- rweibull(n,beta,eta)
  t[t>c]=c
  EMV.est[i,1]=EMV(t,c)$beta_EMV
  EMV.est[i,2]=EMV(t,c)$eta_EMV
}

var(EMV.est[,1])
var(EMV.est[,2])

mean(EMV.est[,1])-beta
mean(EMV.est[,2])-eta

hist(EMV.est[,1])
hist(EMV.est[,2])
plot(density(EMV.est[,1]))
abline(v=beta,col=2)
plot(density(EMV.est[,2]))
abline(v=eta,col=2)