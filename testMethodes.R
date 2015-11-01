#comparaison des deux methodes
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
                     upper=10,ech=x,N=n)$root
  eta.est = (1/n * sum(x^beta.est))^(1/beta.est)
  return(list(beta_EMV_SEM=beta.est,eta_EMV_SEM=eta.est))
}

SEM=function(t,c,beta0=1,eta0=50,epsilon=0.001,iterMax=40)
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

betaVrai=0.5
etaVrai=100
c=40
n=25
k=50 #nb d'echantillons simules
M=c()
EMV.est=c()
SEM.est=c()
Res=matrix(ncol=5,nrow=2)
total_c=0
for (i in 1:k)
{
  t=rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  M=c(M,length(t[t<c]))
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
  if (length(t[t<c])==0) {total_c=total_c+1; next;}  
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
}
  
Res[1,1]=mean(EMV.est[,1])
Res[2,1]=mean(EMV.est[,2])
Res[1,2]=sd(EMV.est[,1])
Res[2,2]=sd(EMV.est[,2])
Res[1,3]=mean(SEM.est[,1])
Res[2,3]=mean(SEM.est[,2])
Res[1,4]=sd(SEM.est[,1])
Res[2,4]=sd(SEM.est[,2])
Res[1,5]=mean(M)


betaVrai=1.2
M=c()
EMV.est=c()
SEM.est=c()
Res2=matrix(ncol=5,nrow=2)
pb2=0
for (i in 1:k)
{
  t=rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  M=c(M,length(t[t<c]))
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
  if (length(t[t<c])==0) {pb2=pb2+1; next;}  
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
  
}

Res2[1,1]=mean(EMV.est[,1])
Res2[2,1]=mean(EMV.est[,2])
Res2[1,2]=sd(EMV.est[,1])
Res2[2,2]=sd(EMV.est[,2])
Res2[1,3]=mean(SEM.est[,1])
Res2[2,3]=mean(SEM.est[,2])
Res2[1,4]=sd(SEM.est[,1])
Res2[2,4]=sd(SEM.est[,2])
Res2[1,5]=mean(M)

EMV=function(t,c)
{
  m=length(t[t<c])
  beta.est = uniroot(f_beta,lower=0.001,upper=20,ech=t,nbdef=m,censure=c)$root
  eta.est = (1/m * sum(t^beta.est))^(1/beta.est)   #si m=0 on a des problèmes!! et cela peut arriver...
  return(list(beta_EMV=beta.est,eta_EMV=eta.est))
}

EMV_SEM=function(x)
{
  n=length(x)
  beta.est = uniroot(f_beta_SEM,lower=0.001,
                     upper=20,ech=x,N=n)$root
  eta.est = (1/n * sum(x^beta.est))^(1/beta.est)
  return(list(beta_EMV_SEM=beta.est,eta_EMV_SEM=eta.est))
}

betaVrai=2
M=c()
EMV.est=c()
SEM.est=c()
Res3=matrix(ncol=5,nrow=2)
pb3=0
for (i in 1:k)
{
  t=rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  M=c(M,length(t[t<c]))
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
  if (length(t[t<c])==0) {pb3=pb3+1; next;}  
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
  
}

Res3[1,1]=mean(EMV.est[,1])
Res3[2,1]=mean(EMV.est[,2])
Res3[1,2]=sd(EMV.est[,1])
Res3[2,2]=sd(EMV.est[,2])
Res3[1,3]=mean(SEM.est[,1])
Res3[2,3]=mean(SEM.est[,2])
Res3[1,4]=sd(SEM.est[,1])
Res3[2,4]=sd(SEM.est[,2])
Res3[1,5]=mean(M)

EMV=function(t,c)
{
  m=length(t[t<c])
  beta.est = uniroot(f_beta,lower=0.001,upper=70,ech=t,nbdef=m,censure=c)$root
  eta.est = (1/m * sum(t^beta.est))^(1/beta.est)   #si m=0 on a des problèmes!! et cela peut arriver...
  return(list(beta_EMV=beta.est,eta_EMV=eta.est))
}

EMV_SEM=function(x)
{
  n=length(x)
  beta.est = uniroot(f_beta_SEM,lower=0.001,
                     upper=70,ech=x,N=n)$root
  eta.est = (1/n * sum(x^beta.est))^(1/beta.est)
  return(list(beta_EMV_SEM=beta.est,eta_EMV_SEM=eta.est))
}

betaVrai=3
M=c()
EMV.est=c()
SEM.est=c()
Res4=matrix(ncol=5,nrow=2)
pb4=0
for (i in 1:k)
{
  t=rweibull(n,betaVrai,etaVrai)
  t[t>c]=c
  M=c(M,length(t[t<c]))
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
  if (length(t[t<c])==0) {pb4=pb4+1; next}  
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
  
}

Res4[1,1]=mean(EMV.est[,1])
Res4[2,1]=mean(EMV.est[,2])
Res4[1,2]=sd(EMV.est[,1])
Res4[2,2]=sd(EMV.est[,2])
Res4[1,3]=mean(SEM.est[,1])
Res4[2,3]=mean(SEM.est[,2])
Res4[1,4]=sd(SEM.est[,1])
Res4[2,4]=sd(SEM.est[,2])
Res4[1,5]=mean(M)
