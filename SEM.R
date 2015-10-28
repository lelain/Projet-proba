#Methode SEM

#fonction permettant de simuler N realisations de loi de Weibull de parametre eta0 et beta0, conditionnee par X>c
simu=function(c,eta0,beta0,N)
{
  u=runif(N)
  return(eta0*((-log(u)+(c/eta0)^beta0)^(1/beta0)))
}

#la log-vraisemblance pour l'?chantillon de d?part
logVrai=function(c,eta,beta,ech)
{
  m=length(tVrai[tVrai<c])
  m*log(beta) - m*log(eta) + (beta-1)*(sum(tVrai[tVrai<c])) 
           -m*(beta-1)*log(eta) - (sum(tVrai^beta))/(eta^beta)                               
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

SEM=function(t,c,beta0=3,eta0=50,epsilon=0.2,iterMax=100)
{
  res=c(beta0,eta0) 
  n=length(t)
  m=length(t[t<c])  
  V=logVrai(c,eta0,beta0,t)   
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
    V=c(V,logVrai(c,eta.est,beta.est,t))
    if (abs(V[i+1]-V[i])<epsilon) {break}
    if (i==100) {break}
    i=i+1
  }
  return(list(beta_SEM=beta.est,eta_SEM=eta.est,iter=i))
}


#Les donnees
betaVrai=1.5
etaVrai=100
c=40
n=5000
tVrai <- rweibull(n,betaVrai,etaVrai)
tVrai[tVrai>c]=c

SEM(tVrai,c)



#On prend eta0 et beta0 arbitrairement
eta0=50
beta0=3

#on stocke la suite des eta_k, beta_k dans une matrice
res=c(beta0,eta0)

n=length(tVrai)
m=length(tVrai[tVrai<c])   #nb de deffaillances observees

V=logVrai(c,eta0,beta0,tVrai)

#un simple compteur
i=1
Epsilon=0.2

repeat
{
  #on simule ce qu'on a pas pu observer
  tVrai[tVrai>=c]=simu(c,eta0,beta0,n-m)
  #calcul des EMV correspondants
  beta0=EMV_SEM(tVrai)$beta_EMV_SEM
  eta0=EMV_SEM(tVrai)$eta_EMV_SEM
  res=rbind(res,c(beta0,eta0))
  V=c(V,logVrai(c,eta0,beta0,tVrai))
  if (abs(V[i+1]-V[i])<Epsilon) {break}
  if (i==100) {break}
  i=i+1
}

x=seq(from=1,to=i+1,by=1)
par(mfrow = c(1,3))
plot(x,res[,1],xlab="nb iterations",ylab="beta estime")
abline(h=betaVrai,col="red")
plot(x,res[,2],xlab="nb iterations",ylab="eta estime")
abline(h=etaVrai,col="red")
plot(x,V)


#a faire : tracer la vraisemblance pour chaque couple (beta_k, eta_k)
#une fois qu'on aura un crit?re d'arr?t correct, on peut voir la convergence de la m?thode en tracant 
#en fonction de n les estimations obtenues



