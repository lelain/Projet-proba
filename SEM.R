#Methode SEM

#Les donnees
betaVrai=1.5
etaVrai=100
c=40
n=500
t <- rweibull(n,betaVrai,etaVrai)
t[t>c]=c
m=length(t[t<c])



#fonction permettant de simuler N realisations de loi de Weibull de parametre eta0 et beta0, conditionnee par X>c
simu=function(c,eta0,beta0,N)
{
  u=runif(N)
  return(eta0*((-log(u)+(c/eta0)^beta0)^(1/beta0)))
}



f_beta_SEM=function(b,ech,N)
{
  1/b + (sum(log(ech))/N) - ((sum((ech^b)*log(ech)))/(sum(ech^b)))
}


#fonction calculant l'estimateur de max de vraisemblance d'une 
#loi de Weibull de parametre eta0 et beta0, conditionnee par X>c
EMV_SEM=function(t)
{
  n=length(t)
  beta.est = uniroot(f_beta_SEM,lower=0.001,upper=10,ech=t,N=n)$root
  eta.est = (1/n * sum(t^beta.est))^(1/beta.est)
  return(list(beta_EMV_SEM=beta.est,eta_EMV_SEM=eta.est))
}

#On prend eta0 et beta0 arbitrairement
eta0=50
beta0=3

#nb d'iteration de la methode
K=1000

#on stocke la suite des eta_k, beta_k dans une matrice
res=matrix(nrow = K+1,ncol = 2)
res[1,1]=beta0
res[1,2]=eta0

n=length(t)
m=length(t[t<c])   #nb de deffaillances observees
for (k in 1:K)
{
  #on simule ce qu'on a pas pu observer
  t[t>=c]=simu(c,eta0,beta0,n-m)
  #calcul des EMV correspondants
  beta0=EMV_SEM(t)$beta_EMV_SEM
  eta0=EMV_SEM(t)$eta_EMV_SEM
  res[k+1,1]=beta0
  res[k+1,2]=eta0
}

x=seq(from=1,to=K+1,by=1)
par(mfrow = c(1,2))
plot(x,res[,1],xlab="nb iterations",ylab="beta estime")
abline(h=betaVrai,col="red")
plot(x,res[,2],xlab="nb iterations",ylab="eta estime")
abline(h=etaVrai,col="red")


#a faire : tracer la vraisemblance pour chaque couple (beta_k, eta_k)
#une fois qu'on aura un critère d'arrêt correct, on peut voir la convergence de la méthode en tracant 
#en fonction de n les estimations obtenues



