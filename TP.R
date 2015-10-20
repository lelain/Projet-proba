#Estimation de durees de vie censurees

curve(dweibull(x,shape=0.5,scale=2),xlim=c(0,8),ylim=c(0,1.3),ylab="f(x)",main="Diverses densités de loi de Weibull et taux de défaillance")
curve(dweibull(x,shape=1,scale=3),add=T,col="green")
curve(dweibull(x,shape=5,scale=1.5),add=T,col="red")

curve((0.5/2)*((x/2)^(0.5-1)),add=T,lty=3)
curve((1/3)*((x/3)^(1-1)),add=T,lty=3,col="green")
curve((5/1.5)*((x/1.5)^(5-1)),add=T,lty=3,col="red")

legend("topright",lty=c(1,1,1,3,3,3),col=c("black","green","red","black","green","red"),legend=c("beta=0.5,eta=2","beta=1,eta=3","beta=5,eta=1.5","beta=0.5,eta=2,défaillance","beta=1,eta=3,défaillance","beta=5,eta=1.5,défaillance"))

#simuler un échantillon 
beta=1.5
eta=100
c=40
n=50
t <- rweibull(n,beta,eta)
t[t>c]=c

#M = nombre de pannes qu'on a observés, loi binomiale(n,1-G(c)), avec G(c)=P(X<=c)=F(c)
#X loi de weibull, T loi de weibull censurée



f_beta=function(beta,t,m,c)
{
  1/beta + (sum(log(t[t<c]))/m) 
            - ((sum((t^beta)*log(t)))/(sum(t^beta)))
}

EMV=function(t,c)
{
  m=length(t[t<c])
  beta.est = uniroot(f_beta,lower=0.001,upper=10,t=t,m=m,c=c)$root
  eta.est = (1/m * sum(t^beta.est))^(1/beta.est)
  return(list(beta_EMV=beta.est,eta_EMV=eta.est))
}

EMV(t,c)

#estimation de la variance et du biais
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



