#comparaison des deux methodes

betaVrai=1.2
etaVrai=100
c=40
n=25
k=50 #nb d'echantillons simules
M=c()
EMV.est=c()
SEM.est=c()
for (i in 1:k)
{
  repeat
  {
    t=rweibull(n,betaVrai,etaVrai)
    t[t>c]=c
    if (length(t[t<c])!=0) {break}  #pour s'assurer de ne pas avoir m=0
  }
  M=c(M,length(t[t<c]))
  EMV.est=rbind(EMV.est,c(EMV(t,c)$beta_EMV,EMV(t,c)$eta_EMV))
  SEM.est=rbind(SEM.est,c(SEM(t,c)$beta_SEM,SEM(t,c)$eta_SEM))
}

Res=matrix(ncol=4,nrow=2)
Res[1,1]=mean(EMV.est[,1])
Res[2,1]=mean(EMV.est[,2])
Res[1,2]=sd(EMV.est[,1])
Res[2,2]=sd(EMV.est[,2])
Res[1,3]=mean(SEM.est[,1])
Res[2,3]=mean(SEM.est[,2])
Res[1,4]=sd(SEM.est[,1])
Res[2,4]=sd(SEM.est[,2])