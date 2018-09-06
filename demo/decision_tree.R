rm(list=ls())

library(rpart)

magasin=read.table(here::here("data","magasin.csv"),
                   sep=";", header=T, dec=",")
attach(magasin)
dt=rpart(Achat~Prix+Flashitude+Branchitude+Qualite, data=magasin)

#plot(dt)
#text(dt, use.n=T)
print(dt)

calcul_impurete=function(x){
  probas=table(x)/length(x)
  impurete=sum(probas*(1-probas))
  return(impurete)
}

probas=table(Achat)/length(Achat)
impurete=calcul_impurete(Achat)

prA=1
IA=calcul_impurete(Achat)
print(IA)

n=length(Achat)
indL=which(Branchitude<85)
indR=which(Branchitude>=85)
prL=length(indL)/n
prR=length(indR)/n
IL=calcul_impurete(Achat[indL])
IR=calcul_impurete(Achat[indR])

perte_impurete=prA*IA-prL*IL-prR*IR
print(perte_impurete)

prA=1
indL=which(Branchitude<65)
indR=which(Branchitude>=65)

prL=length(indL)/n
prR=length(indR)/n
IL=calcul_impurete(Achat[indL])
IR=calcul_impurete(Achat[indR])

perte_impurete=prA*IA-prL*IL-prR*IR
print(perte_impurete)

prA=1
indL=which(Flashitude<65)
indR=which(Branchitude>=65)

n=length(Achat)
prL=length(indL)/n
prR=length(indR)/n
IA=calcul_impurete(Achat)
IL=calcul_impurete(Achat[indL])
IR=calcul_impurete(Achat[indR])

perte_impurete=prA*IA-prL*IL-prR*IR
print(perte_impurete)

## [1] -0.1247