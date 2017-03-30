rm(list=ls())
# setwd("~/Google Drive/Plant Matrix model/Draft/WernerPeacock2017")
source('functions.R')

A<-ImportMatrices()


library(lattice)
library(gplots)

#p.burn<-matrix(c(0,0,1,1,1,1,1,1), nrow=4, ncol=2, byrow=TRUE)
# Werner (NTN) burn data
p.burn<-matrix(c(0,0,0.947,0.531,1,1,0.986,0.932), nrow=4, ncol=2, byrow=TRUE)

recruitment.nofire<-c(1,2,2)
seedlingSurv.nofire<-0.5
for(i in 1:4){
	for(j in 1:2){
		A[[i,j]][1,6:8]<-recruitment.nofire*(1-p.burn[i,j])
		A[[i,j]][2,1]<-seedlingSurv.nofire*(1-p.burn[i,j])
	}
}

######################################################################################
######################################################################################

scenarios<-data.frame(burn.interval=rep(c(2,4,3,5,"2/5"), 2), fire.type=rep(c(rep(c("early", "late"), each=2), "early/late"), 2), fire.sequence=rep(c("EN", "ENNN", "LNN", "LNNNN", "NENEL"), 2), sorghum=rep(c("Y", "N"), each=5), lambda.det=rep(NA, 10), lambda.stoch.mean=rep(NA, 10), lambda.stoch.LI=rep(NA, 10), lambda.stoch.UI=rep(NA, 10), lambda.stoch.pop=rep(NA, 10))

N.all<-list();length(N.all)<-dim(scenarios)[1]

set.seed(327854)
for(m in 1:dim(scenarios)[1]){
	if(scenarios$fire.sequence[m]=="NENEL") stoch.probs.<-rbind(c(0.4, 0.2, 0.4), c(0.5, 0.1, 0.4)) else stoch.probs.=NA
	N.all[[m]]<-sim.N(sorghum=scenarios$sorghum[m], burn.sequence=as.character(scenarios$fire.sequence[m]), stoch.probs=stoch.probs., nt=500, n=1000, A.=A, init=NA, make.smaller=TRUE)
	}

for(m in 1:dim(scenarios)[1]) scenarios[m,5:9]<-N.all[[m]][[2]]

Table6<-scenarios
Table6$lambda.stoch<-paste(scenarios$lambda.stoch.mean, " (", scenarios$lambda.stoch.LI, ", ", scenarios$lambda.stoch.UI, ")", sep="")
Table6$lambda.stoch.range<-	paste(scenarios$lambda.stoch.LI, "-", scenarios$lambda.stoch.UI, sep="")

# write.csv(Table6, file="Table6.csv")
# save.image("scenarios_20170310.RData")