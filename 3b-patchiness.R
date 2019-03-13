rm(list=ls())
# setwd("~/Google Drive/Plant Matrix model/Draft/WernerPeacock2017")
source('functions.R')

A<-ImportMatrices()
		
library(lattice)
library(gplots)
library(popbio)

recruitment.nofire<-3
seedlingSurv.nofire<-0.5
x<-c(1:20)

p.burn2<-c(0.25, 0.5, 0.75, 1)


nt<-500
n<-1000

L.patchy<-list(); length(L.patchy)<-length(p.burn2)*2*2; dim(L.patchy)<-c(length(p.burn2),2,2) #Just for early dry season fires


for(j in 1:2){
	for(i in 1:length(p.burn2)){
		for(k in 1:2)L.patchy[[i,j,k]]<-matrix(NA, length(x), 4)
		}}

#------------------------------

time0<-proc.time()[3]
# f <- 1 	# Just consider early dry season fires 
f <- 3 	# Just consider wet season fires
for(m in 1:length(p.burn2)){
	for(k in 1:2){ #moderate and high seedling survival
		for(s in 1:2){ #for sorghum and non-sorghum
		
		A[[1,s]][1,6:8]<-c(1,2,2)
		A[[1,s]][2,1]<-c(0.5, 0.9)[k]
		
		A[[f+1,s]][1,6:8]<-A[[1,s]][1,6:8]*(1-p.burn2[m])
		A[[f+1,s]][2,1]<-A[[1,s]][2,1]*(1-p.burn2[m])
		
		for(j in 1:length(x)){
				
			burn.seq.num<-c(f+1, rep(1, x[j]-1))
			
			# Stochastic series of burns
			burned<-list(
				deterministic=rep(burn.seq.num, ceiling(nt/x[j])), 
				stochastic=matrix(rbinom(nt*n, size=1, prob=1/x[j]), nrow=nt, ncol=n)*f+1)
	
	

			# Matrix of 8 stage, sorg/non sorg (2), deterministic/stochastic (2), timesteps (nt)
			N<-list(); length(N)<-2
			N[[1]]<-matrix(); length(N[[1]])<-8*(nt+1); dim(N[[1]])<-c(8,(nt+1)) #Deterministic
			N[[2]]<-matrix(); length(N[[2]])<-n*8*(nt+1); dim(N[[2]])<-c(n,8,(nt+1)) #Stochastic
			
			# Initial condition: Starting distirbution is the stable stage for no fire
			init<-1000*stableStage(A[[1,s]])
			N[[1]][,1]<-init
			for(i in 1:n){
				N[[2]][i,,1]<-init
				}

			# Simulation
			for(t in 1:nt){ # for each timestep
				
				# Deterministic
				N[[1]][,t+1]<-A[[burned[[1]][t],s]]%*%N[[1]][,t]
				
				# Stochastic
				for(i in 1:n){
					N[[2]][i,,t+1]<-A[[burned[[2]][t,i],s]]%*%N[[2]][i,,t]
					} #end i
				
			}

		# Stochastic lambda
		log.lambda_s<-numeric(n)
		for(i in 1:n){
			LLsik<-numeric(10)
			for(t in (nt-10):(nt-1)){
				LLsik[t-(nt-10)+1]<-log(sum(N[[2]][i,,t+1]))-log(sum(N[[2]][i,,t]))
			}
			log.lambda_s[i]<-mean(LLsik)
		}
	

		A1<-A[[burn.seq.num[1],s]]
		if(x[j]>1){for(i in 1:(x[j]-1)) A1<-A[[1,s]]%*%A1}
		
		lambda.all<-c(det=round(calcLambda(A1)^(1/x[j]), 3),
		stoch.mean=round(exp(mean(log.lambda_s)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.025)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.975)), 3))

		L.patchy[[m,s,k]][j,]<-lambda.all

		} #end j fire ret int
	} # end s sorg and non sorg
	} # end k seelding survival
} # end m patchiness

cat("run time: ", proc.time()[3]-time0) #688 mins

#------------------
MinFRT.patch<-data.frame(patchiness=rep(c(p.burn2), each=2*2), fire=rep("early", 4*2*2), seedling.surv=rep(rep(c("moderate", "high"), each=2), 4), understory=rep(c("S", "NS"), 4*2), MinFRT_det=rep(NA, 4*2*2), MinFRT_stoch=rep(NA, 4*2*2))

for(m in 1:length(p.burn2)){
	for(f in 1:2){
		for(s in 1:2){
			
			ind<-which(MinFRT.patch[,1]==p.burn2[m]&MinFRT.patch$seedling.surv==c("moderate", "high")[k]&MinFRT.patch$understory==c("S", "NS")[s])
			
			MinFRT.patch$MinFRT_det[ind]<-x[which(L.patchy[[m,s,k]][,1]>1)[1]]
			MinFRT.patch$MinFRT_stoch[ind]<-x[which(L.patchy[[m,s,k]][,2]>1)[1]]
		}
	}
}

#------------------
save.image("patchiness_wet_20170312.RData")