rm(list=ls())

setwd('~/Google Drive/Plant Matrix model/Draft/WernerPeacock2017')

source('functions.R', chdir = TRUE)

A<-ImportMatrices()


library(lattice)
library(gplots)
library(popbio)
#p.burn<-matrix(c(0,0,1,1,1,1,1,1), nrow=4, ncol=2, byrow=TRUE)
# Werner (NTN) burn data
p.burn<-matrix(c(0,0,0.947,0.531,1,1,0.986,0.932), nrow=4, ncol=2, byrow=TRUE)


recruitment.nofire<-cbind(1, 2, 10)
seedlingSurv.nofire<-c(0.1, 0.5, 0.9)
x<-c(1:20)

nt<-500
n<-1000

L.recruit<-list(); length(L.recruit)<-6*3; dim(L.recruit)<-c(3,3,2)
L.seed<-L.recruit

for(i in 1:3){
	for(j in 1:2){
		for(k in 1:3){
			L.recruit[[k,i,j]]<-matrix(NA, length(x), 4)
			L.seed[[k,i,j]]<-matrix(NA, length(x), 4)
		}}}

#------------------------------

time0<-proc.time()[3]
for(u in 1:2){ # First recruitment, then seedling survival
	for(m in 1:3){	# for each level of the parameter (low, moderate, high)
		
		# Create matrices for the simulation
		for(i in 1:4){
			for(j in 1:2){
				
				if(u==1){
					A[[i,j]][1,6:8]<-c(0.5, 1, 1)*recruitment.nofire[m]*(1-p.burn[i,j])
					A[[i,j]][2,1]<-0.5*(1-p.burn[i,j]) #baseline moderate seedling surv
				}else if(u==2){
					A[[i,j]][1,6:8]<-c(1,2,2)*(1-p.burn[i,j]) #baseline moderate recruitment
					A[[i,j]][2,1]<-seedlingSurv.nofire[m]*(1-p.burn[i,j]) 
				}
			}# end j
		}# end i

		for(f in 1:3){ # for each fire type except no fire
			for(s in 1:2){ #for sorghum and non-sorghum
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
					N[[1]][,1] <- init
					N[[2]][,,1] <- init
		
					# Simulation
					for(t in 1:nt){ # for each timestep						
						# Deterministic
						N[[1]][,t+1]<-A[[burned[[1]][t],s]]%*%N[[1]][,t]					
						# Stochastic
						for(i in 1:n) N[[2]][i,,t+1]<-A[[burned[[2]][t,i],s]]%*%N[[2]][i,,t]
						} #end t

				# Stochastic lambda
				log.lambda_s<-numeric(n)
				for(i in 1:n){
					LLsik<-numeric(10)
					for(t in (nt-10):(nt-1)) LLsik[t-(nt-10)+1]<-log(sum(N[[2]][i,,t+1]))-log(sum(N[[2]][i,,t]))
					log.lambda_s[i]<-mean(LLsik)
				} #end i

	A1<-A[[burn.seq.num[1],s]]
	if(x[j]>1){for(i in 1:(x[j]-1)) A1<-A[[1,s]]%*%A1}
	
	lambda.all<-c(det=round(calcLambda(A1)^(1/x[j]), 3),
		stoch.mean=round(exp(mean(log.lambda_s)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.025)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.975)), 3)
		)

	if(u==1) L.recruit[[m,f,s]][j,]<-lambda.all
	if(u==2) L.seed[[m,f,s]][j,]<-lambda.all
	
	} #end j fire ret int
} # end s sorg and non sorg
} # end f fire type

} # end m
} # end u
cat("run time: ", proc.time()[3]-time0)

#-----------------------------------------
m<-1 # low, moderate, high
f<-1 # early, late, wet
s<-1 # sorghum, non-sorghum

#Recruitment effect
data.frame(RetInt=x, lambda_det=L.recruit[[m,f,s]][,1], lambda_stoch=L.recruit[[m,f,s]][,2])

# Seedling survival
data.frame(RetInt=x, lambda_det=L.seed[[m,f,s]][,1], lambda_stoch=L.seed[[m,f,s]][,2])


MinFRT.recruit<-data.frame(recruitment=rep(c("low", "moderate", "high"), each=3*2), fire=rep(rep(c("early", "late", "wet"), each=2), 3), understory=rep(c("S", "NS"), 3*3), MinFRT_det=rep(NA, 3*3*2), MinFRT_stoch=rep(NA, 3*3*2))
MinFRT.seed<-MinFRT.recruit; names(MinFRT.seed)[1]<-"seedling.surv"

for(m in 1:3){
	for(f in 1:3){
		for(s in 1:2){
			
			ind<-which(MinFRT.recruit[,1]==c("low", "moderate", "high")[m]&MinFRT.recruit$fire==c("early", "late", "wet")[f]&MinFRT.recruit$understory==c("S", "NS")[s])
			
			if(f==1&s==2){
				MinFRT.recruit$MinFRT_det[ind]<-x[which(L.recruit[[m,f,s]][,1]<1)[1]-1]
				MinFRT.recruit$MinFRT_stoch[ind]<-x[which(L.recruit[[m,f,s]][,2]<1)[1]-1]
				
				MinFRT.seed$MinFRT_det[ind]<-x[which(L.seed[[m,f,s]][,1]<1)[1]-1]
				MinFRT.seed$MinFRT_stoch[ind]<-x[which(L.seed[[m,f,s]][,2]<1)[1]-1]
			}else{
				MinFRT.recruit$MinFRT_det[ind]<-x[which(L.recruit[[m,f,s]][,1]>1)[1]]
				MinFRT.recruit$MinFRT_stoch[ind]<-x[which(L.recruit[[m,f,s]][,2]>1)[1]]
				
				MinFRT.seed$MinFRT_det[ind]<-x[which(L.seed[[m,f,s]][,1]>1)[1]]
				MinFRT.seed$MinFRT_stoch[ind]<-x[which(L.seed[[m,f,s]][,2]>1)[1]]
			}
		}
	}
}

write.csv(MinFRT.recruit, file="MinFRT_recruitment.csv")
write.csv(MinFRT.seed, file="MinFRT_seedlingSurvival.csv")
