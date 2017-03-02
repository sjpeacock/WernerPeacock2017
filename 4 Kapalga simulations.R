rm(list=ls())
library(lattice)
library(gplots)

# setwd("~/Google Drive/Plant Matrix model/Draft/WernerPeacock2017")
source('functions.R', chdir = TRUE)

A<-ImportMatrices()

######################################################################################
# Predictions from 1982 to 2003
######################################################################################

# Plots were in NON_sorghum (all)
s<-2

# Starting distirbution in 1982
tot<-c(0, 208, 162, 1, 59, 11, 456, 21)
init<-cbind(
	c(0,71,56,0,42,8,133,5),
	c(0,43,34,0,8,1,98,7),
	c(0,26,20,1,5,1,56,6),
	c(0,24,19,0,6,1,52,5)
	)

	
dat2003<-cbind(
	c(NA, NA, NA, 9,0,7,82,4),
	c(NA, NA, NA,7,1,6,60,2),
	c(NA, NA, NA, 1,4,6,20,5),
	c(NA, NA, NA, 1,1,2,34,1)
	)
	
dat1989<-cbind(
	c(NA, NA, NA, 0,0,1,82,2),
	c(NA, NA, NA,0,1,0,60,1),
	c(NA, NA, NA, 0,2,1,20,3),
	c(NA, NA, NA, 0,1,0,33,1)
	)

plots<-c("4S", "4N", "3S", "3N")

# Series of burns for four different plots
burned<-cbind(c(1,1,1,1,3,4,1,2,2,2,2,2,2,2,1,1,2,1,1,1,1),
	c(1,1,3,1,1,2,1,2,2,2,2,2,2,2,1,1,2,1,1,1,1),
	c(1,1,3,1,3,1,3,2,2,2,2,3,1,2,2,1,2,2,2,1,1),
	c(1,1,1,1,1,2,3,2,2,2,2,3,1,2,2,1,2,2,2,1,1), 
	c(1983:2003))
# 1= no fire,2= early fire, 3= late fire, 

nt<-21 # number of timesteps

# Matrix of 4 plots, 8 stage, timesteps (nt)
N<-matrix(); length(N)<-4*8*(nt+1); dim(N)<-c(4,8,(nt+1)) 
for(i in 1:4) N[i,,1]<-init[,i]

# Simulation
for(i in 1:4){	# for each plot	
	for(t in 1:nt){ # for each timestep
		
		# Deterministic
		N[i,,t+1]<-A[[burned[t,i], s]]%*%N[i,,t]
	
	}
}



lambda<-numeric(4)
for(i in 1:4){
	U<-A[[burned[1,i],s]]
	for(t in 2:nt){
	 U<-A[[burned[t,i],s]]%*%U
	 }
	lambda[i]<-calcLambda(U)^(1/nt)
}

# # Write table 6 
burned2<-burned
burned2[burned2==1]<-"N"
burned2[burned2==2]<-"E"
burned2[burned2==3]<-"L"
burned2[burned2==4]<-"W"

# write.csv(burned2[,c(5,1,2,3,4)], file="Table6.csv")
# write.csv(rbind(c("Stage", plots), cbind(c("1. SEEDLING", "2. SMALL JUVENILE", "3. LARGE JUVENILE", "4. SMALL SAPLING", "5. LARGE SAPLING", "6. POLE", "7. ADULT", "8. LARGE ADULT"), init)), file="Table7.csv")

