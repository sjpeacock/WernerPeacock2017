library(here)
library(lattice)
library(gplots)

source('functions.R')

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

quartz(width=6.3, height=6)

# colK<-c(colF, "orange")
layout(matrix(c(1:8), nrow=4, ncol=2, byrow=TRUE))
	par(mar=c(2,2,2,0), oma=c(2,2,3,1), mgp=c(2.5, 0.8, 0))
	for(i in 1:8){
		plot(1:(nt+1), N[1,i,1:(nt+1)], "n", bty="l", xaxt="n", las=1, xlab="", ylab="", ylim=c(0,1.05*max(c(N[,i,1:(nt+1)], dat2003[i,], dat1989[i,]), na.rm=TRUE)), xaxt="n")
		for(j in 1:4){
			points(1:(nt+1), N[j,i,1:(nt+1)], "o", pch=20+j, bg="white")#, col=colK[j])
			points(nt+1, dat2003[i,j], pch=20+j, col="#00000050", bg="#00000050")#, col=colK[j])
			points(8, dat1989[i,j], pch=20+j, col="#00000050",  bg="#00000050")#, col=colK[j])
			# lines(1:(nt+1), N[1,i,1:(nt+1)], col=colK[j])
			}
		
		axis(side=1, at=seq(1,(nt+1),1), labels=FALSE, tck=-0.05)
		axis(side=1, at=seq(3,(nt+1),5), labels=seq(1985, 2003,5), tck=-0.1, mgp=c(2.5, 0.5,0))
		
		mtext(side=3, line=0.5, paste(" ", letters[i], ") ", c("Seedlings", "Small juveniles", "Large juveniles", "Small saplings", "Large saplings", "Poles", "Adults", "Large adults")[i], sep=""), cex=par("cex"), adj=0)
		
		if(i==1){
		legend(15, 2.1*max(c(N[,i,1:(nt+1)], dat2003[i,], dat1989[i,]), na.rm=TRUE), col=c(1, "#00000050"), pch=c(21, 19), pt.bg="white", c("Prediction", "Observation"), bty="n", xpd=NA, bg="white", lwd=c(1, NA))
		legend(25, 2.1*max(c(N[,i,1:(nt+1)], dat2003[i,], dat1989[i,]), na.rm=TRUE), pch=21:25, legend=plots, border=NA, bty="n", xpd=NA, bg="white", ncol=3)
			}
		}
mtext(side=2, outer=TRUE, "Number of individuals",line=0.5)
mtext(side=1, outer=TRUE, "Year",line=0.5)
