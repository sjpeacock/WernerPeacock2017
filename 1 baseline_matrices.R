library(lattice)
library(matrixcalc)

rm(list=ls())

source('~/Google Drive/Plant Matrix model/Draft/Code/functions.R', chdir = TRUE)

A<-ImportMatrices()

# Werner (NTN) burn data
p.burn<-matrix(c(0,0,0.947,0.531,1,1,0.986,0.932), nrow=4, ncol=2, byrow=TRUE)
colnames(p.burn)<-c("sorghum", "non-sorghum")
rownames(p.burn)<-c("no fire", "early", "late", "wet")

# Re-calculate recruitment and seedling survival based on p.burn
recruitment.nofire<-c(1,2,2)
for(i in 1:4){
	for(j in 1:2){
		A[[i,j]][1,6:8]<-recruitment.nofire*(1-p.burn[i,j])
		A[[i,j]][2,1]<-0.5*(1-p.burn[i,j])
	}
}

######################################################################################
# Stable stage
######################################################################################


plot.SS(A)

######################################################################################
# Sensitivity
######################################################################################

S<-list(); length(S)<-8; dim(S)<-c(4,2)
M<-0
S0<-S
valid<-matrix(c(0,1,rep(0,7), rep(1,4), rep(0,4), rep(1,4), rep(0,4), rep(1,5), rep(0,3), rep(1,5), rep(0,2), 1, rep(0,3), rep(1,3), 0, 1, rep(0,4), rep(1,4), rep(0,5), 1, 1), nrow=8, ncol=8, byrow=FALSE)
for(i in 1:4){
	for(j in 1:2){
		S0[[i,j]]<-sensitivity(A[[i,j]])
		S[[i,j]]<-S0[[i,j]]*valid
		M<-max(c(M, S[[i,j]]))
		
		}
	}


quartz(width=5, height=6)
layout(matrix(c(rep(c(1:4), 3), rep(c(5:8), 3), rep(9,4)), nrow=4, ncol=7,byrow=FALSE))
par(mar=c(1,1,0,0), oma=c(5,4,2,5))
for(j in 1:2){
	for(i in 1:4){
	plot(0:8+0.5, 0:8+0.5, "n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n")
		if(j==1) axis(side=2, at=1:8, labels=rev(c("sdlg", "smJuv", "Juv", "smSap", "sap", "pole", "adult", "lgAdult")), las=1)
		if(i==4) axis(side=1, at=1:8, labels=c("sdlg", "smJuv", "Juv", "smSap", "sap", "pole", "adult", "lgAdult"), las=2)
		for(f in 1:8){
			for(t in 1:8){
				Sft<-round(S[[i,j]][t,f], 5)
				if(Sft==0) polygon(x=f+c(-0.5, -0.5, 0.5, 0.5), y=c(8:1)[t]+c(-0.5, 0.5, 0.5, -0.5), border=1, col="white", lwd=0.1) else polygon(x=f+c(-0.5, -0.5, 0.5, 0.5), y=c(8:1)[t]+c(-0.5, 0.5, 0.5, -0.5), border=1, col=rev(cols(100))[round(100*(Sft/M))], lwd=0.1)#grey(1-Sft/sum(S[[i,j]]))
				
			}
		}
		}}
mtext(side=2, outer=TRUE, line=2.5,"To")
mtext(side=1, outer=TRUE, line=4,"From")
#mtext(side=3, "Sensitivity", font=2, outer=TRUE, line=0.5)

par(mar=c(0,1,0,0))
plot(1,1,"n",xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,100))
for(i in 1:100) polygon(c(0,1,1,0), c(i-1, i-1, i, i), border=NA, col=rev(cols(100))[i])
axis(side=4, at=seq(0,100*1.54/1.5,30), labels=seq(0,1.5, 0.5), las=1)
mtext(side=4, "Senstivity", line=2)

S2<-round(rbind(
	cbind(S[[1,1]], S[[1,2]]),
	cbind(S[[2,1]], S[[2,2]]),
	cbind(S[[3,1]], S[[3,2]]),
	cbind(S[[4,1]], S[[4,2]])),3)
	
write.csv(S2, "sensitivity.csv")
######################################################################################
# Elasticity
######################################################################################



E<-list(); length(E)<-8; dim(E)<-c(4,2)
for(i in 1:4){
	for(j in 1:2){
		E[[i,j]]<-elasticity(A.all[[i,j]])
		}
	}

par(mfrow=c(4,2), mar=c(1,1,0,0), oma=c(5,4,2,1))
for(i in 1:4){
	for(j in 1:2){
		plot(0:8+0.5, 0:8+0.5, "n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n")
		if(j==1) axis(side=2, at=1:8, labels=rev(c("sdlg", "smJuv", "Juv", "smSap", "sap", "pole", "adult", "lgAdult")), las=1)
		if(i==4) axis(side=1, at=1:8, labels=c("sdlg", "smJuv", "Juv", "smSap", "sap", "pole", "adult", "lgAdult"), las=2)
		for(f in 1:8){
			for(t in 1:8){
				if(round(E[[i,j]][f,t], 5)==0) polygon(x=f+c(-0.5, -0.5, 0.5, 0.5), y=c(8:1)[t]+c(-0.5, 0.5, 0.5, -0.5), border=1, col="#FF000020", lwd=0.1) else polygon(x=f+c(-0.5, -0.5, 0.5, 0.5), y=c(8:1)[t]+c(-0.5, 0.5, 0.5, -0.5), border=1, col=grey(1-round(E[[i,j]][f,t], 5)), lwd=0.1)
				
			}
		}
		}}
mtext(side=2, outer=TRUE, line=2.5,"To")
mtext(side=1, outer=TRUE, line=4,"From")
mtext(side=3, "Elasticity", font=2, outer=TRUE, line=0.5)
#------------------
