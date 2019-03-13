library(here)
library(lattice)
library(gplots)
library(caTools)

source('functions.R')
source('SimN_revised.R')

# Color palette for plotting stages
colF<-c("#2c7bb6", "#FF6347", "#228B22")

# Load transition matrices
A<-ImportMatrices()


# Werner (NTN) burn data for proportion of ground burned in each fire type
p.burn<-matrix(c(0,0,0.947,0.531,1,1,0.986,0.932), nrow=4, ncol=2, byrow=TRUE)

# Adjust recruitment and seeding survival based on p.burn
recruitment.nofire<-c(1,2,2)
seedlingSurv.nofire<-0.5
for(i in 1:4){
	for(j in 1:2){
		A[[i,j]][1,6:8]<-recruitment.nofire*(1-p.burn[i,j])
		A[[i,j]][2,1]<-seedlingSurv.nofire*(1-p.burn[i,j])
	}
}


######################################################################################
# Environment transition matrices
######################################################################################

env.mat_all<-list(
	EN.neg=matrix(c(0.25, 0.75, 0.75, 0.25), 2, 2),
	EN.0=matrix(rep(0.5, 4), 2, 2),
	EN.pos=matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2),
	ENNN.neg=matrix(c(0, 1, 1/3, 2/3), 2, 2),
	ENNN.0=matrix(c(0.25, 0.75, 0.25, 0.75), 2, 2),
	ENNN.pos=matrix(c(5/12, 7/12, 1/6, 5/6), 2, 2),
	LNN.neg=matrix(c(0, 1, 0.5, 0.5), 2, 2),
	LNN.0=matrix(c(1/3, 2/3, 1/3, 2/3), 2, 2),
	LNN.pos=matrix(c(2/3, 1/3, 1/6, 5/6), 2, 2),
	LNNNN.neg=matrix(c(0, 1, 0.25, 0.75), 2, 2),
	LNNNN.0=matrix(c(1/5, 4/5, 1/5, 4/5), 2, 2),
	LNNNN.pos=matrix(c(2/5, 3/5, 3/20, 17/20), 2, 2),
	NENEL.neg=matrix(c(0, 0.5, 0.5, 0, 0, 1, 1, 0, 0), 3,3),
	NENEL.0=matrix(c(rep(0.4, 3), rep(0.2, 3), rep(0.4, 3)), 3,3),
	NENEL.pos=NA
	)
dim(env.mat_all)<-c(3,5)

# Define parameters for simulation:
nt<-100 # Number of timesteps (transient dynamics - short 100 years)
n<-1000 # Number of MC trials for each scenario

# Transition matrices for each scenario
A._all<-list(
	list(A[[2,1]], A[[1,1]]), list(A[[2,2]], A[[1,2]]),
	list(A[[2,1]], A[[1,1]]), list(A[[2,2]], A[[1,2]]),
	list(A[[3,1]], A[[1,1]]), list(A[[3,2]], A[[1,2]]),
	list(A[[3,1]], A[[1,1]]), list(A[[3,2]], A[[1,2]]),
	list(A[[2,1]], A[[3,1]], A[[1,1]]), list(A[[2,2]], A[[3,2]], A[[1,2]]))
dim(A._all)<-c(2,5)

burn.seq_all=list(c(1,2), c(1, 2, 2, 2), c(1, 2, 2), c(1,2,2,2,2), c(3,1,3,1,2))
	

######################################################################################
# Simulation
######################################################################################

# Transient dynamics are extremely sensitive to starting distribution
# We start with the stable stage for early-dry season fires in sorghum and 
# Stable-stage for no fire in non-sorghum.

start.state<-rbind(1000*stableStage(A[[2,1]]), 1000*stableStage(A[[1,2]]))

# Create lists to store output of all simulations:
N_det<-list(); length(N_det)<-10; dim(N_det)<-c(2,5)
N_all<-list(); length(N_all)<-10; dim(N_all)<-c(2,5)
for(i in 1:2){for(j in 1:5){N_all[[i,j]]<-list(stoch.neg=NA, stoch.0=NA, stoch.pos=NA)}}

t0<-proc.time()[3] # takes ~1 minute
set.seed(32784) # Change to run in parallel?
for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
		N_det[[i,j]]<-sim.det(nt=nt, burn.seq=burn.seq_all[[j]], A.=A._all[[i,j]], start.state = start.state[i,])
		for(k in 1:3){ # For 3 different autocorrelations
			if(j==5&k==3){}else{
				N_all[[i,j]][[k]]<-sim.stoch_n(env.matrix=env.mat_all[[k,j]], nt=nt, n=n, A.=A._all[[i,j]], start.state=start.state[i,], transient=TRUE, sorghum=c(TRUE, FALSE)[i]) 
				}
			}#end k
	} # end j scenarios
}# end sog/nonsorg
time.tot<-(proc.time()[3]-t0)/60
cat("Run time: ", time.tot, "mins")

# Calculate proportion of 1000 simulation that experienced a increase (avg. popn growth rate >=1) over the 100 years (Proportion Declining)
PI<-matrix(rep(NA,30)); dim(PI)<-c(2, 5, 3)
 for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
		for(k in 1:3){
			if(j==5&k==3){}else{
				PI[i,j,k]<-length(which(N_all[[i,j]][[k]]$log.lambda_s>=0))/n
}}}}

######################################################################################
# Figures and summary
######################################################################################

# #--------------------------------------------------------------------------------
# # Compare two starting distributions
# par(mfrow=c(2,1), mar=c(2,4,2,1), oma=c(3,0,0,0))
# barplot(start.state[1,], beside=TRUE, names.arg=1:8, ylab="Number of trees", ylim=c(0,700), xlab="Stage", main="Sorghum (Early dry season)")
# barplot(start.state[2,], beside=TRUE, ylab="Number of trees", ylim=c(0,700), xlab="Stage", main="Non-sorghum (No fire)")

#--------------------------------------------------------------------------------
# Table of avg. population change and PI
#--------------------------------------------------------------------------------
tableX<-data.frame(scenario=c(1:10), sorghum=c(rep("Y", 5), rep("N", 5)), avg_growth_rate=rep(NA, 10), PI=rep(NA, 10))
for(i in 1:2){
	for(j in 1:5){
		ind<-(i-1)*5+j
		
		tableX$avg_growth_rate[ind]<-paste(sprintf("%.3f", round(exp(mean(N_all[[i,j]][[2]]$log.lambda_s)), 3)), " (", sprintf("%.3f" , round(exp(quantile(N_all[[i,j]][[2]]$log.lambda_s, 0.025)), 3)), ", ", sprintf("%.3f", round(exp(quantile(N_all[[i,j]][[2]]$log.lambda_s, 0.975)), 3)), ")", sep="")
		
		tableX$PI[ind]<-sprintf("%.3f" , round(PI[i,j,2], 3))
		}}

write.csv(tableX, "TableX.csv")

#--------------------------------------------------------------------------------
# Figure of growth rate distributions for 10 scenarios (Fig. 6)
#--------------------------------------------------------------------------------

plot.corr<-FALSE
xlims<-matrix(rep(NA,20)); dim(xlims)<-c(2,5,2)
xlims[,1,1]<-c(-0.05, 0.06)
xlims[,1,2]<-c(-0.03, 0.02)
xlims[,2,1]<-c(-0.02, 0.08)
xlims[,2,2]<-c(-0.03, 0.01)
xlims[,3,1]<-c(-0.08, 0.04)
xlims[,3,2]<-c(-0.04, -0.005)
xlims[,4,1]<-c(-0.04, 0.06)
xlims[,4,2]<-c(-0.03, -0.005)
xlims[,5,1]<-c(-0.1, 0.02)
xlims[,5,2]<-c(-0.06, 0)
n.smooth<-50

quartz(width=5, height=6, pointsize=12)
par(mfcol=c(5,2), mar=c(2,1,2,0), oma=c(3,4,2,6))

for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
		
		xrange<-range(N_all[[i,j]][[2]]$log.lambda_s)
		
		h<-hist(N_all[[i,j]][[2]]$log.lambda_s, lab="", ylab="", main="", ylim=c(0,300), yaxt="n", breaks=seq(min(xrange), max(xrange), length.out=15), col=grey(0.8), yaxs="i")
		u<-par('usr');segments(u[1],u[3], u[1], u[4], xpd=NA);segments(u[1],u[3], u[2], u[3],xpd=NA)
		
		if(i==1) axis(side=2, at=c(0,100,200,300), las=1) else axis(side=2, at=c(0,100,200,300), labels=FALSE)
		# if(plot.corr==TRUE){
			# lines(density(N_all[[i,j]][[1]]$log.lambda_s, n= n.smooth), col="#FF000070")
			# if(j!=5){
				
				# lines(density(N_all[[i,j]][[3]]$log.lambda_s, n= n.smooth), col="#00FF0070")
			# }
			# }
		# arrows(0, 300, 0, 0, length=0.08, lwd=1.5)
		abline(v=0, lty=2, lwd=1.5)
		text(u[2], u[4], adj=1,paste("PI = ", sprintf("%.3f", round(PI[i,j,2], 3))), xpd=NA)
		
		
		mtext(side=3, adj=0, line=0, paste("  ", letters[(j-1)*2+i], ")", sep=""), cex=0.8)
		
		if(j==1&i==1) mtext(side=3, line=2, "Sorghum", cex=0.8)
		if(j==1&i==2) mtext(side=3, line=2, "Non-sorghum", cex=0.8)
		
	if(i==2&j==1){
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),180, "Early dry \nseason fire", xpd=NA)
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),80, expression(italic(x)==2), xpd=NA)
		}
		
	if(i==2&j==2){
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),180, "Early dry \nseason fire", xpd=NA)
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),80, expression(italic(x)==4), xpd=NA)
		}
		
	if(i==2&j==3){
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),180, "Late dry \nseason fire", xpd=NA)
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),80, expression(italic(x)==3), xpd=NA)
		}
		
	if(i==2&j==4){
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),180, "Late dry \nseason fire", xpd=NA)
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),80, expression(italic(x)==5), xpd=NA)
		}
		
	if(i==2&j==5){
		text((xrange[2]+0.3*(xrange[2]-xrange[1])),150, "Early and late \n dry season \nfires", xpd=NA)
		}

		}}
mtext(side=1, outer=TRUE, "Log short-term population growth rate", line=1, cex=0.8)		
mtext(side=2, outer=TRUE, "Number of simulations out of 1000", line=2, cex=0.8)		

#--------------------------------------------------------------------------------
# Example of 100 year stage dynamics (Fig. 7)
#--------------------------------------------------------------------------------
colstageLight<-c(colors()[c(558, rep(105,4), 432, 432, 432)])
# pdf(file="~/Google Drive/Plant Matrix model/Draft/Review Ecol Mono/Figures/AppendixS4_all_transient%03d.pdf",width=6.3, height=7, onefile=TRUE)
# quartz(width=6.3, height=7)
png(filename="~/Google Drive/Plant Matrix model/Draft/Review Ecol Mono/Aug26/AppendixS4_all_transient%03d.png", width=1258, height=1400, pointsize=34)
for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
		
		par(mfrow=c(4,2), mar=c(3,3,1,0), oma=c(1,1.5,3,1))
		for(k in 1:8){
			plot(1:nt, N_det[[i,j]]$N[k,2:101], "n", bty="l", xlab="", ylab="", las=1, ylim=c(0, max(N_det[[i,j]]$N[k,2:101])*1.5))
			polygon(x=c(1:100, 100:1), y=c(N_all[[i,j]][[2]]$rangeN[1,k,1:100], rev(N_all[[i,j]][[2]]$rangeN[2,k,1:100])), border=NA, col=colstageLight[k])
			lines(1:nt, N_det[[i,j]]$N[k,2:101], lwd=2.5)
			
			mtext(side=3, line=0, adj=0, cex=par('cex'), c(" a) Seedlings", " b) Small juveniles", " c) Large juveniles", " d) Small saplings", " e) Large saplings", " f) Poles", " g) Adults", " h) Large adults")[k], col=colstage[k])		
			}
		mtext(side=1, outer=TRUE, "Time (yrs)", line=0, cex=par("cex"))	
		mtext(side=2, outer=TRUE, "Number of individuals", cex=par("cex"))	
		
		mtext(side=3, outer=TRUE, paste("Fire scenario ", j, " in ", c("sorghum", "non-sorghum")[i], sep=""), line=1)

}}
dev.off()


#--------------------------------------------------------------------------------
# Example of 100 year stable stage distributions 
#--------------------------------------------------------------------------------

pdf(file="~/Google Drive/Plant Matrix model/Draft/Review Ecol Mono/Figures/AppendixS4_all_stagedist%03d.pdf",width=2.5, height=6, onefile=TRUE)
# quartz(width=2.5, height=6)

for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios

	par(mfrow=c(4,1), mar=c(2,2,2,1), oma=c(2,2,3,0))
	for(k in 1:4){
		bp<-barplot(N_all[[i,j]][[2]]$SS[,k], col=colstage, las=1, ylim=c(0,0.7))
		axis(side=1, at=bp, tck=0, labels=c(1:8), lwd=1)
			
			mtext(side=3, adj=0, substitute(paste(L, ") ", S, Y), list(L=letters[k], S=c(0, 30, 60, 100)[k], Y=c(" year", rep(" years", 3))[i])), cex=par("cex"), line=0.5)
				}
	mtext(side=2, "Proportion of individuals", outer=TRUE, cex=par('cex'), line=1)
	mtext(side=1, "Stage", outer=TRUE, cex=par('cex'))
	
	mtext(side=3, outer=TRUE, paste("Fire scenario ", j, " in ", c("sorghum", "non-sorghum")[i], sep=""), line=1, cex=0.8)
}}

dev.off()

# Double-check these results:
i<-1
j<-2

pdf(file="~/Google Drive/Plant Matrix model/Draft/Review Ecol Mono/Figures/Check_stage_dist_20180712.pdf",width=6.3, height=7, onefile=TRUE)

for(i in 1:2){
	for(j in 1:5){
	SS.ij<-matrix(nrow=nt+1, ncol=8)
	for(t in 1:(nt+1)){
		ss<-apply(N_all[[i,j]][[2]]$N[,,t], 2, mean)
		SS.ij[t,]<-ss/sum(ss)
	}
	
	# quartz(width=6.3, height=7)
	par(mfrow=c(4,2), mar=c(3,3,1,0), oma=c(1,1.5,3,1))
	for(k in 1:8){
		plot(1:(nt+1), SS.ij[,k], "l", bty="l", xlab="", ylab="", las=1, ylim=c(0, 0.7), col=colstage[k], lwd=1.5)
		abline(v=c(0,30,60,100), lty=3)
		
		mtext(side=3, line=0, adj=0, cex=par('cex'), c(" a) Seedlings", " b) Small juveniles", " c) Large juveniles", " d) Small saplings", " e) Large saplings", " f) Poles", " g) Adults", " h) Large adults")[k], col=colstage[k])	
	}
	
	mtext(side=1, outer=TRUE, "Time (yrs)", line=0, cex=par("cex"))	
	mtext(side=2, outer=TRUE, "Mean proportion of individuals in stochastic sims", cex=par("cex"))	
	mtext(side=3, outer=TRUE, paste("Fire scenario ", j, " in ", c("sorghum", "non-sorghum")[i], sep=""), line=1)
}}

dev.off()		
#------------------------------------------------------------------------------------------------------------
# Figure 10: stage dist at year 100
#------------------------------------------------------------------------------------------------------------
quartz(width=6.3, height=7, pointsize=14)
stage.t<-c(1,30,60,100)
k<-c(1,3,5,7,9,2,4,6,8,10)
par(mfcol=c(5,2), mar=c(2,1,1,0), oma=c(1,4,2,8))

for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
	
	# if(scenarios$lambda.det[m]>1&scenarios$lambda.stoch.mean[m]>1){
		# densityM=c(NA,NA); angleM<-c(NA,NA)
		# colM<-c(colF[3], paste(colF[3], "50", sep=""))
		# textM<-c(expression(lambda[d]>1), expression(lambda[s]>1))
	# }else if(scenarios$lambda.det[m]<1&scenarios$lambda.stoch.mean[m]>1){
		# densityM=c(50,NA); angleM<-c(45,NA)
		# colM<-c(colF[2], paste(colF[3], "50", sep=""))
		# textM<-c(expression(lambda[d]<1), expression(lambda[s]>1))
	# }else if(scenarios$lambda.det[m]<1&scenarios$lambda.stoch.mean[m]<1){
		# densityM=c(50,50); angleM<-c(45,45)
		# colM<-c(colF[2], paste(colF[2], "50", sep=""))
		# textM<-c(expression(lambda[d]<1), expression(lambda[s]<1))
		# }
		
	bp<-barplot(N_all[[i,j]][[2]]$SS[,4], beside=TRUE, col=colstage, yaxt="n", ylim=c(0, 0.6))
	if(i==1) axis(side=2, at=c(0, 0.2, 0.4, 0.6), las=1) else axis(side=2, at=c(0, 0.2, 0.4, 0.6), labels=FALSE) 
	axis(side=1, at=bp, tck=0, labels=c(1:8), lwd=1)
					
	mtext(side=3, adj=0, paste("  ", letters[k[(i-1)*5+j]], ")", sep=""), cex=par("cex"))
	# text(20, 0.55, textM[1], col=1, cex=1.2)
	# text(20, 0.40, textM[2], col=grey(0.5), cex=1.2)
	
	
	if(j==1) mtext(side=3, c("Sorghum","Non-sorghum")[i], line=2, cex=par('cex')) 
	
	if(i==2&j==1){
		text(max(bp)*1.35, 0.35, "Early dry \nseason fire", xpd=NA)
		text(max(bp)*1.35, 0.2, expression(italic(x)==2), xpd=NA)
		}
		
	if(i==2&j==2){
		text(max(bp)*1.35, 0.35, "Early dry \nseason fire", xpd=NA)
		text(max(bp)*1.35, 0.2, expression(italic(x)==4), xpd=NA)
		}
		
	if(i==2&j==3){
		text(max(bp)*1.35, 0.35, "Late dry \nseason fire", xpd=NA)
		text(max(bp)*1.35, 0.2, expression(italic(x)==3), xpd=NA)
		}
		
	if(i==2&j==4){
		text(max(bp)*1.35, 0.35, "Late dry \nseason fire", xpd=NA)
		text(max(bp)*1.35, 0.2, expression(italic(x)==5), xpd=NA)
		}
		
	if(i==2&j==5){
		text(max(bp)*1.35, 0.35, "Early and late \n dry season \nfire (NENEL)", xpd=NA)
		}
}}
mtext(side=2, "Proportion of individuals", outer=TRUE, cex=par('cex'), line=2)
mtext(side=1, "Stage", outer=TRUE, cex=par('cex'), line=0)


# Proportion of trees in larger stages
propL<-matrix(NA, nrow=dim(scenarios)[1], ncol=2)
colnames(propL)<-c("det", "stoch")
for(m in 1:dim(scenarios)[1]){
	SS.t<-N.all[[m]][['SS']]
	propL[m,]<-round(apply(SS.t[,6:8,4], 1, sum), 2)
}



