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
# Run long-term simulations for 40,000 timesteps
######################################################################################

# Autocorrelation between fire types in stochastic simulations
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
nt<-40000 # Number of time steps
n<-1 # Number of independent simulations

# Burn sequence for deterministic calculations
burn.seq_all=list(c(1,2), c(1, 2, 2, 2), c(1, 2, 2), c(1,2,2,2,2), c(3,1,3,1,2))

# Composite transition matrices for each burn sequence
A._all<-list(
	list(A[[2,1]], A[[1,1]]), list(A[[2,2]], A[[1,2]]),
	list(A[[2,1]], A[[1,1]]), list(A[[2,2]], A[[1,2]]),
	list(A[[3,1]], A[[1,1]]), list(A[[3,2]], A[[1,2]]),
	list(A[[3,1]], A[[1,1]]), list(A[[3,2]], A[[1,2]]),
	list(A[[2,1]], A[[3,1]], A[[1,1]]), list(A[[2,2]], A[[3,2]], A[[1,2]]))
dim(A._all)<-c(2,5)

# Initial starting distribution is steady-state for fire type start.state_all
start.state_all<-c(2,2,2,2,3)

# Create lists to store output of all simulations:
N_all<-list(); length(N_all)<-10; dim(N_all)<-c(2,5)
for(i in 1:2){for(j in 1:5){N_all[[i,j]]<-list(det=NA, stoch.neg=NA, stoch.0=NA, stoch.pos=NA)}}
lambda_all<-N_all

t0<-proc.time()[3] # Time simulatiuon
set.seed(32784) # Change to run in parallel?
for(i in 1:2){ #For sorghum and non-sorghum
	for(j in 1:5){ # For each of 5 different scenarios
		# Deterministic:
		N_all[[i,j]][['det']]<-sim.det(nt=nt, burn.seq=burn.seq_all[[j]], A.=A._all[[i,j]], start.state=start.state_all[j])
		lambda_all[[i,j]][['det']]<-N_all[[i,j]][['det']]$lambda_d
		
		#Stochastic
		for(k in 1:3){
			if(j==5&k==3){}else{
				N_all[[i,j]][[k+1]]<-sim.stoch(env.matrix=env.mat_all[[k,j]], nt=nt, A.=A._all[[i,j]], start.state=start.state_all[j], include.N=TRUE)
				lambda_all[[i,j]][[k+1]]<-N_all[[i,j]][[k+1]]$log.lambda_s
				}
			}#end k
	} # end j scenarios
}# end sog/nonsorg
time.tot<-(proc.time()[3]-t0)/60
cat("Run time: ", time.tot, "mins")

#-------------------------------------------------------
# Table 5

table5<-data.frame(scenario=c(1:10), sorghum=c(rep("Y", 5), rep("N", 5)), lambda_d=rep(NA, 10), lambda_s=rep(NA, 10), Pe1=rep(NA, 10), Pe2=rep(NA, 10), Pe3=rep(NA, 10))
for(i in 1:2){
	for(j in 1:5){
		ind<-(i-1)*5+j
		table5$lambda_d[ind]<-sprintf("%.3f", round(lambda_all[[i,j]]$det, 3))
		table5$lambda_s[ind]<-paste(sprintf("%.3f", round(exp(N_all[[i,j]]$stoch.0$log.lambda_s), 3)), " (", sprintf("%.3f" , round(exp(N_all[[i,j]]$stoch.0$log.lambda_s.CI[1]), 3)), ", ", sprintf("%.3f", round(exp(N_all[[i,j]]$stoch.0$log.lambda_s.CI[2]), 3)), ")", sep="")
		
		for(k in 1:3){
			if(N_all[[i,j]]$stoch.0$P.quasi[k]<0.001) table5[ind,4+k]<-"<0.001" else table5[ind,4+k]<-sprintf("%.3f", round(N_all[[i,j]]$stoch.0$P.quasi[k], 3))}
		
}}

write.csv(table5, "Table5.csv")

#-------------------------------------------------------
# Length of stoch sim
lT<-matrix(rep(NA, 30)); dim(lT)<-c(2,5,3)
for(i in 1:2){
	for(j in 1:5){
		for(k in 1:3){
			if(j==5&k==3){}else lT[i,j,k]<-N_all[[i,j]][[k+1]]$T.max
			}}}
		
i<-1; j<-1
k<-2
X<-N_all[[i,j]][[k]]


LLs<-numeric(nt)
NT<-seq(1, nt, 10)
for(i in 1:nt) LLs[i]<-1/i*sum(log(X$N.tot[i+1]/X$N.tot[i]))
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(c(1:nt)[NT], LLs[NT], "l", xlab="Timesteps (T)", ylab="", bty="l", las=1)
plot(c(1:nt)[NT], LLs[NT], "l", xlab="Timesteps (T)", ylab="", bty="l", las=1, ylim=range(LLs[3000:nt]))
mtext(outer=TRUE, side=2, expression(paste("Stochastic population growth rate (", widehat(log(lambda[s])), ")", sep="")), line=-1)


# Change in log lambda over T Just for no autocorrelation
par(mfcol=c(5,2), mar=c(2,4,1,1), oma=c(2,2,4,0))
thin<-seq(1,40000,10)
X<-list(); length(X)<-10; dim(X)<-c(2,5)
for(i in 1:2){
	for(j in 1:5){
		X[[i,j]]<-numeric(40000)
		for(t in 1:40000) X[[i,j]][t]<-1/t*sum(log(N_all[[i,j]][[3]]$N.tot[2:(t+1)]/N_all[[i,j]][[3]]$N.tot[1:t]))	
}}


for(i in 1:2){
	for(j in 1:5){					
		if(sum(is.na(X[[i,j]])>0)){ind<-min(which(is.na(X[[i,j]])==TRUE|X[[i,j]]==Inf))-1}else ind<-40000
				 
		plot(c(1:40000)[thin],  X[[i,j]][thin], "l", xlab="", ylab="", bty="l", las=1, col=colF[2], ylim=range(X[[i,j]][1000:ind]))
		abline(h=X[[i,j]][1000], lty=3, lwd=2)
		abline(h=X[[i,j]][10000], lty=2, lwd=1.5)
		abline(h=X[[i,j]][ind])
		
		u<-par('usr')
		if(i==1&j==1) legend(u[1], u[4]+0.8*(u[4]-u[3]), lty=c(3,2,1), lwd=c(2,1.5,1), title=expression(paste(widehat(log(lambda[s])), "estimated at")), c("T=1000", "T=10,000", "T=40,000 (or max)"), xpd=NA, bty="n")
		
}}
mtext(side=1, outer=TRUE, "Timsteps (T)", line=1)
mtext(outer=TRUE, side=2, expression(paste("Stochastic population growth rate (", widehat(log(lambda[s])), ")", sep="")), line=0)
mtext(side=3, outer=TRUE, "Calculated from entire series 1:T", line=1)


#-------------------------------------------------------------------------------------------------------------------
# SLiding window method for log lambda
par(mfcol=c(5,2), mar=c(2,4,1,1), oma=c(2,2,4,0))
thin<-seq(1,40000,10)
for(i in 1:2){
	for(j in 1:5){
		
		thin<-seq(1000,40000,10)
		Y<-numeric(39000)
		for(t in 1:39000) Y[t]<-1/1000*sum(log(N_all[[i,j]][[3]]$N.tot[(t+1):(t+1001)]/N_all[[i,j]][[3]]$N.tot[t:(t+1000)]))
			
		if(sum(is.na(Y)>0)){ind<-min(which(is.na(Y)==TRUE|Y==Inf))-1}else ind<-40000
				 
		plot(c(1000:40000)[thin], Y[thin], "l", xlab="", ylab="", bty="l", las=1, col=colF[2])
		abline(h=Y[1000], lty=3, lwd=2)
		abline(h=Y[10000], lty=2, lwd=1.5)
		abline(h=Y[ind-1000])

		u<-par('usr')
		if(i==1&j==1) legend(u[1], u[4]+0.8*(u[4]-u[3]), lty=c(3,2,1), lwd=c(2,1.5,1), title=expression(paste(widehat(log(lambda[s])), "estimated at")), c("T=1000", "T=10,000", "T=40,000 (or max)"), xpd=NA, bty="n")
		
}}
mtext(side=1, outer=TRUE, "Timsteps (T)", line=1)
mtext(side=3, outer=TRUE, "Sliding window from (T-1000):T", line=1)
mtext(outer=TRUE, side=2, expression(paste("Stochastic population growth rate (", widehat(log(lambda[s])), ")", sep="")), line=0)


#-------------------------------------------------------
# Appendix: effect of autocorrelation
par(mfcol=c(5,2), mar=c(2,4,1,1), oma=c(2,2,4,0))
for(i in 1:2){
	for(j in 1:5){
		if(j==5){
			pts<-c(N_all[[i,j]]$stoch.neg$log.lambda_s, N_all[[i,j]]$stoch.0$log.lambda_s)
			li<-c(N_all[[i,j]]$stoch.neg$log.lambda_s.CI[1], N_all[[i,j]]$stoch.0$log.lambda_s.CI[1])
			ui<-c(N_all[[i,j]]$stoch.neg$log.lambda_s.CI[2], N_all[[i,j]]$stoch.0$log.lambda_s.CI[2])
			
		}else{
			pts<-c(N_all[[i,j]]$stoch.neg$log.lambda_s, N_all[[i,j]]$stoch.0$log.lambda_s, N_all[[i,j]]$stoch.pos$log.lambda_s)
			li<-c(N_all[[i,j]]$stoch.neg$log.lambda_s.CI[1], N_all[[i,j]]$stoch.0$log.lambda_s.CI[1], N_all[[i,j]]$stoch.pos$log.lambda_s.CI[1])
			ui<-c(N_all[[i,j]]$stoch.neg$log.lambda_s.CI[2], N_all[[i,j]]$stoch.0$log.lambda_s.CI[2], N_all[[i,j]]$stoch.pos$log.lambda_s.CI[2])
			}
		
		yaxp.all<-list(c(0, 0.01, 0.02), c(0.03, 0.035, 0.04), c(-0.02, -0.01, 0), c(0.004, 0.010, 0.016), c(-0.025, -0.020, -0.015), c(0.006, 0.010, 0.014), c(-0.001, 0.000, 0.001, 0.002), c(-0.02, -0.018, -0.016), c(-0.0145, -0.0140, -0.0135), c(-0.012, -0.008, -0.004))
		plotCI(1:length(pts), pts, li=li, ui=ui, xlim=c(0,4), xaxt="n", bty="l", xlab="", gap=0, ylim=range(c(li, ui, log(lambda_all[[i,j]]$det), yaxp.all[[(i-1)*5+j]])), pch=19, ylab="", col=colF, yaxt="n", cex=1.2)
		axis(side=2, at=yaxp.all[[(i-1)*5+j]], las=1)
		abline(h=log(lambda_all[[i,j]]$det), lty=2)
		mtext(side=3, line=0.5, adj=0, paste(letters[(i-1)*5+j], ")", sep=""))
		if(j==5){
			axis(side=1, at=1, labels="Neg", col=colF[1], col.axis=colF[1])
			axis(side=1, at=2, labels="Zero", col=colF[2], col.axis=colF[2])
			axis(side=1, at=3, labels="Pos", col=colF[3], col.axis=colF[3])
			}
			
		abline(h=0)
		# if(i==1) mtext(side=2, line=3.5, paste("Scenario ", j), cex=1.2)
		if(j==1) mtext(side=3, line=3, c("Sorghum", "Non-sorghum")[i])
	
		u<-par('usr')
		# if(j==1&i==1) legend(u[2], u[4]*2.2, col=c(colF,1), lwd=1, pch=c(rep(19,3), NA), lty=c(1,1,1,3), c("Neg. auto.", "Zero auto", "Pos. auto", "Deterministic"), bty="n", xpd=NA)
		}}

mtext(side=1, outer=TRUE, "Autocorrelation in environment", line=0.5)
mtext(side=2, outer=TRUE, expression(paste("Stochastic population growth rate (", widehat(lambda[s]), ")", sep="")))


#-------------------------------------------------------
# Probability of quasi-extinction for each level of autocorrelation
par(mfcol=c(5,2), mar=c(2,4,1,1), oma=c(2,2,4,0))
for(i in 1:2){
	for(j in 1:5){
		if(j==5) Pe<-list(N_all[[i,j]]$stoch.neg$P.quasi.long, N_all[[i,j]]$stoch.0$P.quasi.long) else Pe<-list(N_all[[i,j]]$stoch.neg$P.quasi.long, N_all[[i,j]]$stoch.0$P.quasi.long, N_all[[i,j]]$stoch.pos$P.quasi.long)				
		
		plot(seq(0.01, 0.99, 0.01), Pe[[1]], "l", col=colF[1], bty="l", xlab="", ylim=c(0,1), ylab="", las=1, lwd=2, xaxt="n")
		lines(seq(0.01, 0.99, 0.01), Pe[[2]], col=colF[2])
		if(j!=5) lines(seq(0.01, 0.99, 0.01), Pe[[3]], col=colF[3], lty=2)
		
		mtext(side=3, line=0.5, adj=0, paste(letters[(i-1)*5+j], ")", sep=""))
		if(j==5) axis(side=1) else axis(side=1, labels=FALSE)
			
		if(j==1) mtext(side=3, line=3, c("Sorghum", "Non-sorghum")[i])
	
		u<-par('usr')
		if(j==1&i==1) legend(0.15, 1.1, col=c(colF,1), lwd=c(2,1,1), lty=c(1,1,2), c("Neg. auto.", "Zero auto", "Pos. auto"), bty="n", xpd=NA, title="Autocorrelation in\nthe environment")
		}}

mtext(side=1, outer=TRUE, expression(paste("Proportion of initial population size (", theta, ")")), line=1)
mtext(side=2, outer=TRUE, expression(paste("Probability of quasi-extinction ", P[italic(e)](theta), sep="")))


