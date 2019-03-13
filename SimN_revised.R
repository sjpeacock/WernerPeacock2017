# Sim scenario
# function that takes the burn sequence and sorghum/non-sorghum and returns N
# stochastic and deterministic

sim.stoch<-function(env.matrix, nt=40000, A.=list(A[[2,1]], A[[1,1]]), start.state=2, include.N=FALSE){

# 1) Create a sequence for environment
	burn.seq<-numeric(nt)
	burn.seq[1]<-start.state # All start in unburned state (no fire)
	
	n.states<-dim(env.matrix)[1] #Number of environment states included
	if(n.states!=length(A.)) cat("Error: dimension of env.matrix != number of trans. matrices provided in A.")
	
	for(t in 2:nt) burn.seq[t]<-sample(c(1:n.states), size=1, replace=FALSE, prob=env.matrix[,burn.seq[t-1]])
	
# 2) Simulate population dynamics given burn sequence
	# Matrix of 8 stage, timesteps (nt+1)
	N<-matrix(rep(NA, 8*(nt+1)), nrow=8, ncol=(nt+1)) #Stochastic
	
	# Initial condition: Starting distirbution is the stable stage for no fire
	N[,1]<-1000*stableStage(A.[[start.state]])
		
	# Simulation
	for(t in 1:nt) N[,t+1]<-A.[[burn.seq[t]]]%*%N[,t]
	
	# Sum stages to get total population size at each timestep
	N.tot<-apply(N, 2, sum)
	
	if(sum(is.na(N.tot))>0){
		T.max<-min(which(is.na(N.tot)==TRUE|N.tot==Inf))-2 
		warning(paste("N=NaN or Inf, T.max=", T.max))
	}else{T.max<-nt}
	
	
# 3) Calculate log lambda_s at increasing timestep T=nt from eq'n (25) in Caswell and Kaye (2001)
	# Stochastic lambda
	log.lambda_s<-1/T.max*sum(log(N.tot[2:(T.max+1)]/N.tot[1:T.max]))	
			
	# 95% CI on stochastic lambda from eq'n (26) in Caswell and Kaye (2001)
	V.log<-var(log(N.tot[2:(T.max+1)]/N.tot[1:T.max]))/T.max
	log.lambda_s.CI<-c(log.lambda_s-1.96*sqrt(V.log), log.lambda_s+1.96*sqrt(V.log))
	
# 4) Calculate Prob quasi-extinction from eq'n (30) in Caswell and Kaye (2001)
	n.stages<-dim(A.[[1]])[1]
	# Calculate average matrix from eq'n (29)
	if(dim(env.matrix)[1]==2){
		A.big<-cbind(
			rbind(A.[[1]], matrix(rep(0, n.stages*n.stages), nrow=n.stages)),
			rbind(matrix(rep(0, n.stages*n.stages), nrow=n.stages), A.[[2]]))
		}else if(dim(env.matrix)[1]==3){
		A.big<-cbind(
			rbind(A.[[1]], matrix(rep(0, n.stages*2*n.stages), nrow=n.stages*2)),
			rbind(matrix(rep(0, n.stages*n.stages), nrow=n.stages), A.[[2]], matrix(rep(0, n.stages*n.stages), nrow=n.stages)),
			rbind(matrix(rep(0, n.stages*2*n.stages), nrow=n.stages*2), A.[[3]]))			
			}
		
	A.avg<-A.big%*%(env.matrix%x%diag(n.stages))
	
	# Calculate log mu
	log.mu<-as.numeric(log(eigen(A.avg)$values[1]))

	# Calculate sigma^2
	sigma2<-2*(log.mu-log.lambda_s)

	# Threshold for extinction (% of original population size; here I use three values)
	theta<-c(0.1, 0.5, 0.9)
	if(log.lambda_s<=0) P.quasi<-rep(1, length(theta)) else P.quasi<-exp(2*log.lambda_s*log(theta)/sigma2)

	theta.long<-seq(0.01, 0.99, 0.01)
	if(log.lambda_s<=0) P.quasi.long<-rep(1, length(theta.long)) else P.quasi.long<-exp(2*log.lambda_s*log(theta.long)/sigma2)


	# "Average" Stage distribution
	ASD<-numeric(n.stages)
	for(i in 1:n.stages) ASD[i]<-mean(N[i,]/N.tot)


	if(include.N==TRUE){return(list(log.lambda_s=log.lambda_s, log.lambda_s.CI= log.lambda_s.CI, P.quasi=P.quasi, theta=theta,T.max=T.max, N=N, N.tot=N.tot, stage.dist=ASD))}else{return(list(log.lambda_s=log.lambda_s, log.lambda_s.CI= log.lambda_s.CI, P.quasi=P.quasi, theta=theta, T.max=T.max, P.quasi.long=P.quasi.long, stage.dist=ASD))}

}


#-----------------------------------------------------------------
# for n>1
sim.stoch_n<-function(env.matrix, nt=1000, n=100, A.=list(A[[2,1]], A[[1,1]]), start.state=2, transient=FALSE, sorghum=TRUE){

	stage.t <- c(0,30,60,100)
	
	burn.seq<-matrix(NA, nrow=n, ncol=nt)
	if(length(start.state)>1){
		if(sorghum==TRUE) burn.seq[,1]<-2 else burn.seq[,1]<-1
	}else{
		burn.seq[,1]<-start.state 
		}
	
	n.states<-dim(env.matrix)[1] #Number of environment states included
	if(n.states!=length(A.)) cat("Error: dimension of env.matrix != number of trans. matrices provided in A.")
	
	for(i in 1:n){
		for(t in 2:nt){
			burn.seq[i,t]<-sample(c(1:n.states), size=1, replace=FALSE, prob=env.matrix[,burn.seq[i,t-1]])
		}}
	
	#====
	# Matrix of 1000 sims, 8 stage, timesteps (nt+1)
	N<-matrix(); length(N)<-n*8*(nt+1); dim(N)<-c(n,8,(nt+1)) #Stochastic
	
	# Initial condition: Starting distirbution is the stable stage for no fire
	if(length(start.state)==1){
		for(i in 1:n) N[i,,1]<-1000*stableStage(A.[[start.state]])
		}else if(length(start.state==8)){
			for(i in 1:n) N[i,,1]<-start.state
		}
		
	# Simulation
	for(i in 1:n){ # for each simulation
		for(t in 1:nt){ # for each timestep
			N[i,,t+1]<-A.[[burn.seq[i,t]]]%*%N[i,,t]
			}}
	
	N.tot<-matrix(NA, nrow=n, ncol=(nt+1))
	for(i in 1:n) N.tot[i,]<-apply(N[i,,], 2, sum)
	
	# # range in stochastic
	rangeN<-matrix(rep(NA, 2*8*nt)); dim(rangeN)<-c(2,8, nt)
	for(k in 1:8){
		for(t in 1:nt){
			rangeN[,k,t]<-quantile(N[,k,t+1], c(0.025, 0.975))
		}
	}
	
	# Stochastic lambda
	log.lambda_all<-matrix(NA, nrow=n, ncol=nt)
	for(i in 1:n){
		for(t in 1:nt){	
			log.lambda_all[i,t]<-1/t*sum(log(N.tot[i,2:(t+1)])-log(N.tot[i,1:t]))
		}}
	
	log.lambda_s<-log.lambda_all[,nt]
	
	# hist(log.lambda_s[i,], breaks=seq(min(log.lambda_s[i,])-0.001, max(log.lambda_s[i,])+0.001, 0.001), xlim=c(-0.01, 0.03))
	# abline(v=mean(log.lambda_s[i,]), lwd=2, col=2)
	
	# plot(1:nt, log.lambda_s[1,], "l", col="#00000030", lwd=0.5, ylim=range(log.lambda_s))
	# for(i in 2:nt) lines(1:nt, log.lambda_s[i,], "l", col="#00000030", lwd=0.5)
	
	# hist(log.lambda_s[,900:1000])
	# abline(v=mean(log.lambda_s[,900:1000]), col=2, lwd=2)
	
		
	#--------------------------------------------------
	# Stable stage
	
	SS.t<-matrix(rep(NA, 8*length(stage.t))); dim(SS.t)<-c(8, length(stage.t))
	for(i in 1:length(stage.t)){
			Nts<-apply(N[,,stage.t[i]+1], 2, mean)
			SS.t[,i]<-Nts/sum(Nts)	
		}
	
	
	return(list(N=N, N.tot=N.tot, log.lambda_s=log.lambda_s, log.lambda_all=log.lambda_all, SS=SS.t, rangeN=rangeN, burn.seq=burn.seq))
	
}

###########################################################################################################################
sim.det<-function(nt=1000, burn.seq=c(1,2,2), A.=list(A[[2,1]], A[[1,1]]), start.state=2){

	if(length(unique(burn.seq))!=length(A.)) cat("Error: num of env states in burn.seq != number of trans. matrices provided in A.")
	
	# What is the average transition matrix?
	A1<-A.[[burn.seq[1]]]
	for(i in 1:(length(burn.seq)-1)) A1<-A.[[burn.seq[i+1]]]%*%A1
	
	# Calculate lambda from that
	lambda <- calcLambda(A1)^(1/length(burn.seq))
	
	stage.t <- c(1,30,60,100)
	
	burn.seq<-rep(burn.seq, length.out=nt)
	
	#====
	# Matrix of 8 stage, timesteps (nt+1)
	N<-matrix(); length(N)<-8*(nt+1); dim(N)<-c(8,(nt+1)) #Deterministic
	
	# Initial condition: Starting distirbution as number or actual distribution
	if(length(start.state)==8) N[,1]<-start.state else if(length(start.state)==1) N[,1]<-1000*stableStage(A.[[start.state]])  	
		
	# Simulation
	for(t in 1:nt){ # for each timestep
		N[,t+1]<-A.[[burn.seq[t]]]%*%N[,t]
		}
	
	N.tot<-apply(N, 2, sum)
	
	#--------------------------------------------------
	# Stable stage
	
	SS.t<-matrix(rep(NA, 8*length(stage.t))); dim(SS.t)<-c(8, length(stage.t))
	for(i in 1:length(stage.t)){
			Ntd<-N[,stage.t[i]]
			SS.t[,i]<-Ntd/sum(Ntd)			
		}
	
	
	return(list(N=N, lambda_d=lambda, SS=SS.t, burn.seq=burn.seq))
	
}
