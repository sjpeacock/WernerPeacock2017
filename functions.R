cols<-colorRampPalette(rev(c("#0000FF80", "#228B2280", "#FFD70080", "#FF000080")))
	
	
stableStage<-function(A){
	Re(eigen(A)$vector[,1]/sum(eigen(A)$vector[,1])) # stable stage
}


plot.SS<-function(A, returnSS=FALSE){
	lambda<-matrix(NA, nrow=4, ncol=2)
	for(i in 1:4){
		for(j in 1:2){
			lambda[i,j]<-calcLambda(A[[i,j]])
		}
	}
	
	SS<-list(); length(SS)<-8; dim(SS)<-c(4,2)
for(i in 1:4){
	for(j in 1:2){
		SS[[i,j]]<-stableStage(A[[i,j]])
		}
	}
	quartz(width=5, height=6)
	par(mfrow=c(4,2), mar=c(1,1,1,0), oma=c(4.5,3,3,1))
	for(i in 1:4){
		for(j in 1:2){
			bp<-barplot(SS[[i,j]], yaxt="n", ylim=c(0,1), space=0, border="white", col=cols(8))
			if(j==1) axis(side=2, las=1) else axis(side=2, labels=FALSE)
			abline(h=0)
			#text(bp[1], 0.75, substitute(paste(lambda, "=", L), list(L=round(lambda[i,j],3))), adj=0)
			#text(bp[1], 0.95, paste(c("no", "early", "late", "wet")[i], " fire \n", c("sorghum", "no sorghum")[j], sep=""), adj=0, xpd=NA)	
			if(i==4) axis(side=1, at=c(bp, 8.5)-0.5, labels=FALSE)
			if(i==4) text(bp+0.25, -0.05, c("seedling", "sm. juv.", "juv.", "sm. sap", "sap.", "pole", "adult", "lg. adult"), xpd=NA, col=cols(8), font=2, srt=90, pos=2, adj=c(0.5, 0))
			mtext(side=3, adj=0, substitute(paste(A, ") ", FT, " fire, ", ST, ", ", lambda, "=", L), list(A=letters[(i-1)*2+j], L=sprintf("%.3f", round(lambda[i,j],3)), FT = c("no", "early", "late", "wet")[i], ST=c("sorghum", "no sorghum")[j])) , cex=par("cex"))	
		}
	}
	mtext(side=1, outer=TRUE, "Stage", line=3, cex=0.8)
	mtext(side=2, outer=TRUE, "Proportion of individuals in stage", line=1.5, cex=0.8)
	mtext(side=3, outer=TRUE, "Stable stage distributions", font=2, line=1)
	
	if(returnSS==TRUE) return(list(lambda, SS))
}
calcLambda<-function(A){
	Re(eigen(A)$value[1])
}

simPop<-function(A, init, tmax){
	N<-matrix(NA, nrow=length(init), ncol=tmax+1)
	N[,1]<-init
	
	for(i in 1:tmax){
		N[,i+1]<-A%*%N[,i]
	}
}

ImportMatrices<-function(){
	Amat<-read.csv("matrices.csv", na.string="")
	Amat[,11]<-as.numeric(rep(NA, 44))
	Amat[which(Amat[,c(3,12)]=="sdlng"|Amat[,c(3,12)]=="From"),]<-NA
	
	A<-as.matrix(Amat[c(3:10,14:21,25:32,36:43),3:19])
	A[which(is.na(A)==TRUE, arr.ind=TRUE)]<-0
	A.num<-matrix(nrow=dim(A)[1], ncol=dim(A)[2])
	for(i in 1:dim(A)[1]) A.num[i,]<-as.numeric(A[i,])
	
	A.all<-list(); length(A.all)<-8; dim(A.all)<-c(4,2)
	
	for(i in 1:4){
		first<-8*(i-1)+1
		A.all[[i,1]]<-A.num[first:(first+7),1:8]
		A.all[[i,2]]<-A.num[first:(first+7),10:17]
		
	}
	
	for(i in 1:4){
		for(j in 1:2){
			A.all[[i,j]]<-as.matrix(A.all[[i,j]])
			rownames(A.all[[i,j]])<-1:8
			colnames(A.all[[i,j]])<-1:8
		}
	}
	return(A.all)

	}
	
sensitivity <- function(A) {
  d <- eigen(A)$values   # eigen values
  w <- eigen(A)$vectors  # right eigen vectors
  v <- Conj(solve(w))    # complex conjugate of left eigen vectors
  # output of eigenvalues is decreasingly sorted.
  Re(v[1,] %*% t(w[,1]))
}
	
elasticity<-function(A){
	s<-sensitivity(A)	
	s*A/calcLambda(A)
}	
	
	
# Sim scenario
# function that takes the burn sequence and sorghum/non-sorghum and returns N
# stochastic and deterministic

sim.N<-function(sorghum=Y, burn.sequence="NE", stoch.probs=NA, nt=100, n=1000, A.=A, init=NA, stage.t=c(1,30,60,100), make.smaller=FALSE){

	#====
	s<-as.numeric(sorghum=="N")+1
	burn.seq<-strsplit(burn.sequence, split="")[[1]]
	burn.seq.num<-as.numeric(factor(burn.seq, levels=c("N", "E", "L")))
	burn.seq.length<-length(burn.seq.num)
	
	if(is.na(stoch.probs)==TRUE){
		if(length(unique(burn.seq.num))<3) stoch.probs<-length(which(burn.seq.num!=1))/burn.seq.length else if(length(unique(burn.seq.num))==3){
			stoch.probs<-c(length(which(burn.seq.num==1))/burn.seq.length, length(which(burn.seq.num==2))/burn.seq.length, length(which(burn.seq.num==3))/burn.seq.length)
			stoch.probs<-rbind(stoch.probs, stoch.probs)
			}
		}
		
	
	# Stochastic series of burns
	burned<-list(
		deterministic=rep(burn.seq.num, ceiling(nt/burn.seq.length)), 
		stochastic=matrix(NA, nrow=nt, ncol=n))
	
	if(length(stoch.probs)>1){ #If there three fire types
		burned[[2]][1,]<-burn.seq.num[1]
		for(i in 1:n){
			for(t in 2:nt){
				if(burned[[2]][t-1,i]==2){
					burned[[2]][t,i]<-which(rmultinom(1, size=1, prob=stoch.probs[2,])==1)
				}else{
					burned[[2]][t,i]<-which(rmultinom(1, size=1, prob=stoch.probs[1,])==1) 
					}
			}
		}
	} else{
		burned[[2]]<-matrix(rbinom(nt*n, size=1, prob=stoch.probs), nrow=nt, ncol=n)*(max(burn.seq.num)-1)
		burned[[2]]<-burned[[2]]+1
		}
	
	
	# Matrix of 8 stage, sorg/non sorg (2), deterministic/stochastic (2), timesteps (nt)
	N<-list(); length(N)<-2
	N[[1]]<-matrix(); length(N[[1]])<-8*(nt+1); dim(N[[1]])<-c(8,(nt+1)) #Deterministic
	N[[2]]<-matrix(); length(N[[2]])<-n*8*(nt+1); dim(N[[2]])<-c(n,8,(nt+1)) #Stochastic
	
	# Initial condition: Starting distirbution is the stable stage for no fire
	if(is.na(init[1])==TRUE) init<-1000*stableStage(A.[[1,s]])
	N[[1]][,1]<-init
	for(i in 1:n){
		N[[2]][i,,1]<-init
		}
		
	# Simulation
	for(t in 1:nt){ # for each timestep
		
		# Deterministic
		N[[1]][,t+1]<-A.[[burned[[1]][t],s]]%*%N[[1]][,t]
		
		# Stochastic
		for(i in 1:n){
			N[[2]][i,,t+1]<-A.[[burned[[2]][t,i],s]]%*%N[[2]][i,,t]
			} #end i
		
	}
	
	
	# range in stochastic
	rangeN<-matrix(rep(NA, 2*8*nt)); dim(rangeN)<-c(2, 8, nt)
	for(k in 1:8){
		for(t in 1:nt){
			rangeN[,k,t]<-quantile(N[[2]][,k,t+1], c(0.025, 0.975))
		}
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
	
	
	A1<-A.[[burn.seq.num[1],s]]
	for(i in 1:(length(burn.seq)-1))A1<-A.[[burn.seq.num[i+1],s]]%*%A1
	
	lambda.all<-c(det=round(calcLambda(A1)^(1/burn.seq.length), 3),
		stoch.mean=round(exp(mean(log.lambda_s)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.025)), 3),
		stoch=round(exp(quantile(log.lambda_s, 0.975)), 3)
		)
	
		
	#--------------------------------------------------
	# Stable stage
	
	SS.t<-matrix(rep(NA, 2*8*length(stage.t))); dim(SS.t)<-c(2, 8, length(stage.t))
	for(t in 1:length(stage.t)){
			Ntd<-N[[1]][,stage.t[t]]
			SS.t[1,,t]<-Ntd/sum(Ntd)
			
			Nts<-apply(N[[2]][,,stage.t[t]], 2, mean)
			SS.t[2,,t]<-Nts/sum(Nts)
			
		}
	
	
	if(make.smaller==TRUE) N[[2]]<-NA
	return(list(N=N, lambda=lambda.all, SS=SS.t, rangeN=rangeN))
	
}
