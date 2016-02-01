#####################################################################
# IPM SCRIPT FOR PRETTY FIGURES WITH REGIONAL SUMMARIES

# this script should be run right after
# "04.1_basicIPMdemograpy.R"

#####################################################################


## Check parameters	
Regions <- c("Central", "North", "South")
pars <-c("m.parC", "m.parN", "m.parS")
m.parC; m.parN; m.parS
rm(m.par); rm(P); rm(K); rm(G); rm(S)
rm(ProbRepVec); rm(RecMatrix)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1D - Construct kernels and define P, F & K pdfs
## ~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KERNEL PROPERTIES - shared for all regions
		# number of classes, or points for midpoint rule approximation
		n.size = 100
		# tolerance for iterations
		tol = 1.e-8
		# minimum and maximum sizes (0.9*min & 1.1*max size from data)
		minsize = 0.9*(min(rbind(dfclean$z1, dfclean$z), na.rm=TRUE)) # 0.9 to go slightly lower than smallest size
		maxsize = 1.1*(max(rbind(dfclean$z1, dfclean$z), na.rm=TRUE)) # 1.1 to go slightly above largest size
				
		# calculate h - the width of the size classes in the matrix (discretized kernel)
		L= minsize; U= maxsize; n <- n.size
		# b - boundary points, limits of the size bins defining the kernels 
		b = L+c(0:n)*(U-L)/n 
		# y - mesh points, the centers of these bins, and the points at which the kernel is evaluated, under the midpoint rule of numerical integration
		y = 0.5*(b[1:n]+b[2:(n+1)])
		# h - step size, width of the bins 
		h = y[2]-y[1]	

	#--------------------
	# kernel components
	#--------------------
		par(mfrow=c(1,3))
		## runs Grow_z1z for each combination of possible sizes
		G_N = h*outer(y,y,Grow_z1z,m.par=m.parN); dim(G); str(G); class(G); image.plot(G_N, main="N")
		G_C = h*outer(y,y,Grow_z1z,m.par=m.parC); dim(G); str(G); class(G); image.plot(G_C, main="C")
		G_S = h*outer(y,y,Grow_z1z,m.par=m.parS); dim(G); str(G); class(G); image.plot(G_S, main="S")

		S_N = SurvFrom_z(y,m.par=m.parN) # survival vector 
		S_C = SurvFrom_z(y,m.par=m.parC) # survival vector 
		S_S = SurvFrom_z(y,m.par=m.parS) # survival vector 

		ProbRepVec_N= Prob_repr(y,m.par=m.parN) # vector of probability of flowering 
		ProbRepVec_C= Prob_repr(y,m.par=m.parC) # vector of probability of flowering 
		ProbRepVec_S= Prob_repr(y,m.par=m.parS) # vector of probability of flowering 
		
		RecMatrix_N=h*outer(y,y,RecrDist,m.par=m.parN) # reproduction matrix
		RecMatrix_C=h*outer(y,y,RecrDist,m.par=m.parC) # reproduction matrix
		RecMatrix_S=h*outer(y,y,RecrDist,m.par=m.parS) # reproduction matrix

		image.plot(RecMatrix_N, main="N"); image.plot(RecMatrix_C, main="C"); image.plot(RecMatrix_S, main="S")

	# Growth matrix 
		P_N=G_N						     # placeholder; redefine P on the next line
		P_C=G_C						     # placeholder; redefine P on the next line
		P_S=G_S						     # placeholder; redefine P on the next line

		for(i in 1:n) P_N[,i]=G_N[,i]*S_N[i]	 # growth/survival matrix
		for(i in 1:n) P_C[,i]=G_C[,i]*S_C[i]	 # growth/survival matrix
		for(i in 1:n) P_S[,i]=G_S[,i]*S_S[i]	 # growth/survival matrix

		# Fecundity matrix 
		F_N=RecMatrix_N
		F_C=RecMatrix_C
		F_S=RecMatrix_S

		for(i in 1:n) F_N[,i]=RecMatrix_N[,i]*ProbRepVec_N[i]	 # flower/fecundity matrix
		for(i in 1:n) F_C[,i]=RecMatrix_C[,i]*ProbRepVec_C[i]	 # flower/fecundity matrix
		for(i in 1:n) F_S[,i]=RecMatrix_S[,i]*ProbRepVec_S[i]	 # flower/fecundity matrix

		K_N=P_N+F_N # full matrix, add matrices together 
		K_C=P_C+F_C # full matrix, add matrices together 
		K_S=P_S+F_S # full matrix, add matrices together 
		
		
		image.plot(K_N, main="N"); image.plot(K_C, main="C"); image.plot(K_S, main="S")
	
setwd(path.fig)
pdf(file="04.1_IPM_NCS2.pdf", width=11, height=8.5)
	
	# basic lambda estimate	
		lambdas1 <- c((lam <- Re(eigen(K_N)$values[1])), #  0.939
		(lam <- Re(eigen(K_C)$values[1])), #  0.939
		(lam <- Re(eigen(K_S)$values[1]))) #  0.939
	par(mfrow=c(1,1))
		barplot(lambdas1, col=c("darkblue", "darkgreen", "darkred"),
			names.arg=c("North", "Center", "South"),
			ylab="lambda", main="Population Growth Rate Across Regions")
	
##################################################		
	ipmmats <- c("K_N","K_C","K_S")
	letter <- c("N", "C", "S")
	
	# Make elasticity and sensitivty matrix 
	for(i in 1:length(ipmmats)){
		K <- get(ipmmats[i])
		w.eigen <- Re(eigen(K)$vectors[,1])
		stable.dist <- w.eigen/sum(w.eigen)
		assign(paste("stable.dist_",letter[i],sep=""), stable.dist)
		v.eigen <- Re(eigen(t(K))$vectors[,1])
		repro.val <- v.eigen/v.eigen[1]
		assign(paste("repro.val_",letter[i],sep=""), repro.val)
		v.dot.w=sum(stable.dist*repro.val)*h
		sens=outer(repro.val,stable.dist)/v.dot.w
		assign(paste("sens_",letter[i],sep=""), sens)
		elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
		assign(paste("elas_",letter[i],sep=""), elas)
		rm(K); rm(elas); rm(sens); rm(repro.val); rm(stable.dist)
		}
		
		
	# plot out the IPM components 
		par(mfrow=c(1,3),mar=c(4,5,2,2))
		image.plot(y,y,t(K_N), xlab="Size (t)",ylab="Size (t+1)", main="N", col=topo.colors(500))
			contour(y,y,t(K_N), add = TRUE, drawlabels = TRUE)
		image.plot(y,y,t(K_C), xlab="Size (t)",ylab="Size (t+1)", main="C", col=topo.colors(500))
			contour(y,y,t(K_C), add = TRUE, drawlabels = TRUE)
		image.plot(y,y,t(K_S), xlab="Size (t)",ylab="Size (t+1)", main="S", col=topo.colors(500))
			contour(y,y,t(K_S), add = TRUE, drawlabels = TRUE)
				
		plot(y,stable.dist_N,xlab="Size",type="l",main="Stable Size Dist N", ylab=c(0,0.05))
		plot(y,stable.dist_C,xlab="Size",type="l",main="Stable Size Dist C", ylab=c(0,0.05))
		plot(y,stable.dist_S,xlab="Size",type="l",main="Stable Size Dist S", ylab=c(0,0.05))
				
		plot(y,repro.val_N,xlab="Size",type="l",main="RVs N", ylab=c(0,2000))
		plot(y,repro.val_C,xlab="Size",type="l",main="RVs C", ylab=c(0,2000))
		plot(y,repro.val_S,xlab="Size",type="l",main="RVs S", ylab=c(0,2000))

		image.plot(y,y,t(elas_N),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity N")
		image.plot(y,y,t(elas_C),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity C")
		image.plot(y,y,t(elas_S),xlab="Size (t)",ylab="Size (t+1)", main="Elasticity S")

		image.plot(y,y,t(sens_N),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity N")
		image.plot(y,y,t(sens_C),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity C")
		image.plot(y,y,t(sens_S),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity S")

#########################################
# FIX EVICTION FROM MODELS 
		
		# fix eviction for offspring
			for(i in 1:(n/2)) {
			G_N[1,i]<-G_N[1,i]+1-sum(G_N[,i])
			P_N[,i]<-G_N[,i]*S_N[i] # growth/survival matrix
			}
			for(i in 1:(n/2)) {
			G_C[1,i]<-G_C[1,i]+1-sum(G_C[,i])
			P_C[,i]<-G_C[,i]*S_C[i] # growth/survival matrix
			}
			for(i in 1:(n/2)) {
			G_S[1,i]<-G_S[1,i]+1-sum(G_S[,i])
			P_S[,i]<-G_S[,i]*S_S[i] # growth/survival matrix
			}
		

		# fix eviction of large adults
			for(i in (n/2+1):n) {
			G_N[n,i]<-G_N[n,i]+1-sum(G_N[,i])
			P_N[,i]<-G_N[,i]*S_N[i]
			}				
			for(i in (n/2+1):n) {
			G_C[n,i]<-G_C[n,i]+1-sum(G_C[,i])
			P_C[,i]<-G_C[,i]*S_C[i]
			}	
			for(i in (n/2+1):n) {
			G_S[n,i]<-G_S[n,i]+1-sum(G_S[,i])
			P_S[,i]<-G_S[,i]*S_S[i]
			}		
		
		for(i in 1:n) F_N[,i]=RecMatrix_N[,i]*ProbRepVec_N[i]	 # flower/fecundity matrix	
		K_N=P_N+F_N # full matrix
		for(i in 1:n) F_C[,i]=RecMatrix_C[,i]*ProbRepVec_C[i]	 # flower/fecundity matrix	
		K_C=P_C+F_C 
		for(i in 1:n) F_S[,i]=RecMatrix_S[,i]*ProbRepVec_S[i]	 # flower/fecundity matrix	
		K_S=P_S+F_S 		
		
		lambdas2 <- c((lam=Re(eigen(K_N)$values[1])),(lam=Re(eigen(K_C)$values[1])),(lam=Re(eigen(K_S)$values[1]))) # new population growth rate
		lambdas1 - lambdas2
		setwd(path.obj)
		write.csv(lambdas2, file="CIs_for_lambda.csv") # save for later
		setwd(path.fig)

		
		library(IPMpack)
		# Will convert our objects to work with IPMpack functions 
		# keep all parameters the same as above 
		Pmat_N = new("IPMmatrix", nDiscrete = 0, nEnvClass = 0,
			nBigMatrix = n, nrow = n, ncol = n, meshpoints = y,
			env.index = 0)
		Pmat_C = new("IPMmatrix", nDiscrete = 0, nEnvClass = 0,
			nBigMatrix = n, nrow = n, ncol = n, meshpoints = y,
			env.index = 0)
		Pmat_S = new("IPMmatrix", nDiscrete = 0, nEnvClass = 0,
			nBigMatrix = n, nrow = n, ncol = n, meshpoints = y,
			env.index = 0)		
		
		Pmat[, ] = P_N
		plot(y,meanLifeExpect(Pmat), xlab="Size (t)",ylab="Time", main="N", ylim=c(0,5))
		Pmat[, ] = P_C
		plot(y,meanLifeExpect(Pmat), xlab="Size (t)",ylab="Time", main="C", ylim=c(0,5))
		Pmat[, ] = P_S
		plot(y,meanLifeExpect(Pmat), xlab="Size (t)",ylab="Time", main="S", ylim=c(0,5))
		
	# Plot 
		par(mfrow=c(1,3),mar=c(4,5,2,2))
		image.plot(y,y,t(P_N), xlab="Size (t)",ylab="Size (t+1)",
			col=topo.colors(100), main="IPM matrix N")
			contour(y,y,t(P_N), add = TRUE, drawlabels = TRUE)
			abline(0,1,lwd=3,lty=2)
		image.plot(y,y,t(P_C), xlab="Size (t)",ylab="Size (t+1)",
			col=topo.colors(100), main="IPM matrix C")
			contour(y,y,t(P_C), add = TRUE, drawlabels = TRUE)
			abline(0,1,lwd=3,lty=2)
		image.plot(y,y,t(P_S), xlab="Size (t)",ylab="Size (t+1)",
			col=topo.colors(100), main="IPM matrix S")
			contour(y,y,t(P_S), add = TRUE, drawlabels = TRUE)
			abline(0,1,lwd=3,lty=2)
				
		dfcleanN <- subset(dfclean, Region == "N")
		dfcleanC <- subset(dfclean, Region == "C")
		dfcleanS <- subset(dfclean, Region == "S")

			plot(density(dfcleanN$z1[!is.na(dfcleanN$z1)]),xlab="Size(t+1)",
		main="Observed distribution of sizes N", ylim=c(0,0.3))
				plot(density(dfcleanC$z1[!is.na(dfcleanC$z1)]),xlab="Size(t+1)",
		main="Observed distribution of sizes C", ylim=c(0,0.3))
				plot(density(dfcleanS$z1[!is.na(dfcleanS$z1)]),xlab="Size(t+1)",
		main="Observed distribution of sizes S", ylim=c(0,0.3))
		
		dev.off()
############	