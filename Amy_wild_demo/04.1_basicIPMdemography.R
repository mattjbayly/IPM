############################################################################################
##
##      M A I N     I P M    S C R I P T 
##		
##		Modified for Mimulus cardinalis dataset (April 2015)
## 
###########################################################################################

## INDEX: 
	## Section 1 - Fit statistical models to simulated data and plot these where appropriate
		## 1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)
		## 1.2 - growth (continuous variable 'z1')
		## 1.3A - reproduction (binary indicator 'Repr' = 1 if reproduced)
		## 1.3B - recruitment (binary indicator 'Recr' = 1 if recruitment)
		## 1.4 - recruit size (continuous variable 'Rcsz')



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Ungulate IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(doBy)
require(car)
require(mgcv)
library(coefplot2)

rm(list=ls(all=TRUE))

set.seed(270875)

## working directory must be set here, so the source()'s below run
	setwd("C:/Users/DW/Desktop/transplant_analysis/Planning_Docs/2.IPM_tutorials/Rees_2014_how to IPM/Reese example")

## run the utility functions
	source("./Standard Graphical Pars.R")

## run the ungulate IBM
	source("./Ungulate Demog Funs.R") # but will not use these. 
	

## set the simulation parameters
	init.pop.size <- 500
	n.yrs <- 400
	m.par <- m.par.true # fixed values from "real population"

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## IMPORT CARDINALIS DATA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	### _Set directories for computer_ ###
	path.root="C:/Users/DW/Desktop/transplant_analysis"	
	setwd(path.root)
			path.dat=paste(path.root, "/Data/Amy_wild_demo_data", sep="") # MAIN DATA FILES
			path.code=paste(path.root, "/Rcode/Amy_wild_demo", sep="") # Source Code
			path.obj=paste(path.root, "/Robjects/Amy_wild_demo", sep="") # To store temporary objects
			path.fig=paste(path.root, "/Figures/Amy_wild_demo", sep="") # Polished & exploratory figures
			setwd(path.dat); setwd(path.code); setwd(path.fig); setwd(path.obj)
			path.funct = paste(path.root,"/Rcode/Amy_wild_demo/functions", sep="")
			setwd(path.dat); dir()
	#######################################
	dfclean <- read.csv(file="demo_noMajOutlie.csv")
	# log transformation of basic size data 
	dfclean$z <- log(dfclean$z + 1)
	dfclean$z1 <- log(dfclean$z1 + 1)

	# convert regions down to just three levels N, C & S and make year a facot variable 
	library(plyr)
	dfclean$Region <- revalue(dfclean$Region, c("C1"="C", "C2"="C", "C3"="C", "N1"="N", "N2"="N", "S1"="S", "S2"="S"))
	dfclean$Year <- as.factor(dfclean$Year)
	# also want to round of fruit counts for plants to whole numbers so that they work with poisson distribution 
	dfclean$Fec <- round(dfclean$Fec, digits = 0)
	dfclean$Fec <- as.integer(dfclean$Fec)
	plot.range <- c(0, 9)
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1A - Load vital rate regression models & coefficients from scripts 04.2.1 - 04.2.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## sort by size and print a sample to the screen
	dfclean <- dfclean[order(dfclean$z),]					# sort smallest -> largest (z) 

	m.par <- c(
	surv.int = NA,
	surv.z =  NA,
	grow.int = NA,
	grow.z = NA,
	grow.sd =  NA,
	repr.int = NA, 
	repr.slope = NA,
	repr.slope2 = NA, 
	FecCountComp.int = NA,
	FecCountComp.slope = NA, 
	FecZIComp.int = NA,
	SeedPerFruit = NA,
	RecrCoef.slope = NA,
	mean.lnrom.RecrSize =  NA,
	stdev.lnrom.RecrSize =  NA
	)

# All IPMs to run 
	m.parN <- m.par
	m.parC <- m.par
	m.parS <- m.par

	
#----------------------------------------
## 1A.1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)
	# 04.2_basicIPM_Surv.R
	# source cardinalis survivorship functions
	setwd(path.code)
	#source("04.2.1_Surv.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	surv_model <- get(load("SurvMod6.rda")) # load in top survival mixed model 
	surv_model_boot_coef <- read.csv("SurvModBootCoef.csv") # load in boot strapped coefficients estimates for survivorship model 
	# prepare these parameters to enter the IPM in the 'm.par.est' function  below
	fixef(surv_model); surv_model_boot_coef[,4] # difference between raw & bootstrapped values 
	m.par.est_Surv <- fixef(surv_model)
	m.par.est_Surv[1:6] <- surv_model_boot_coef[1:6,4] # replace values with bootstrapped estimates
	# assign values 
		m.parC["surv.int"] <- m.par.est_Surv["(Intercept)"]; m.parC["surv.z"] <- m.par.est_Surv["z"]
		m.parN["surv.int"] <- m.par.est_Surv["(Intercept)"]+m.par.est_Surv["RegionN"]; m.parN["surv.z"] <- m.par.est_Surv["z"]+m.par.est_Surv["z:RegionN"]
		m.parS["surv.int"] <- m.par.est_Surv["(Intercept)"]+m.par.est_Surv["RegionS"]; m.parS["surv.z"] <- m.par.est_Surv["z"]+m.par.est_Surv["z:RegionS"]
	
	
#----------------------------------------
## 1A.2 - growth (continuous variable 'z1')
	# Growth (given survival), conditional on survival 
	#source(".R") # only run this once on the first time (takes a while)
	setwd(path.code)
	#source("04.2.2_Grow.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	grow_model <- get(load("GrowMod.rda")) # load in top growth mixed model 
	grow_model_boot_coef <- read.csv("GrowModBootCoef.csv") # load in boot strapped coefficients estimates for growth model 
	# prepare these parameters to enter the IPM in the 'm.par.est' function below 
	fixef(grow_model); grow_model_boot_coef[,4] # difference between raw & bootstrapped values 
	m.par.est_Grow <- fixef(grow_model); m.par.est_Grow[1:6] <- grow_model_boot_coef[1:6,4] # replace values with bootstrapped estimates
	# assign values 
		m.parC["grow.sd"]<- sigma(grow_model); m.parN["grow.sd"]<- sigma(grow_model); m.parS["grow.sd"]<- sigma(grow_model)
		m.parC["grow.int"] <- m.par.est_Grow["(Intercept)"]; m.parC["grow.z"] <- m.par.est_Grow["z"]
		m.parN["grow.int"] <- m.par.est_Grow["(Intercept)"]+m.par.est_Grow["RegionN"]; m.parN["grow.z"] <- m.par.est_Grow["z"]+m.par.est_Grow["z:RegionN"]
		m.parS["grow.int"] <- m.par.est_Grow["(Intercept)"]+m.par.est_Grow["RegionS"]; m.parS["grow.z"] <- m.par.est_Grow["z"]+m.par.est_Grow["z:RegionS"]


#----------------------------------------	
## 1A.5 - recruitment model (probability that a fruit or seed dropped from a plant 
	# becomes a new recruit the next year 
	setwd(path.code)
	#source("03_RecruitmentModelComp.R") # only run this once on the first time 
	setwd(path.obj)
	Rec_model <- read.csv(file="RecruitmentModCoefs.csv")
	# assign values 
		m.parC["SeedPerFruit"]<- Rec_model$meanSeedcount[1]; m.parN["SeedPerFruit"]<- Rec_model$meanSeedcount[2]; m.parS["SeedPerFruit"]<- Rec_model$meanSeedcount[3]
		m.parC["RecrCoef.slope"] <- Rec_model$est[1]; m.parN["RecrCoef.slope"] <- Rec_model$est[2]; m.parS["RecrCoef.slope"] <- Rec_model$est[3]


#----------------------------------------	
## 1A.6 - Size distribution of new recruits. 
	setwd(path.code)
	#source("02_basicRecruitmentMod.R") # only run this once on the first time 
	setwd(path.obj)
	Rec_dist <- read.csv(file="lnormFecKern.csv")
	require(stringr)
	Rec_dist$Reg <- str_extract(Rec_dist$tempy,"[[:upper:]]")
	# assign values 
	library(plyr)
	newRec_dist <- ddply(Rec_dist, .(Reg), summarize, meanlog=mean(meanlog), sdlog=mean(sdlog))
	# assign values 
		m.parC["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[1]; m.parN["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[2]; m.parS["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[3]
		m.parC["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[1]; m.parN["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[2]; m.parS["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[3]	
	
###################
# will do IPM matrices for now averaging over all yearly coefficient estimates rather than a separate
# model for each year

#----------------------------------------
## 1A.3 - Flowering (Reproductive?) (binary variable 'Repr')
	# Was plant reproductive (0/1) given survival and size? 
	setwd(path.code)
	#source("04.2.3_Flower.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Repr_model <- get(load("Reprmod3.rda")) # load in top reproductive glm mixed model 
	Repr_model_boot_coef <- read.csv("ReprModBootCoef.csv") # load in boot strapped coefficients estimates for reproductive model 
		# prepare these parameters to enter the IPM in the 'm.par.est' function below 
	fixef(Repr_model); Repr_model_boot_coef[,4] # difference between raw & bootstrapped values 
	m.par.est_Repr <- fixef(Repr_model); m.par.est_Repr[1:9] <- Repr_model_boot_coef[1:9,4] # replace values with bootstrapped estimates
	# Average over slope & int 
	AvgInt <- mean(c(m.par.est_Repr["(Intercept)"], 
	(m.par.est_Repr["(Intercept)"]+m.par.est_Repr["Year1112"]),
	(m.par.est_Repr["(Intercept)"]+m.par.est_Repr["Year1213"])))
	AvgSlope <- mean(c(m.par.est_Repr["z"], 
	(m.par.est_Repr["z"]+m.par.est_Repr["z:Year1112"]),
	(m.par.est_Repr["z"]+m.par.est_Repr["z:Year1213"])))
	# assign values 
		m.parC["repr.int"] <- AvgInt; m.parC["repr.slope"] <- AvgSlope; m.parC["repr.slope2"] <- m.par.est_Repr["I(z^2)"]
		m.parN["repr.int"] <- AvgInt + m.par.est_Repr["RegionN"]; m.parN["repr.slope"] <- AvgSlope; m.parN["repr.slope2"] <- m.par.est_Repr["I(z^2)"]
		m.parS["repr.int"] <- AvgInt + m.par.est_Repr["RegionS"]; m.parS["repr.slope"] <- AvgSlope; m.parS["repr.slope2"] <- m.par.est_Repr["I(z^2)"]

	
#----------------------------------------
## 1A.4 - Fecundity (Given reproductive) (count variable, number of fruits 'Fec')
	# How many fruits were produced (given plant was reproductive 0, 1, 2 .... to max fruits). 
	setwd(path.code)
	#source("04.2.4_Fecundity.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Fec_model <- get(load("FecMod.rda")) # load in top reproductive glm mixed model 
	# don't swap estimated values as we did above. 
	library(coefplot2)
	fecTable <- coeftab(Fec_model,ptype="fixef")
	# Average over slope & int for years
	estimates <- unlist(fecTable[1]); names(estimates) <- rownames(fecTable)
	AvgInt <- mean(c(estimates["Sol.traitFec"], 
	(estimates["Sol.traitFec"]+estimates["Sol.at.level(trait, 1):Year1112"]),
	(estimates["Sol.traitFec"]+estimates["Sol.at.level(trait, 1):Year1213"])))
	AvgSlope <- mean(c(estimates["Sol.at.level(trait, 1):z"], 
	(estimates["Sol.at.level(trait, 1):z"]+estimates["Sol.at.level(trait, 1):Year1112:z"]),
	(estimates["Sol.at.level(trait, 1):z"]+estimates["Sol.at.level(trait, 1):Year1213:z"])))
	# assign values 
		m.parS["FecCountComp.int"] <- AvgInt; m.parS["FecCountComp.slope"] <- AvgSlope; m.parS["FecZIComp.int"] <- estimates["Sol.traitzi_Fec"]
		m.parN["FecCountComp.int"] <- AvgInt + estimates["Sol.at.level(trait, 1):RegionN"]; m.parN["FecCountComp.slope"] <- AvgSlope; m.parN["FecZIComp.int"] <- estimates["Sol.traitzi_Fec"]
		# watch out mcmc dosen't have the same dummy variables as above 
		m.parC["FecCountComp.int"] <- AvgInt + estimates["Sol.at.level(trait, 1):RegionC"]; m.parC["FecCountComp.slope"] <- AvgSlope; m.parC["FecZIComp.int"] <- estimates["Sol.traitzi_Fec"]



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1B - store the estimate parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#DELETE THIS ~ QUICK TEST WITH CENTRAL REGION 
	m.par <- c(
	surv.int = -2.1485265,
	surv.z = 0.6766322,
	grow.int = 1.46555964,
	grow.z = 0.69909474,
	grow.sd =  1.02206,
	repr.int = -5.2466345, 
	repr.slope = 1.0494966,  # 1.3494966
	FecCountComp.int = -1.198446,
	FecCountComp.slope = 0.30549, # 0.60549
	FecZIComp.int = -3.66113,
	SeedPerFruit = 1600,
	RecrCoef.slope = 0.07,
	mean.lnrom.RecrSize =  -0.009299, # -0.016299
	stdev.lnrom.RecrSize =   0.9747318 # 0.7747318
	)


Regions <- c("Central", "North", "South")
pars <-c("m.parC", "m.parN", "m.parS")
	
	# save results in a figure
	setwd(path.fig)
	pdf(file="04.1_IPMbasicALL.pdf", width=11, height=8.5)


for(i in 1:length(pars)){
	m.par <- get(pars[i])
	

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1C - Define basic life history functions for a manual IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# code modified from:
		# Reese et al 2014 - cited above. Some steps also taken from Reese et al 2006
		# Also some steps followed from Merow IPM pack appendices S1 section 6.4.2 "Define functions to describe life history" (pg.153).
		# Separate kernel run for each unique environmental covariate combination.  

		#----------------------------------------
		## 1B.1 - SURVIVAL function:  
			SurvFrom_z <- function(z, m.par)
		{
			linear.p <- m.par["surv.int"] + m.par["surv.z"] * z       # linear predictor (linear.p = intercept + B(initial size))
			p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale) 
			return(p)												  # returns p on logit scale 
		}

			
		#----------------------------------------
		## 1B.2 - GROWTH function:
			Grow_z1z <- function(z1, z, m.par)
		{
			mu <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year (mu = intercept + B(initial size))
			sig <- m.par["grow.sd"]                                 # sd about mean (sig = stdev of growth)
			p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # probability density fxn (pdf) that you are size z1 given you were size z
			return(p.den.grow)										# returns probability density function of growth 
		}
		# end of P-kernel functions, start the F-kernel functions 
			
			
		#----------------------------------------
		## 1B.3 - FLOWERING function (becomes reproductive?):	

		Prob_repr <-function(z,m.par) {
			# Probability of flowering based on size (0 to 1)
				z2 <- z*z # for quadratic term 
				linear.p <- m.par["repr.int"] + (m.par["repr.slope"] * z) + (m.par["repr.slope2"]*z2)  # linear predictor (linear.p = intercept + B(initial size))
				p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale)	
				# p = response (the actual probability of flowering)
				return(p)
		}	
		
	###################################################	
	#  RECRUITMENT SIZE DISTRIBUTION FUNCTION: 
			# RecrDist() is a large function with three components: 
		RecrDist <- function(z1,z,m.par) {	
		#----------------------------------------
		## 1B.4 - FECUNDITY function (how many fruits produced?):			
				# Number of fruits given individual flowered: 
				# Recall in the original model we fit a ZIPGLMM (zero inflated have a binary & count component)
				# Making direct count predictions from estimates requires multiplying the probability that 
				# the value wasn't zero (1-pi) and then multiplying that by the count part (1-pi) * mu
				# http://www.ats.ucla.edu/stat/stata/faq/predict_zip.htm
				count_linear = m.par["FecCountComp.int"] + m.par["FecCountComp.slope"]*z # linear prediction for the count part of the model 
				zi_linear = m.par["FecZIComp.int"]	# linear prediction for the zero inflated part of the model 
				zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
				count_number = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 
				# count of fruits for plant of size z 
				#count_number <- count_number[count_number<0]<-0 # manually set fruit count to zero if predicted values are under 0 otherwise individuals will 'absorb' fruits 
		#----------------------------------------
		## 1B.5 - RECRUITMENT function (how many recruits are produced from these fruits):
			# this function is nested with the the fecundity function 
			# recall we had to log transform both the fruit count & the count of recruits in 
			# the previous script '03_RecruitModelComp.r' so the formula looked like this:
			# log(number of recruits + 1) = (intercept set to zero) + [log(seed count + 1)*log(fruit count + 1)]*coefficient
			# so to get the recruit number we have to exponentialte and take log of other values
			# Number of recruits = exp((log(seeds per fruit +1)*log(numb of fruits +1)*(coef)))	
			recr_count = exp((log(m.par["SeedPerFruit"] + 1)*log(count_number + 1)) * m.par["RecrCoef.slope"])
	#----------------------------------------
		## 1B.6 - SIZE DISTRIBUTION OF NEW RECRUITS	
			# density from the log-normal distribution with
			# z1=vector, mean and standard deviation of the distribution on the log scale)
			NewRecruits =dlnorm(z1,m.par["mean.lnrom.RecrSize"],m.par["stdev.lnrom.RecrSize"])/ # 
						(1-plnorm(minsize,m.par["mean.lnrom.RecrSize"],m.par["stdev.lnrom.RecrSize"]))#
			NewRecruitDist = NewRecruits*recr_count
			return(NewRecruitDist) # size distribution of the new recruits 
		}	

		
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1D - Construct kernels and define P, F & K pdfs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(fields)
	
	# KERNEL PROPERTIES  
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
		par(mfrow=c(2,2))
		G = array(0,dim=c(n.size,n.size)); dim(G); class(G); #image.plot(G)
		## runs Grow_z1z for each combination of possible sizes
		G = h*outer(y,y,Grow_z1z,m.par=m.par); dim(G); str(G); class(G); image.plot(G)
		# the function 'outer()' evaluates the matrix at all pairwise combinations of the two vectors y & y and returns matrices representing the kernel components of the growth 
		# make Vectors from starting size: survive, prob repro, number rec
		S = SurvFrom_z(y,m.par=m.par) # survival vector 
		ProbRepVec= Prob_repr(y,m.par=m.par) # vector of probability of flowering 
		###
		RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
			dim(RecMatrix); class(RecMatrix); image.plot(RecMatrix)
	
	# Growth matrix 
		P=G 						     # placeholder; redefine P on the next line
		for(i in 1:n) P[,i]=G[,i]*S[i]	 # growth/survival matrix
	# Fecundity matrix 
		F=RecMatrix
		for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix
			# this is equivalent to adding "*Prob_repr(z,m.par)" to the last line of the RecrDist function 
			# recall in the fecundity model we only included observations from individuals that were already flowering 
		K=P+F # full matrix, add matrices together 
		image.plot(K)
	# basic lambda estimate	
		(lam <- Re(eigen(K)$values[1])) #  0.939


# next steps are basic components with script from 
# IPMpack S1: Appendix 1
##################################################
# BASIC ANALYSIS 
	# Asymmetric population growth rate
		(lam1 <- Re(eigen(K)$values[1])) # 1.04

	# Make elasticity and sensitivty matrix 
		w.eigen <- Re(eigen(K)$vectors[,1])
		stable.dist <- w.eigen/sum(w.eigen)
		v.eigen <- Re(eigen(t(K))$vectors[,1])
		repro.val <- v.eigen/v.eigen[1]
		v.dot.w=sum(stable.dist*repro.val)*h
		sens=outer(repro.val,stable.dist)/v.dot.w
		elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)
	
	# plot out the IPM components 
		par(mfrow=c(2,3),mar=c(4,5,2,2))
		image.plot(y,y,t(K), xlab="Size (t)",ylab="Size (t+1)",
		col=topo.colors(100), main=paste("IPM matrix", Regions[i], sep=" "))
		contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
		plot(y,stable.dist,xlab="Size",type="l",main=paste("Stable Size Dist", Regions[i], sep=" "))
		plot(y,repro.val,xlab="Size",type="l",main=paste("RVs ", Regions[i], sep=" "))
		image.plot(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main=paste("Elasticity", Regions[i], sep=" "))
		image.plot(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main=paste("Sensitivity", Regions[i], sep=" "))
	
##################################################
# EVICTION	
	# Check for eviction 	
		par(mfrow=c(1,1))
		plot(y,SurvFrom_z(y,m.par),xlab="Size",type="l",
		ylab="Survival Probability", lwd=8, main="Quick Check for eviction")
		points(y,apply(P,2,sum),col="red",lwd=3,cex=.1,pch=19) 
		text(6, 0.85, "preds from survival model", col="black")
		text(6, 0.55, "col sums of P matrix", col="red")
		# for smaller individuals 
		plot(y,SurvFrom_z(y,m.par),xlab="Size",type="l", xlim=c(0,1.5), ylim=c(0,0.2), ylab="Survival Probability", lwd=8, main="Quick Check for eviction")
		points(y,apply(P,2,sum),col="red",lwd=3,cex=.1,pch=19) 
		# One way to correct for eviction is to return the evicted individuals to the cells at 
		# the boundaries where they were evicted. All individuals smaller than the lower integration
		# limit are assigned to the smallest size class. This mainly affects the offspring size distribution; 
		# it sets a lower limit on offspring size.
		# To correct for eviction we need to adjust the matrices before combining them.  
			# fix eviction of offspring
			
		G<-h*outer(y,y,Grow_z1z,m.par=m.par)		# growth matrix 
		S=SurvFrom_z(y,m.par=m.par) 				# survival vector 
		ProbRepVec= Prob_repr(y,m.par=m.par) 		# flowering vector 
		RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
		P=G; F=RecMatrix 							# place holders 
				
		# fix eviction for offspring
			for(i in 1:(n/2)) {
			G[1,i]<-G[1,i]+1-sum(G[,i])
			P[,i]<-G[,i]*S[i] # growth/survival matrix
			}
		# fix eviction of large adults
			for(i in (n/2+1):n) {
			G[n,i]<-G[n,i]+1-sum(G[,i])
			P[,i]<-G[,i]*S[i]
			}				
		F=RecMatrix
		for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix	
		K=P+F # full matrix
		(lam=Re(eigen(K)$values[1])) # new population growth rate
		lam1 - lam # difference in lambda after correction 
		# Re-plot eviction diagnostic plot to see if changes were effective 
		par(mfrow=c(1,2)); plot(y,SurvFrom_z(y,m.par),xlab="Size",type="l", ylab="Survival Probability", lwd=8, main="Quick Check for eviction")
		points(y,apply(P,2,sum),col="red",lwd=3,cex=.1,pch=19) 
		plot(y,SurvFrom_z(y,m.par),xlab="Size",type="l", xlim=c(0,1.5), ylim=c(0,0.2), ylab="Survival Probability", lwd=8, main="Quick Check for eviction")
		points(y,apply(P,2,sum),col="red",lwd=3,cex=.1,pch=19) 
		# looks good... 
		
#####################################
# PASSAGE TIME 
		library(IPMpack)
		# Will convert our objects to work with IPMpack functions 
		# keep all parameters the same as above 
		Pmat = new("IPMmatrix", nDiscrete = 0, nEnvClass = 0,
		nBigMatrix = n, nrow = n, ncol = n, meshpoints = y,
		env.index = 0)
		Pmat[, ] = P
		str(Pmat)		
		(mle=meanLifeExpect(Pmat))
		plot(y,meanLifeExpect(Pmat), xlab="Size (t)",ylab="Time")
		
		# Plot 
		par(mfrow=c(2,1),mar=c(4,5,2,2))
		image.plot(y,y,t(P), xlab="Size (t)",ylab="Size (t+1)",
		col=topo.colors(100), main=paste("IPM matrix", Regions[i], sep=" "))
		contour(y,y,t(P), add = TRUE, drawlabels = TRUE)
		abline(0,1,lwd=3,lty=2)
		plot(density(dfclean$z1[!is.na(dfclean$z1)]),xlab="Size(t+1)",
		main="Observed distribution of sizes")
		
#####################################
# TRANSIENT DYNAMICS
		(lam=Re(eigen(K)$values[1])) # asymptotic growth rate
		
		# Given the current size distribution, how long would it take to return to steady state under
		# The time scale for reaching the steady state can be estimated from the damping ratio (Caswell 2001, p. 95)
		#  - The ratio of the dominant eigenvalue to the second largest one.
		# Small values (near 1) mean slow transients; higher values
		# means the dominant eigenvalue dictates the dynamics.
		(damp=Re(eigen(K)$values[1])/Re(eigen(K)$values[2])) # damping ratio
		
		initial.pop=runif(100) #random starting population structure
		initial.pop=initial.pop/sum(initial.pop)
		nyears=20
		size.dist=matrix(NA,n,nyears)
		lambda=rep(NA,nyears)
		xold=initial.pop
		for(i in 1:nyears){
			xnew=K%*%xold
			lambda[i]=sum(xnew)/sum(xold)
			size.dist[,i]=xnew
			xold=xnew
		}
		lambda
	
	}	
	# close plotting device
	dev.off()
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	