#05.1_IPM_MicrositeMoisture.R


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use code originally developed for the the Ungulate IBM to illustrate the construction of an IPM
## this code was from Rees et al 2014: Building integral projection models: a user’s guide
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(doBy)
require(car)
require(mgcv)
library(coefplot2)
library(plyr)
rm(list=ls(all=TRUE))
set.seed(270875)

## working directory must be set here, so the source()'s below run
	setwd("C:/Users/DW/Desktop/transplant_analysis/Planning_Docs/2.IPM_tutorials/Rees_2014_how to IPM/Reese example")
## run the utility functions
	source("./Standard Graphical Pars.R")
## run the ungulate IBM
	source("./Ungulate Demog Funs.R") # but will not use these. 
		
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## IMPORT CARDINALIS DATA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SET DIRECTORIES & LOAD FILES ###
	path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
	setwd(path.set)
	source("00_SetDirectories.R") # directory script (edit for your own computer). 
	setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)
	# Open 2015/2014 plant datafiles 
	setwd(path.dat); dir()
	dfclean <- read.csv(file="IPMData_transplant.csv")
	dim(dfclean); colnames(dfclean)
	head(dfclean) # looks good ;)
	
	# log transformation of basic size data 
	dfclean$z <- log(dfclean$z ); summary(dfclean$z)
	dfclean$z1 <- log(dfclean$z1); summary(dfclean$z1)
	# round down fecundity data to actual fruit counts:
	dfclean$Fec <- round(dfclean$Fec, digits = 0) # round 
	dfclean$Fec <- as.integer(dfclean$Fec) # make integer for Poisson regression later.
	plot.range <- c(0, 9) # range of size data to be plotted (~ 0 - 9 log(total stem length))
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ADD ON ENVIRO DATA 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###
setwd(path.dat.raw)
enviro <- read.csv(file="environmental_variables.csv")
dim(enviro);
colnames(enviro)
	enviro <- enviro[ , c("SITE","PLOT","SUBSTRATE","Asp_Expo","Densi_cover","Particle_mean","Particle_even","moist_score","Vert._dist_wat","Horiz_dist_to_wat","tot_area_Seasonal","plot_slope","WatAREA__to_bankAREA")] 
	summary(enviro) # do all variable values seem ok?
	# TRANSFORMATIONS NEEDED FOR NORMALITY
	enviro$Particle_mean <- log(enviro$Particle_mean + 1)
	enviro$Particle_even <- sqrt(enviro$Particle_even + 1)
	enviro$Vert._dist_wat <- log(enviro$Vert._dist_wat + 1)
	enviro$tot_area_Seasonal <- log(enviro$tot_area_Seasonal + 1)
	enviro$WatAREA__to_bankAREA <- log(enviro$WatAREA__to_bankAREA + 1)
	enviro <- subset(enviro, SITE != "WILD")
	levels(enviro$SITE)[levels(enviro$SITE)=="HUNT"] <- "HUNTER"
	enviro$Join <- paste0(enviro$SITE, enviro$PLOT)

	dfcleanCHECK <- merge(dfclean, enviro, by.x='PlotID', by.y='Join', all.x=TRUE, all.y=TRUE)
	dfclean <- dfcleanCHECK 
		
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1A - Load vital rate regression models & coefficients from scripts 04.2.1 - 04.2.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## sort by size and print a sample to the screen
	dfclean <- dfclean[order(dfclean$z),]					# sort smallest -> largest (z) 

	m.par <- c( # m.par <- will become the list of all parameter values from the final models to be used in the IPMs 
	surv.int = NA,
	surv.z =  NA,
	surv.int.Site = NA,
	grow.int = NA,
	grow.z = NA,
	grow.sd =  NA,
	grow.int.Site = NA, 
	grow.z.Site = NA, 
	repr.int = NA, 
	repr.slope = NA,
	repr.int.Site = NA, 
	FecCountComp.int = NA,
	FecCountComp.slope = NA, 
	FecCountComp.int.site = NA, 
	FecCountComp.slope.site = NA, 
	SeedPerFruit = NA,
	RecrCoef.slope = NA,
	mean.lnrom.RecrSize =  NA,
	stdev.lnrom.RecrSize =  NA,
	Surv_moist_score = NA,
	Grow_moist_score = NA,
	Repr_moist_score = NA,
	Fec_moist_score = NA
	)	

	m.parN <- m.par # will duplicate this list & as with the demography dataset with use recruitment estimates from 
	# demography sites 
		m.parC <- m.par; m.parS <- m.par # won't end up using m.parS (South) or m.parC (Center) except for exploratory purposes. 

		

		
#----------------------------------------
## 1A.1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)
	# 04.2_basicIPM_Surv.R
	# source cardinalis survivorship functions
	setwd(path.code)
	#source("04.2.1_Surv.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	surv_model <- get(load("SurvMod6.rda")) # load in top survival mixed model 
	surv_model_boot_coef <- read.csv("SurvModBootCoef.csv") # load in boot strapped coefficients estimates for survivorship model 
	surv_model_boot_coef
	
	# REBUILD MODEL SWAP PLOT INT WITH MICROSITE VAR
	dfclean2 <- dfclean[complete.cases(dfclean$z),]
	surv2 <- glm(Surv ~ z + SiteID + moist_score, family = binomial(link="logit"), data=dfclean2, na.action=na.omit)
	# mean moist_score = 3.5
	
	# prepare these parameters to enter the IPM in the 'm.par.est' function  below
	m.par.est_Surv <- coef(surv2)
	names(m.par.est_Surv) <- paste0("Surv_", names(m.par.est_Surv))
	
	
#----------------------------------------
## 1A.2 - growth (continuous variable 'z1')
	# Growth (given survival), conditional on survival 
	#source(".R") # only run this once on the first time (takes a while)
	setwd(path.code)
	#source("04.2.2_Grow.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	grow_model <- get(load("GrowMod.rda")) # load in top growth mixed model 
	grow_model_boot_coef <- read.csv("GrowModBootCoef.csv") # load in boot strapped coefficients estimates for growth model 
	
	# REBUILD MODEL SWAP PLOT INT WITH MICROSITE VAR
	dfclean2 <- dfclean[complete.cases(dfclean$z),]
	dfclean2 <- dfclean2[complete.cases(dfclean2$z1),]

>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>

	grow2 <- lm(z1 ~ z + SiteID + z*SiteID +  moist_score, data=dfclean2, na.action=na.omit)
	mySd <- summary(grow2)$sigma
	# mean moist_score = 3.5
	
	m.par.est_Grow <- coef(grow2)
	# m.par.est_Grow[1:6] <- grow_model_boot_coef[1:6,4] # replace values with bootstrapped estimates
	# assign values 
		names(m.par.est_Grow) <- paste0("Grow_", names(m.par.est_Grow))

		
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
	# also have to use the size distribution of new recruits 
	# form the demography dataset 
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
	
	
#----------------------------------------
## 1A.3 - Flowering (Reproductive?) (binary variable 'Repr')
	# Was plant reproductive (0/1) given survival and size? 
	setwd(path.code)
	#source("04.2.3_Flower.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Repr_model <- get(load("Reprmod9L.rda")) # load in top reproductive glm mixed model 
	Repr_model_boot_coef <- read.csv("ReprModBootCoef.csv") # load in boot strapped coefficients estimates for reproductive model 
	
	# REBUILD MODEL SWAP PLOT WITH MICROSITE VAR
	dfclean2 <- dfclean[which(dfclean$Year=="2014to2015"),]
	dfclean2 <- dfclean2[complete.cases(dfclean2$Repr),]
	dfclean3 <- dfclean[which(dfclean$Year=="2015to2016"),]
	dfclean3 <- dfclean3[complete.cases(dfclean3$Repr),]
	dfclean3$z <- dfclean3$z1 # neeed to swap years
	dfclean2 <- rbind(dfclean2, dfclean3) # bind years
	repr2 <- glm(Repr ~ z + SiteID + moist_score, family = binomial(link="logit"), data=dfclean2, na.action=na.omit)
	# mean moist_score = 3.5
	# Average over slope & int 
		m.par.est_Repr <- coef(repr2)
		names(m.par.est_Repr) <- paste0("Repr_", names(m.par.est_Repr))


#----------------------------------------
## 1A.4 - Fecundity (Given reproductive) (count variable, number of fruits 'Fec')
	# How many fruits were produced (given plant was reproductive 0, 1, 2 .... to max fruits). 
	setwd(path.code)
	#source("04.2.4_Fecundity.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Fec_model <- get(load("FecMod_NB.rda")) # load in top reproductive glm mixed model 
	
	# REBUILD MODEL SWAP PLOT WITH MICROSITE VAR
	dfclean2 <- dfclean[which(dfclean$Year=="2014to2015"),]
	dfclean2 <- dfclean2[complete.cases(dfclean2$Fec),]
	dfclean3 <- dfclean[which(dfclean$Year=="2015to2016"),]
	dfclean3 <- dfclean3[complete.cases(dfclean3$Fec),]
	dfclean3$z <- dfclean3$z1 # neeed to swap years
	dfclean2 <- rbind(dfclean2, dfclean3) # bind years


>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>

	Fec2 <- glm(Fec ~ z + SiteID + z*SiteID+ moist_score, family = poisson(link="log"), data=dfclean2, na.action=na.omit)
	
	# don't swap estimated values as we did above. 
	m.par.est_Fec <- coef(Fec2)
	names(m.par.est_Fec) <- paste0("Fec_", names(m.par.est_Fec))

	
########################################################################################	
########################################################################################
# Fill m.par lists with site values.
 
	# will have to add values for "CALAPOOIA" first
	m.par.CALAPOOIA <- m.parN # load missing values from northern region 
	m.par.CALAPOOIA["surv.int"] =  m.par.est_Surv["Surv_(Intercept)"]
	m.par.CALAPOOIA["surv.z"] = m.par.est_Surv["Surv_z"]
	m.par.CALAPOOIA["surv.int.Site"] = 0
	m.par.CALAPOOIA["grow.int"] = m.par.est_Grow["Grow_(Intercept)"]
	m.par.CALAPOOIA["grow.z"] = m.par.est_Grow["Grow_z"]
	m.par.CALAPOOIA["grow.int.Site"] = 0
	m.par.CALAPOOIA["grow.z.Site"] = 0
	m.par.CALAPOOIA["grow.sd"] = summary(grow2)$sigma
	m.par.CALAPOOIA["repr.int"] = m.par.est_Repr["Repr_(Intercept)"]
	m.par.CALAPOOIA["repr.slope"] = m.par.est_Repr["Repr_z"]
	m.par.CALAPOOIA["repr.int.Site"] = 0
	m.par.CALAPOOIA["FecCountComp.slope.site"] = 0
	m.par.CALAPOOIA["FecCountComp.int.site"] = 0
	m.par.CALAPOOIA["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
	m.par.CALAPOOIA["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
	m.par.CALAPOOIA["Surv_moist_score"] = m.par.est_Surv["Surv_moist_score"]
	m.par.CALAPOOIA["Grow_moist_score"] = m.par.est_Grow["Grow_(Intercept)"]
	m.par.CALAPOOIA["Repr_moist_score"] = m.par.est_Repr["Repr_(Intercept)"]
	m.par.CALAPOOIA["Fec_moist_score"] = m.par.est_Fec["Fec_moist_score"]
		
	to_loop <- levels(dfclean$SiteID)[2:length(levels(dfclean$SiteID))] #all sites but the one used as a dummy variable. 
	to_loop
	# loop to fill m.par values for sites 
	for(i in 1:length(to_loop)){
		temp <- m.parN # temporary obj to fill 
		m.par.temp <- m.parN # load missing values from northern region 
		m.par.temp["surv.int"] =  m.par.est_Surv["Surv_(Intercept)"]
		m.par.temp["surv.z"] = m.par.est_Surv["Surv_z"]
		m.par.temp["surv.int.Site"] = m.par.est_Surv[paste0("Surv_SiteID", to_loop[i])]
		m.par.temp["grow.int"] = m.par.est_Grow["Grow_(Intercept)"]
		m.par.temp["grow.z"] = m.par.est_Grow["Grow_z"]
		m.par.temp["grow.int.Site"] = m.par.est_Grow[paste0("Grow_SiteID", to_loop[i])]
		m.par.temp["grow.z.Site"] = m.par.est_Grow[paste0("Grow_z:SiteID", to_loop[i])]
		m.par.temp["grow.sd"] = summary(grow2)$sigma
		m.par.temp["repr.int"] = m.par.est_Repr["Repr_(Intercept)"]
		m.par.temp["repr.slope"] = m.par.est_Repr["Repr_z"]
		m.par.temp["repr.int.Site"] = m.par.est_Repr[paste0("Repr_SiteID", to_loop[i])]
		m.par.temp["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
		m.par.temp["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
		m.par.temp["Surv_moist_score"] = m.par.est_Surv["Surv_moist_score"]
		m.par.temp["Grow_moist_score"] = m.par.est_Grow["Grow_(Intercept)"]
		m.par.temp["Repr_moist_score"] = m.par.est_Repr["Repr_(Intercept)"]
		m.par.temp["Fec_moist_score"] = m.par.est_Fec["Fec_moist_score"]
	# m.par.temp
		assign(paste0("m.par.", to_loop[i]), m.par.temp)
	}
	#**** suspect values: FecCountComp.slope.site values for (MOSBY & LOOK) 
	#**** suspect values: FecCountComp.int.site values for (MOSBY & LOOK) 



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
			linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  +  m.par["surv.int.Site"] + (m.par["Surv_moist_score"]*0)   # linear predictor (linear.p = intercept + B(initial size))
			p <- 1/(1+exp(-linear.p))                                 						 # logistic transformation to probability (p - logit scale) 
			return(p)																	     # returns p on logit scale 
		}

			
		#----------------------------------------
		## 1B.2 - GROWTH function:
			Grow_z1z <- function(z1, z, m.par)
		{
			mu <- m.par["grow.int"] + m.par["grow.z"] * z  +  m.par["grow.int.Site"] +  m.par["grow.z.Site"]*z*1 + (m.par["Grow_moist_score"]*0)       # mean size next year (mu = intercept + B(initial size))
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
				linear.p <- m.par["repr.int"] + (m.par["repr.slope"] * z) + m.par["repr.int.Site"] + (m.par["Repr_moist_score"]*0)  # linear predictor (linear.p = intercept + B(initial size) + plus site level intercept)
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
				count_number= exp((m.par["surv.int"]) + (m.par["FecCountComp.int.site"]) +
									(m.par["surv.z"]*z) + (m.par["FecCountComp.slope.site"]*z*1) + (m.par["Fec_moist_score"]*0))					

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
	setwd(path.fig)
	pdf(file="05.1_IPMbasicALL.pdf", width=11, height=8.5)

	# Run a loop through all sites & separate IPM model for each site 
	all_sites <-  c("CALAPOOIA", to_loop) # add on calapooia as a site 
	siteLambdas <- matrix(NA, nrow = length(all_sites), ncol = 4)

		for(p in 1:length(all_sites)){
			m.par <- get(paste0("m.par.", all_sites[p])) # set parameter list according to site

	# KERNEL PROPERTIES  
		# number of classes, or points for midpoint rule approximation
		n.size = 100
		# tolerance for iterations
		tol = 1.e-8
		# minimum and maximum sizes (0.9*min & 1.1*max size from data)
		minsize = 0.9*(min(rbind(dfclean$z1, dfclean$z), na.rm=TRUE)) # 0.9 to go slightly lower than smallest size
		minsize <- ifelse(minsize==0, 0.1, minsize) #is minsize zero, cant be  
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
		G = array(0,dim=c(n.size,n.size)); dim(G); class(G); image.plot(G, main=paste0(all_sites[p], " growth"))
		## runs Grow_z1z for each combination of possible sizes
		G = h*outer(y,y,Grow_z1z,m.par=m.par); dim(G); str(G); class(G); image.plot(G, main=paste0(all_sites[p], " growth"))
		# the function 'outer()' evaluates the matrix at all pairwise combinations of the two vectors y & y and returns matrices representing the kernel components of the growth 
		# make Vectors from starting size: survive, prob repro, number rec
		S = SurvFrom_z(y,m.par=m.par) # survival vector 
		ProbRepVec= Prob_repr(y,m.par=m.par) # vector of probability of flowering 
		###
		RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
			dim(RecMatrix); class(RecMatrix); image.plot(RecMatrix, main=paste0(all_sites[p], " Rec"))
	
	# Growth matrix 
		P=G 						     # placeholder; redefine P on the next line
		for(i in 1:n) P[,i]=G[,i]*S[i]	 # growth/survival matrix
	# Fecundity matrix 
		F=RecMatrix
		for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix
			# this is equivalent to adding "*Prob_repr(z,m.par)" to the last line of the RecrDist function 
			# recall in the fecundity model we only included observations from individuals that were already flowering 
		K=P+F # full matrix, add matrices together 
		image.plot(K, main=paste0(all_sites[p], " F & P matrix K "))
	# basic lambda estimate	
		(all_sites[p])
		siteLambdas[p,1] <- all_sites[p] 
		(lam <- Re(eigen(K)$values[1])) #  0.939
		siteLambdas[p,2] <- Re(eigen(K)$values[1])
		
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
		col=topo.colors(100), main=paste("IPM matrix", all_sites[p],  sep=" "))
		contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
		plot(y,stable.dist,xlab="Size",type="l",main=paste("Stable Size Dist",  all_sites[p], sep=" "))
		plot(y,repro.val,xlab="Size",type="l",main=paste("RVs ",  all_sites[p], sep=" "))
		image.plot(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main=paste("Elasticity",  all_sites[p], sep=" "))
		image.plot(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main=paste("Sensitivity",  all_sites[p], sep=" "))
		plot.new()
	
		
	
	
##################################################
# EVICTION	
	# Check for eviction 	
		par(mfrow=c(1,2))
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
		siteLambdas[p,3] <- Re(eigen(K)$values[1])

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
		col=topo.colors(100), main=paste("IPM matrix", all_sites[p], sep=" "))
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
		siteLambdas[p,4] <- damp	
		
	}

siteLambdas
dev.off()

setwd(path.obj)
write.csv(siteLambdas, file="siteLambdas_MoistScore0_noInt.csv")
