# 03.1_IPM_bootstrapped_CIs.R
############################################################################################
##
##      I P M   CONFIDENCE INTERVALS FROM BOOTSTRAP REPLICATES 
##		
## 
###########################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use code originally developed for the the Ungulate IBM to illustrate the construction of an IPM
## this code was from Rees et al 2014: Building integral projection models: a user’s guide
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(doBy)
require(car)
require(mgcv)
library(coefplot2)
library(plyr)
library(IPMpack)
library(lme4)
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
	stdev.lnrom.RecrSize =  NA
	)	

	m.parN <- m.par # will duplicate this list & as with the demography dataset with use recruitment estimates from 
	# demography sites 
	#m.parC <- m.par; m.parS <- m.par # won't end up using m.parS (South) or m.parC (Center) except for exploratory purposes. 

########################################################################################
	# LOAD BOOTSRAP REPLICATE COEFFICENTS
	setwd(path.obj)
	grow_boot <- get(load('GrowModBootReplicates500.rda'))
	surv_boot <- get(load('SurvModBootReplicates500.rda'))
	repr_boot <- get(load('ReprModBootReplicates500.rda'))
########################################################################################	
	# MAKE MATRIX FOR LAMBDA VALUES 
	to_loop <- levels(dfclean$SiteID)[2:length(levels(dfclean$SiteID))] #all sites but the one used as a dummy variable. 
	to_loop
	all_sites <-  c("CALAPOOIA", to_loop) # add on calapooia as a site 
	siteLambdas <- matrix(NA, nrow = length(repr_boot), ncol = 9)
	colnames(siteLambdas) <- c("REP", "CALAPOOIA", "COAST", "HUNTER", "LOOK", "MOSBY", "ROCK", "THOMAS", "WILEY")
	siteLambdas[,1] <- c(seq(1:5000))
	head(siteLambdas)
#####################################################################################
	# load in bootstrap replicates 
	setwd(path.obj)
	BootStrapData <- read.csv(file="BootReplicateDataSets.csv")
########################################################################################

for(u in 1:length(surv_boot)){
		thisRep <- BootStrapData[which(BootStrapData$Replicate==u),]
		SitesInRep <- unique(as.factor(thisRep$SiteID))
		SitesInRep <- sort(as.character(SitesInRep)) # make alphabetical to match dummy variable 
		to_loop <- SitesInRep[2:length(SitesInRep)]
		
#----------------------------------------
## 1A.2 - growth (continuous variable 'z1')
	# Growth (given survival), conditional on survival 
	#source(".R") # only run this once on the first time (takes a while)
	setwd(path.code)
	#source("04.2.2_Grow.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	grow_model <- get(load("GrowMod.rda")) # load in top growth mixed model 

#----------------------------------------	
## 1A.5 - recruitment model (probability that a fruit or seed dropped from a plant 
	# becomes a new recruit the next year 
	setwd(path.code)
	#source("03_RecruitmentModelComp.R") # only run this once on the first time 
	setwd(path.obj)
	Rec_model <- read.csv(file="RecruitmentModCoefs.csv")
	# draw from distribution and assign values 
		m.parN["SeedPerFruit"]<- Rec_model$meanSeedcount[2]
		m.parN["SeedPerFruit"]<- rnorm(1, mean = Rec_model$meanSeedcount[2], sd = Rec_model$Seedsd[2])
		# value must be positive
		m.parN["SeedPerFruit"] <- ifelse(m.parN["SeedPerFruit"] <= 0, Rec_model$meanSeedcount[2], m.parN["SeedPerFruit"])
		# slope for recruitment coefficient 
		m.parN["RecrCoef.slope"] <- runif(1, Rec_model$LL[2], Rec_model$UL[2])
		assign(paste0("m.parN_", u), m.parN)	

	
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
		m.parN["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[2]
		m.parN["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[2]
		assign(paste0("m.parN_", u), m.parN)	

#----------------------------------------
## 1A.4 - Fecundity (Given reproductive) (count variable, number of fruits 'Fec')
	# How many fruits were produced (given plant was reproductive 0, 1, 2 .... to max fruits). 
	setwd(path.code)
	#source("04.2.4_Fecundity.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Fec_model <- get(load("FecMod_NB.rda")) # load in top reproductive glm mixed model 
	# don't swap estimated values as we did above. 
	fixef(Fec_model)
	library(coefplot2)
	m.par.est_Fec <- fixef(Fec_model)
	names(m.par.est_Fec) <- paste0("Fec_", names(m.par.est_Fec))	
		
 	SitesInRep
	thisSurv <- surv_boot[[u]]
   	thisGrow <- grow_boot[[u]]
   	thisRepr <- repr_boot[[u]]

	
	# Fill m.par lists with site values.
	# will have to add values for "DUMMY" first ~ DUMMY, BUT CHANGE NAME 
	m.par.DUMMY <- m.parN # load missing values from northern region 
	m.par.DUMMY["surv.int"] =  thisSurv["(Intercept)"]
	m.par.DUMMY["surv.z"] = thisSurv["z"]
	m.par.DUMMY["surv.int.Site"] = 0
	m.par.DUMMY["grow.int"] = thisGrow["(Intercept)"]
	m.par.DUMMY["grow.z"] = thisGrow["z"]
	m.par.DUMMY["grow.int.Site"] = 0
	m.par.DUMMY["grow.z.Site"] = 0
	m.par.DUMMY["grow.sd"] = sigma(grow_model)
	m.par.DUMMY["repr.int"] = thisRepr["(Intercept)"]
	m.par.DUMMY["repr.slope"] = thisRepr["z"]
	m.par.DUMMY["repr.int.Site"] = 0
		if(SitesInRep[1]=="CALAPOOIA") {
			m.par.DUMMY["FecCountComp.slope.site"] = 0
			m.par.DUMMY["FecCountComp.int.site"] = 0
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				} else {
			m.par.DUMMY["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				}
	m.par.DUMMY
	assign(paste0("m.par.", SitesInRep[1], "_", u), m.par.DUMMY)
	
	# loop to fill m.par values for sites 
	for(i in 1:length(to_loop)){
		temp <- m.parN # temporary obj to fill 
		m.par.temp <- m.parN # load missing values from northern region 
		m.par.temp["surv.int"] =  thisSurv["(Intercept)"]
		m.par.temp["surv.z"] = thisSurv["z"]
		m.par.temp["surv.int.Site"] = thisSurv[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.int"] = thisGrow["(Intercept)"]
		m.par.temp["grow.z"] = thisGrow["z"]
		m.par.temp["grow.int.Site"] = thisGrow[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.z.Site"] = thisGrow[paste0("z:SiteID", to_loop[i])]
		m.par.temp["grow.sd"] = sigma(grow_model)
		m.par.temp["repr.int"] =thisRepr["(Intercept)"]
		m.par.temp["repr.slope"] = thisRepr["z"]
		m.par.temp["repr.int.Site"] = thisRepr[paste0("SiteID", to_loop[i])]
		m.par.temp["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
		m.par.temp["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
		# m.par.temp
		assign(paste0("m.par.", to_loop[i], "_", u), m.par.temp)
	}	
	(print(u))
	}

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1C - Define basic life history functions for a manual IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------
	## 1B.1 - SURVIVAL function:  
			SurvFrom_z <- function(z, m.par)
		{
			linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  +  m.par["surv.int.Site"]   # linear predictor (linear.p = intercept + B(initial size))
			p <- 1/(1+exp(-linear.p))                                 						 # logistic transformation to probability (p - logit scale) 
			return(p)																	     # returns p on logit scale 
		}		
#----------------------------------------
	## 1B.2 - GROWTH function:
			Grow_z1z <- function(z1, z, m.par)
		{
			mu <- m.par["grow.int"] + m.par["grow.z"] * z  +  m.par["grow.int.Site"] +  m.par["grow.z.Site"]*z*1       # mean size next year (mu = intercept + B(initial size))
			sig <- m.par["grow.sd"]                                 # sd about mean (sig = stdev of growth)
			p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # probability density fxn (pdf) that you are size z1 given you were size z
			return(p.den.grow)										# returns probability density function of growth 
		}	
#----------------------------------------
	## 1B.3 - FLOWERING function (becomes reproductive?):	
		Prob_repr <-function(z,m.par) {
			# Probability of flowering based on size (0 to 1)
				z2 <- z*z # for quadratic term 
				linear.p <- m.par["repr.int"] + (m.par["repr.slope"] * z) + m.par["repr.int.Site"]  # linear predictor (linear.p = intercept + B(initial size) + plus site level intercept)
				p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale)	
				# p = response (the actual probability of flowering)
				return(p)
		}	
#----------------------------------------
	#  RECRUITMENT SIZE DISTRIBUTION FUNCTION: 
		RecrDist <- function(z1,z,m.par) {	
		#----------------------------------------
		## 1B.4 - FECUNDITY function (how many fruits produced?):			
				# Number of fruits given individual flowered: 		
				count_number= exp((m.par["surv.int"]) + (m.par["FecCountComp.int.site"]) +
									(m.par["surv.z"]*z) + (m.par["FecCountComp.slope.site"]*z*1))					

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
	# Run a loop through all sites & separate IPM model for each site 
	siteLambdas <- data.frame(siteLambdas)
	head(siteLambdas)
	siteLambdas <- as.matrix(siteLambdas)
	testObject <- function(object)
			{exists(as.character(substitute(object)))}

	#
	##
	###
	####
	#####
	######
	#######
for(j in 1:5000){	
	for(p in 1:length(all_sites)){
		# test if site replicate was even successfull?
		# also are there any missing values in m.par?
		moot = exists(paste0("m.par.", all_sites[p], "_", j))
	
	if(exists(paste0("m.par.", all_sites[p], "_", j))) {
					# if everythings ok then, 
					m.par <- get(paste0("m.par.", all_sites[p], "_", j)) # set parameter list according to site
					if('TRUE' %in%  is.na(get(paste0("m.par.", all_sites[p], "_", j)))){
								siteLambdas[j,all_sites[p]] <- NA 
						} else {
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
					G = array(0,dim=c(n.size,n.size))
					## runs Grow_z1z for each combination of possible sizes
					G = h*outer(y,y,Grow_z1z,m.par=m.par)
					# the function 'outer()' evaluates the matrix at all pairwise combinations of the two vectors y & y and returns matrices representing the kernel components of the growth 
					# make Vectors from starting size: survive, prob repro, number rec
					S = SurvFrom_z(y,m.par=m.par) # survival vector 
					ProbRepVec= Prob_repr(y,m.par=m.par) # vector of probability of flowering 
					###
					RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
						#dim(RecMatrix); class(RecMatrix)
				
				# Growth matrix 
					P=G 						     # placeholder; redefine P on the next line
					for(i in 1:n) P[,i]=G[,i]*S[i]	 # growth/survival matrix
				# Fecundity matrix 
					F=RecMatrix
					for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix
						# this is equivalent to adding "*Prob_repr(z,m.par)" to the last line of the RecrDist function 
						# recall in the fecundity model we only included observations from individuals that were already flowering 
					K=P+F # full matrix, add matrices together 
				# basic lambda estimate	
					lam <- Re(eigen(K)$values[1]) #  0.939
					siteLambdas[j,all_sites[p]] <- lam 
				}	
				} else {
				siteLambdas[j,all_sites[p]] <- NA 
	
				}	
	
	}
	(print(j))
	}

	
head(siteLambdas, 15)
tail(siteLambdas, 15)
siteLambdas <- data.frame(siteLambdas)
slamsTrim <- siteLambdas
slamsTrim$CALAPOOIA[slamsTrim$CALAPOOIA > 15] <- NA
slamsTrim$COAST[slamsTrim$COAST > 15] <- NA
slamsTrim$HUNTER[slamsTrim$HUNTER > 15] <- NA
slamsTrim$LOOK[slamsTrim$LOOK > 15] <- NA
slamsTrim$MOSBY[slamsTrim$MOSBY > 15] <- NA
slamsTrim$ROCK[slamsTrim$ROCK > 15] <- NA
slamsTrim$THOMAS[slamsTrim$THOMAS > 15] <- NA
slamsTrim$WILEY[slamsTrim$WILEY > 15] <- NA
summary(slamsTrim)

#siteLambdas
dev.off()

setwd(path.obj)
write.csv(siteLambdas, file="siteLambdasBootCI.csv")


#
##
###
####
######
#######
########
#########
##########
###########
###########
############
############# # RUN BOOTSTRAP REPLICATES FOR SITES
############
###########
###########
##########
#########
########
#######
######
#####
####
###
##
#



	
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
	stdev.lnrom.RecrSize =  NA
	)	

	m.parN <- m.par # will duplicate this list & as with the demography dataset with use recruitment estimates from 
	# demography sites 
	#m.parC <- m.par; m.parS <- m.par # won't end up using m.parS (South) or m.parC (Center) except for exploratory purposes. 

########################################################################################
	# LOAD BOOTSRAP REPLICATE COEFFICENTS
	setwd(path.obj)
	grow_boot <- get(load('GrowModBootReplicates500_sites.rda'))
	surv_boot <- get(load('SurvModBootReplicates500_sites.rda'))
	repr_boot <- get(load('ReprModBootReplicates500_sites.rda'))
########################################################################################	
	# MAKE MATRIX FOR LAMBDA VALUES 
	to_loop <- levels(dfclean$SiteID)[2:length(levels(dfclean$SiteID))] #all sites but the one used as a dummy variable. 
	to_loop
	all_sites <-  c("CALAPOOIA", to_loop) # add on calapooia as a site 
	siteLambdas <- matrix(NA, nrow = length(repr_boot), ncol = 9)
	colnames(siteLambdas) <- c("REP", "CALAPOOIA", "COAST", "HUNTER", "LOOK", "MOSBY", "ROCK", "THOMAS", "WILEY")
	siteLambdas[,1] <- c(seq(1:5000))
	head(siteLambdas)
#####################################################################################
	# load in bootstrap replicates 
	setwd(path.obj)
	BootStrapData <- read.csv(file="BootReplicateDataSets_Sites.csv")
########################################################################################

for(u in 1:length(surv_boot)){
		thisRep <- BootStrapData[which(BootStrapData$Replicate==u),]
		SitesInRep <- unique(as.factor(thisRep$SiteID))
		SitesInRep <- sort(as.character(SitesInRep)) # make alphabetical to match dummy variable 
		to_loop <- SitesInRep[2:length(SitesInRep)]
		
#----------------------------------------
## 1A.2 - growth (continuous variable 'z1')
	# Growth (given survival), conditional on survival 
	#source(".R") # only run this once on the first time (takes a while)
	setwd(path.code)
	#source("04.2.2_Grow.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	grow_model <- get(load("GrowMod.rda")) # load in top growth mixed model 

#----------------------------------------	
## 1A.5 - recruitment model (probability that a fruit or seed dropped from a plant 
	# becomes a new recruit the next year 
	setwd(path.code)
	#source("03_RecruitmentModelComp.R") # only run this once on the first time 
	setwd(path.obj)
	Rec_model <- read.csv(file="RecruitmentModCoefs_sites.csv")
	# draw from distribution and assign values 
		m.parN["SeedPerFruit"]<- Rec_model$meanSeedcount[2]
		m.parN["SeedPerFruit"]<- rnorm(1, mean = Rec_model$meanSeedcount[2], sd = Rec_model$Seedsd[2])
		# value must be positive
		m.parN["SeedPerFruit"] <- ifelse(m.parN["SeedPerFruit"] <= 0, Rec_model$meanSeedcount[2], m.parN["SeedPerFruit"])
		# slope for recruitment coefficient 
		m.parN["RecrCoef.slope"] <- runif(1, Rec_model$LL[2], Rec_model$UL[2])
		assign(paste0("m.parN_", u), m.parN)	

	
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
		m.parN["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[2]
		m.parN["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[2]
		assign(paste0("m.parN_", u), m.parN)	

#----------------------------------------
## 1A.4 - Fecundity (Given reproductive) (count variable, number of fruits 'Fec')
	# How many fruits were produced (given plant was reproductive 0, 1, 2 .... to max fruits). 
	setwd(path.code)
	#source("04.2.4_Fecundity.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Fec_model <- get(load("FecMod_NB.rda")) # load in top reproductive glm mixed model 
	# don't swap estimated values as we did above. 
	fixef(Fec_model)
	library(coefplot2)
	m.par.est_Fec <- fixef(Fec_model)
	names(m.par.est_Fec) <- paste0("Fec_", names(m.par.est_Fec))	
		
 	SitesInRep
	thisSurv <- surv_boot[[u]]
   	thisGrow <- grow_boot[[u]]
   	thisRepr <- repr_boot[[u]]

	
	# Fill m.par lists with site values.
	# will have to add values for "DUMMY" first ~ DUMMY, BUT CHANGE NAME 
	m.par.DUMMY <- m.parN # load missing values from northern region 
	m.par.DUMMY["surv.int"] =  thisSurv["(Intercept)"]
	m.par.DUMMY["surv.z"] = thisSurv["z"]
	m.par.DUMMY["surv.int.Site"] = 0
	m.par.DUMMY["grow.int"] = thisGrow["(Intercept)"]
	m.par.DUMMY["grow.z"] = thisGrow["z"]
	m.par.DUMMY["grow.int.Site"] = 0
	m.par.DUMMY["grow.z.Site"] = 0
	m.par.DUMMY["grow.sd"] = sigma(grow_model)
	m.par.DUMMY["repr.int"] = thisRepr["(Intercept)"]
	m.par.DUMMY["repr.slope"] = thisRepr["z"]
	m.par.DUMMY["repr.int.Site"] = 0
		if(SitesInRep[1]=="CALAPOOIA") {
			m.par.DUMMY["FecCountComp.slope.site"] = 0
			m.par.DUMMY["FecCountComp.int.site"] = 0
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				} else {
			m.par.DUMMY["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				}
	m.par.DUMMY
	assign(paste0("m.par.", SitesInRep[1], "_", u), m.par.DUMMY)
	
	# loop to fill m.par values for sites 
	for(i in 1:length(to_loop)){
		temp <- m.parN # temporary obj to fill 
		m.par.temp <- m.parN # load missing values from northern region 
		m.par.temp["surv.int"] =  thisSurv["(Intercept)"]
		m.par.temp["surv.z"] = thisSurv["z"]
		m.par.temp["surv.int.Site"] = thisSurv[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.int"] = thisGrow["(Intercept)"]
		m.par.temp["grow.z"] = thisGrow["z"]
		m.par.temp["grow.int.Site"] = thisGrow[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.z.Site"] = thisGrow[paste0("z:SiteID", to_loop[i])]
		m.par.temp["grow.sd"] = sigma(grow_model)
		m.par.temp["repr.int"] =thisRepr["(Intercept)"]
		m.par.temp["repr.slope"] = thisRepr["z"]
		m.par.temp["repr.int.Site"] = thisRepr[paste0("SiteID", to_loop[i])]
		m.par.temp["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
		m.par.temp["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
		# m.par.temp
		assign(paste0("m.par.", to_loop[i], "_", u), m.par.temp)
	}	
	(print(u))
	}

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1C - Define basic life history functions for a manual IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------
	## 1B.1 - SURVIVAL function:  
			SurvFrom_z <- function(z, m.par)
		{
			linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  +  m.par["surv.int.Site"]   # linear predictor (linear.p = intercept + B(initial size))
			p <- 1/(1+exp(-linear.p))                                 						 # logistic transformation to probability (p - logit scale) 
			return(p)																	     # returns p on logit scale 
		}		
#----------------------------------------
	## 1B.2 - GROWTH function:
			Grow_z1z <- function(z1, z, m.par)
		{
			mu <- m.par["grow.int"] + m.par["grow.z"] * z  +  m.par["grow.int.Site"] +  m.par["grow.z.Site"]*z*1       # mean size next year (mu = intercept + B(initial size))
			sig <- m.par["grow.sd"]                                 # sd about mean (sig = stdev of growth)
			p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # probability density fxn (pdf) that you are size z1 given you were size z
			return(p.den.grow)										# returns probability density function of growth 
		}	
#----------------------------------------
	## 1B.3 - FLOWERING function (becomes reproductive?):	
		Prob_repr <-function(z,m.par) {
			# Probability of flowering based on size (0 to 1)
				z2 <- z*z # for quadratic term 
				linear.p <- m.par["repr.int"] + (m.par["repr.slope"] * z) + m.par["repr.int.Site"]  # linear predictor (linear.p = intercept + B(initial size) + plus site level intercept)
				p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale)	
				# p = response (the actual probability of flowering)
				return(p)
		}	
#----------------------------------------
	#  RECRUITMENT SIZE DISTRIBUTION FUNCTION: 
		RecrDist <- function(z1,z,m.par) {	
		#----------------------------------------
		## 1B.4 - FECUNDITY function (how many fruits produced?):			
				# Number of fruits given individual flowered: 		
				count_number= exp((m.par["surv.int"]) + (m.par["FecCountComp.int.site"]) +
									(m.par["surv.z"]*z) + (m.par["FecCountComp.slope.site"]*z*1))					

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
	# Run a loop through all sites & separate IPM model for each site 
	siteLambdas <- data.frame(siteLambdas)
	head(siteLambdas)
	siteLambdas <- as.matrix(siteLambdas)
	testObject <- function(object)
			{exists(as.character(substitute(object)))}

	#
	##
	###
	####
	#####
	######
	#######
for(j in 1:5000){	
	for(p in 1:length(all_sites)){
		# test if site replicate was even successfull?
		# also are there any missing values in m.par?
		moot = exists(paste0("m.par.", all_sites[p], "_", j))
	
	if(exists(paste0("m.par.", all_sites[p], "_", j))) {
					# if everythings ok then, 
					m.par <- get(paste0("m.par.", all_sites[p], "_", j)) # set parameter list according to site
					if('TRUE' %in%  is.na(get(paste0("m.par.", all_sites[p], "_", j)))){
								siteLambdas[j,all_sites[p]] <- NA 
						} else {
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
					G = array(0,dim=c(n.size,n.size))
					## runs Grow_z1z for each combination of possible sizes
					G = h*outer(y,y,Grow_z1z,m.par=m.par)
					# the function 'outer()' evaluates the matrix at all pairwise combinations of the two vectors y & y and returns matrices representing the kernel components of the growth 
					# make Vectors from starting size: survive, prob repro, number rec
					S = SurvFrom_z(y,m.par=m.par) # survival vector 
					ProbRepVec= Prob_repr(y,m.par=m.par) # vector of probability of flowering 
					###
					RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
						#dim(RecMatrix); class(RecMatrix)
				
				# Growth matrix 
					P=G 						     # placeholder; redefine P on the next line
					for(i in 1:n) P[,i]=G[,i]*S[i]	 # growth/survival matrix
				# Fecundity matrix 
					F=RecMatrix
					for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix
						# this is equivalent to adding "*Prob_repr(z,m.par)" to the last line of the RecrDist function 
						# recall in the fecundity model we only included observations from individuals that were already flowering 
					K=P+F # full matrix, add matrices together 
				# basic lambda estimate	
					lam <- Re(eigen(K)$values[1]) #  0.939
					siteLambdas[j,all_sites[p]] <- lam 
				}	
				} else {
				siteLambdas[j,all_sites[p]] <- NA 
	
				}	
	
	}
	(print(j))
	}

	
head(siteLambdas, 15)
tail(siteLambdas, 15)
siteLambdas <- data.frame(siteLambdas)
slamsTrim <- siteLambdas
slamsTrim$CALAPOOIA[slamsTrim$CALAPOOIA > 15] <- NA
slamsTrim$COAST[slamsTrim$COAST > 15] <- NA
slamsTrim$HUNTER[slamsTrim$HUNTER > 15] <- NA
slamsTrim$LOOK[slamsTrim$LOOK > 15] <- NA
slamsTrim$MOSBY[slamsTrim$MOSBY > 15] <- NA
slamsTrim$ROCK[slamsTrim$ROCK > 15] <- NA
slamsTrim$THOMAS[slamsTrim$THOMAS > 15] <- NA
slamsTrim$WILEY[slamsTrim$WILEY > 15] <- NA
summary(slamsTrim)

#siteLambdas
dev.off()

setwd(path.obj)
write.csv(siteLambdas, file="siteLambdasBootCI_sites.csv")




#
##
###
####
######
#######
########
#########
##########
###########
###########
############
############# # RUN BOOTSTRAP REPLICATES FOR PLOTS
############
###########
###########
##########
#########
########
#######
######
#####
####
###
##
#


	
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
	stdev.lnrom.RecrSize =  NA
	)	

	m.parN <- m.par # will duplicate this list & as with the demography dataset with use recruitment estimates from 
	# demography sites 
	#m.parC <- m.par; m.parS <- m.par # won't end up using m.parS (South) or m.parC (Center) except for exploratory purposes. 

########################################################################################
	# LOAD BOOTSRAP REPLICATE COEFFICENTS
	setwd(path.obj)
	grow_boot <- get(load('GrowModBootReplicates500_plots.rda'))
	surv_boot <- get(load('SurvModBootReplicates500_plots.rda'))
	repr_boot <- get(load('ReprModBootReplicates500_plots.rda'))
########################################################################################	
	# MAKE MATRIX FOR LAMBDA VALUES 
	to_loop <- levels(dfclean$SiteID)[2:length(levels(dfclean$SiteID))] #all sites but the one used as a dummy variable. 
	to_loop
	all_sites <-  c("CALAPOOIA", to_loop) # add on calapooia as a site 
	siteLambdas <- matrix(NA, nrow = length(repr_boot), ncol = 9)
	colnames(siteLambdas) <- c("REP", "CALAPOOIA", "COAST", "HUNTER", "LOOK", "MOSBY", "ROCK", "THOMAS", "WILEY")
	siteLambdas[,1] <- c(seq(1:5000))
	head(siteLambdas)
#####################################################################################
	# load in bootstrap replicates 
	setwd(path.obj)
	BootStrapData <- read.csv(file="BootReplicateDataSets_Plots.csv")
########################################################################################

for(u in 1:length(surv_boot)){
		thisRep <- BootStrapData[which(BootStrapData$Replicate==u),]
		SitesInRep <- unique(as.factor(thisRep$SiteID))
		SitesInRep <- sort(as.character(SitesInRep)) # make alphabetical to match dummy variable 
		to_loop <- SitesInRep[2:length(SitesInRep)]
		
#----------------------------------------
## 1A.2 - growth (continuous variable 'z1')
	# Growth (given survival), conditional on survival 
	#source(".R") # only run this once on the first time (takes a while)
	setwd(path.code)
	#source("04.2.2_Grow.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	grow_model <- get(load("GrowMod.rda")) # load in top growth mixed model 

#----------------------------------------	
## 1A.5 - recruitment model (probability that a fruit or seed dropped from a plant 
	# becomes a new recruit the next year 
	setwd(path.code)
	#source("03_RecruitmentModelComp.R") # only run this once on the first time 
	setwd(path.obj)
	Rec_model <- read.csv(file="RecruitmentModCoefs_plots.csv")
	# draw from distribution and assign values 
		m.parN["SeedPerFruit"]<- Rec_model$meanSeedcount[2]
		m.parN["SeedPerFruit"]<- rnorm(1, mean = Rec_model$meanSeedcount[2], sd = Rec_model$Seedsd[2])
		# value must be positive
		m.parN["SeedPerFruit"] <- ifelse(m.parN["SeedPerFruit"] <= 0, Rec_model$meanSeedcount[2], m.parN["SeedPerFruit"])
		# slope for recruitment coefficient 
		m.parN["RecrCoef.slope"] <- runif(1, Rec_model$LL[2], Rec_model$UL[2])
		assign(paste0("m.parN_", u), m.parN)	

	
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
		m.parN["mean.lnrom.RecrSize"]<- newRec_dist$meanlog[2]
		m.parN["stdev.lnrom.RecrSize"]<- newRec_dist$sdlog[2]
		assign(paste0("m.parN_", u), m.parN)	

#----------------------------------------
## 1A.4 - Fecundity (Given reproductive) (count variable, number of fruits 'Fec')
	# How many fruits were produced (given plant was reproductive 0, 1, 2 .... to max fruits). 
	setwd(path.code)
	#source("04.2.4_Fecundity.R") # only run this once on the first time (takes a while)
	setwd(path.obj)
	Fec_model <- get(load("FecMod_NB.rda")) # load in top reproductive glm mixed model 
	# don't swap estimated values as we did above. 
	fixef(Fec_model)
	library(coefplot2)
	m.par.est_Fec <- fixef(Fec_model)
	names(m.par.est_Fec) <- paste0("Fec_", names(m.par.est_Fec))	
		
 	SitesInRep
	thisSurv <- surv_boot[[u]]
   	thisGrow <- grow_boot[[u]]
   	thisRepr <- repr_boot[[u]]

	
	# Fill m.par lists with site values.
	# will have to add values for "DUMMY" first ~ DUMMY, BUT CHANGE NAME 
	m.par.DUMMY <- m.parN # load missing values from northern region 
	m.par.DUMMY["surv.int"] =  thisSurv["(Intercept)"]
	m.par.DUMMY["surv.z"] = thisSurv["z"]
	m.par.DUMMY["surv.int.Site"] = 0
	m.par.DUMMY["grow.int"] = thisGrow["(Intercept)"]
	m.par.DUMMY["grow.z"] = thisGrow["z"]
	m.par.DUMMY["grow.int.Site"] = 0
	m.par.DUMMY["grow.z.Site"] = 0
	m.par.DUMMY["grow.sd"] = sigma(grow_model)
	m.par.DUMMY["repr.int"] = thisRepr["(Intercept)"]
	m.par.DUMMY["repr.slope"] = thisRepr["z"]
	m.par.DUMMY["repr.int.Site"] = 0
		if(SitesInRep[1]=="CALAPOOIA") {
			m.par.DUMMY["FecCountComp.slope.site"] = 0
			m.par.DUMMY["FecCountComp.int.site"] = 0
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				} else {
			m.par.DUMMY["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", SitesInRep[1])]
			m.par.DUMMY["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
			m.par.DUMMY["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
				}
	m.par.DUMMY
	assign(paste0("m.par.", SitesInRep[1], "_", u), m.par.DUMMY)
	
	# loop to fill m.par values for sites 
	for(i in 1:length(to_loop)){
		temp <- m.parN # temporary obj to fill 
		m.par.temp <- m.parN # load missing values from northern region 
		m.par.temp["surv.int"] =  thisSurv["(Intercept)"]
		m.par.temp["surv.z"] = thisSurv["z"]
		m.par.temp["surv.int.Site"] = thisSurv[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.int"] = thisGrow["(Intercept)"]
		m.par.temp["grow.z"] = thisGrow["z"]
		m.par.temp["grow.int.Site"] = thisGrow[paste0("SiteID", to_loop[i])]
		m.par.temp["grow.z.Site"] = thisGrow[paste0("z:SiteID", to_loop[i])]
		m.par.temp["grow.sd"] = sigma(grow_model)
		m.par.temp["repr.int"] =thisRepr["(Intercept)"]
		m.par.temp["repr.slope"] = thisRepr["z"]
		m.par.temp["repr.int.Site"] = thisRepr[paste0("SiteID", to_loop[i])]
		m.par.temp["FecCountComp.slope.site"] = m.par.est_Fec[paste0("Fec_z:SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int.site"] = m.par.est_Fec[paste0("Fec_SiteID", to_loop[i])]
		m.par.temp["FecCountComp.int"] = m.par.est_Fec["Fec_(Intercept)"]
		m.par.temp["FecCountComp.slope"] = m.par.est_Fec["Fec_z"]
		# m.par.temp
		assign(paste0("m.par.", to_loop[i], "_", u), m.par.temp)
	}	
	(print(u))
	}

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1C - Define basic life history functions for a manual IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------
	## 1B.1 - SURVIVAL function:  
			SurvFrom_z <- function(z, m.par)
		{
			linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  +  m.par["surv.int.Site"]   # linear predictor (linear.p = intercept + B(initial size))
			p <- 1/(1+exp(-linear.p))                                 						 # logistic transformation to probability (p - logit scale) 
			return(p)																	     # returns p on logit scale 
		}		
#----------------------------------------
	## 1B.2 - GROWTH function:
			Grow_z1z <- function(z1, z, m.par)
		{
			mu <- m.par["grow.int"] + m.par["grow.z"] * z  +  m.par["grow.int.Site"] +  m.par["grow.z.Site"]*z*1       # mean size next year (mu = intercept + B(initial size))
			sig <- m.par["grow.sd"]                                 # sd about mean (sig = stdev of growth)
			p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # probability density fxn (pdf) that you are size z1 given you were size z
			return(p.den.grow)										# returns probability density function of growth 
		}	
#----------------------------------------
	## 1B.3 - FLOWERING function (becomes reproductive?):	
		Prob_repr <-function(z,m.par) {
			# Probability of flowering based on size (0 to 1)
				z2 <- z*z # for quadratic term 
				linear.p <- m.par["repr.int"] + (m.par["repr.slope"] * z) + m.par["repr.int.Site"]  # linear predictor (linear.p = intercept + B(initial size) + plus site level intercept)
				p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability (p - logit scale)	
				# p = response (the actual probability of flowering)
				return(p)
		}	
#----------------------------------------
	#  RECRUITMENT SIZE DISTRIBUTION FUNCTION: 
		RecrDist <- function(z1,z,m.par) {	
		#----------------------------------------
		## 1B.4 - FECUNDITY function (how many fruits produced?):			
				# Number of fruits given individual flowered: 		
				count_number= exp((m.par["surv.int"]) + (m.par["FecCountComp.int.site"]) +
									(m.par["surv.z"]*z) + (m.par["FecCountComp.slope.site"]*z*1))					

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
	# Run a loop through all sites & separate IPM model for each site 
	siteLambdas <- data.frame(siteLambdas)
	head(siteLambdas)
	siteLambdas <- as.matrix(siteLambdas)
	testObject <- function(object)
			{exists(as.character(substitute(object)))}

	#
	##
	###
	####
	#####
	######
	#######
for(j in 1:5000){	
	for(p in 1:length(all_sites)){
		# test if site replicate was even successfull?
		# also are there any missing values in m.par?
		moot = exists(paste0("m.par.", all_sites[p], "_", j))
	
	if(exists(paste0("m.par.", all_sites[p], "_", j))) {
					# if everythings ok then, 
					m.par <- get(paste0("m.par.", all_sites[p], "_", j)) # set parameter list according to site
					if('TRUE' %in%  is.na(get(paste0("m.par.", all_sites[p], "_", j)))){
								siteLambdas[j,all_sites[p]] <- NA 
						} else {
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
					G = array(0,dim=c(n.size,n.size))
					## runs Grow_z1z for each combination of possible sizes
					G = h*outer(y,y,Grow_z1z,m.par=m.par)
					# the function 'outer()' evaluates the matrix at all pairwise combinations of the two vectors y & y and returns matrices representing the kernel components of the growth 
					# make Vectors from starting size: survive, prob repro, number rec
					S = SurvFrom_z(y,m.par=m.par) # survival vector 
					ProbRepVec= Prob_repr(y,m.par=m.par) # vector of probability of flowering 
					###
					RecMatrix=h*outer(y,y,RecrDist,m.par=m.par) # reproduction matrix
						#dim(RecMatrix); class(RecMatrix)
				
				# Growth matrix 
					P=G 						     # placeholder; redefine P on the next line
					for(i in 1:n) P[,i]=G[,i]*S[i]	 # growth/survival matrix
				# Fecundity matrix 
					F=RecMatrix
					for(i in 1:n) F[,i]=RecMatrix[,i]*ProbRepVec[i]	 # flower/fecundity matrix
						# this is equivalent to adding "*Prob_repr(z,m.par)" to the last line of the RecrDist function 
						# recall in the fecundity model we only included observations from individuals that were already flowering 
					K=P+F # full matrix, add matrices together 
				# basic lambda estimate	
					lam <- Re(eigen(K)$values[1]) #  0.939
					siteLambdas[j,all_sites[p]] <- lam 
				}	
				} else {
				siteLambdas[j,all_sites[p]] <- NA 
	
				}	
	
	}
	(print(j))
	}

	
head(siteLambdas, 15)
tail(siteLambdas, 15)
siteLambdas <- data.frame(siteLambdas)
slamsTrim <- siteLambdas
slamsTrim$CALAPOOIA[slamsTrim$CALAPOOIA > 15] <- NA
slamsTrim$COAST[slamsTrim$COAST > 15] <- NA
slamsTrim$HUNTER[slamsTrim$HUNTER > 15] <- NA
slamsTrim$LOOK[slamsTrim$LOOK > 15] <- NA
slamsTrim$MOSBY[slamsTrim$MOSBY > 15] <- NA
slamsTrim$ROCK[slamsTrim$ROCK > 15] <- NA
slamsTrim$THOMAS[slamsTrim$THOMAS > 15] <- NA
slamsTrim$WILEY[slamsTrim$WILEY > 15] <- NA
summary(slamsTrim)

#siteLambdas
dev.off()

setwd(path.obj)
write.csv(siteLambdas, file="siteLambdasBootCI_plots.csv")





