# 03.2_Make lambda plots.R
############################################################################################
##
##      Make basic lambda plots 
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
library(ReporteRs)

set.seed(270875)
	def <- par("mai")
	par(ps=10) # default text size is 10 

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

	# load in bootstrap replicate lambda value predictions:
	setwd(path.obj)
	siteLambdas <- read.csv(file="siteLambdasBootCI_plots.csv")
	MainLams <- read.csv(file="siteLambdas.csv")
	MainLams2 <- read.csv(file="siteLambdas_MoistScore0_noInt.csv") # set moisture to 0
	ENMs <- read.csv(file="site_predsOLD.csv")
	ENMs2 <- read.csv(file="site_predsANNUAL.csv")

	#summary(siteLambdas)
	head(ENMs, 1)
	MainLams
	
	## Add an alpha value to a colour
						add.alpha <- function(col, alpha=1){
						  if(missing(col))
							stop("Please provide a vector of colours.")
						  apply(sapply(col, col2rgb)/255, 2, 
											 function(x) 
											   rgb(x[1], x[2], x[3], alpha=alpha))  
						}


#--------------------------------------------------------------------------------
# LOAD, COMBINE AND STRUCTURE DATA
#--------------------------------------------------------------------------------

#############################################################
# make plots with CI's

library(vioplot)
setwd(path.funct); source("my_vioplot.R")
	
	siteLambdas <- siteLambdas[,c(3:ncol(siteLambdas))]
	# >>>>> need to revisit 
	siteLambdas[siteLambdas > 7] <- NA # remove errorenous bootstrap replications ? 
	
# get 95% quantiles, basic 
	calap <- quantile(siteLambdas$CALAPOOIA, c(0.025, 0.975), na.rm=TRUE) 
	coast <- quantile(siteLambdas$COAST, c(0.025, 0.975), na.rm=TRUE) 
	huntr <- quantile(siteLambdas$HUNTER, c(0.025, 0.975), na.rm=TRUE)
	look <- quantile(siteLambdas$LOOK, c(0.025, 0.975), na.rm=TRUE)
	mosby <- quantile(siteLambdas$MOSBY, c(0.025, 0.975), na.rm=TRUE)
	rock <- quantile(siteLambdas$ROCK, c(0.025, 0.975), na.rm=TRUE)
	thomas <- quantile(siteLambdas$THOMAS, c(0.025, 0.975), na.rm=TRUE)
	wiley <- quantile(siteLambdas$WILEY, c(0.025, 0.975), na.rm=TRUE)

	
	
	
	
	setwd(path.funct)
	source("bootBCa.R")
	
	
	setwd(path.obj)
	siteLambdas <- read.csv(file="siteLambdasBootCI_plots.csv")
	
	tt <- siteLambdas[complete.cases(siteLambdas$THOMAS),c("THOMAS")]
	test_plot <- BCa(tt, NA, mean, alpha=c(0.025,0.975), M=length(tt))
	
	setwd(path.obj)
	siteLambdas <- read.csv(file="siteLambdasBootCI_sites.csv")
	
	tt <- siteLambdas[complete.cases(siteLambdas$THOMAS),c("THOMAS")]
	test_site <- BCa(tt, NA, mean, alpha=c(0.025,0.975))

	setwd(path.obj)
	siteLambdas <- read.csv(file="siteLambdasBootCI.csv")
	
	tt <- siteLambdas[complete.cases(siteLambdas$THOMAS),c("THOMAS")]
	test_raw <- BCa(tt, NA, mean, alpha=c(0.025,0.975))
	
	
	test_plot
	test_site
	test_raw
	
	
	
	
	
	
# load in main estimates for sites
setwd(path.obj)
	MainLams <- MainLams[,c("V1", "V2")]
	MainLams$LL <- c(calap[1], coast[1], huntr[1], look[1], mosby[1], rock[1], thomas[1], wiley[1])
	MainLams$UL <- c(calap[2], coast[2], huntr[2], look[2], mosby[2], rock[2], thomas[2], wiley[2])

	MainLams2 <- MainLams2[,c("V1", "V2")]
	MainLams2$LL <- c(calap[1], coast[1], huntr[1], look[1], mosby[1], rock[1], thomas[1], wiley[1])
	MainLams2$UL <- c(calap[2], coast[2], huntr[2], look[2], mosby[2], rock[2], thomas[2], wiley[2])

	
	par(mfrow=c(1,2))
	barplot(MainLams$V2, 
			names=MainLams$V1,
			ylim=c(0,max(MainLams$V2)), ylab="lambda", main="bootstrapped estimates")
	barplot(MainLams2$V2, 
			names=MainLams$V1,
			ylim=c(0,max(MainLams$V2)), ylab="lambda", main="bootstrapped estimates")
			
	setwd(path.obj)
#load in latitude values 
	enm_el <- ENMs
	ENMs <- ENMs[which(ENMs$ID1=='trans'), ]
	levels(ENMs$ID2)[levels(ENMs$ID2) == "Rock Creek" ] <- "ROCK"
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Coast Fork Willamette'] <- 'COAST'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Mosby Creek'] <- 'MOSBY'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Wiley Creek'] <- 'WILEY'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Thomas'] <- 'THOMAS'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Hunter'] <- 'HUNTER'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Looking Glass'] <- 'LOOK'
	levels(ENMs$ID2)[levels(ENMs$ID2) == 'Calapooia Creek'] <- 'CALAPOOIA'
	binder <- merge(MainLams, ENMs, by.x="V1", by.y="ID2", all.x=TRUE, all.y=FALSE)

	ENMs2$ID1 <- c('trans', 'trans', 'trans', 'trans', 'trans', 'trans', 'trans', 'trans', 'trans', 'elev', 'elev','elev','elev')
	enm_el2 <- ENMs2
	ENMs2 <- ENMs2[which(ENMs2$ID1=='trans'), ]
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == "Rock Creek" ] <- "ROCK"
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Coast Fork Willamette'] <- 'COAST'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Mosby Creek'] <- 'MOSBY'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Wiley Creek'] <- 'WILEY'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Thomas'] <- 'THOMAS'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Hunter'] <- 'HUNTER'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Looking Glass'] <- 'LOOK'
	levels(ENMs2$ID2)[levels(ENMs2$ID2) == 'Calapooia Creek'] <- 'CALAPOOIA'
	binder2 <- merge(MainLams2, ENMs2, by.x="V1", by.y="ID2", all.x=TRUE, all.y=FALSE)

	
	drops <- c("Tmax01", "Tmax02", "Tmax03", "Tmax04",
			"Tmax05", "Tmax06", "Tmax07", "Tmax08", "Tmax09", "Tmax10", 
			"Tmax11", "Tmax12", "Tmin01", "Tmin02", "Tmin03", "Tmin04",
			"Tmin05", "Tmin06", "Tmin07", "Tmin08", "Tmin09", "Tmin10",
			"Tave05", "Tave06", "Tave07", "Tave08", "Tave09", "Tave10",
			"Tave11", "Tave12", "PPT01", "PPT02", "PPT03", "PPT04",
			"PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT10",
			"PPT11", "PPT12", "Tmin11", "Tmin12", "Tave01", "Tave02", "Tave03", "Tave04")
	binder <- binder[,!(names(binder) %in% drops)]
	names(binder)[names(binder)=="V1"] <- "name"
	names(binder)[names(binder)=="V2"] <- "lam"
	
	binder2 <- binder2[,!(names(binder2) %in% drops)]
	names(binder2)[names(binder2)=="V1"] <- "name"
	names(binder2)[names(binder2)=="V2"] <- "lam"


###############################################################
# LOAD IN SITE LEVEL VITAL RATE DATA
setwd(path.obj)
	SurvMod <- read.csv(file="SurvModBootCoef_plots.csv")
	GrowMod <- read.csv(file="GrowModBootCoef_plots.csv")
	ReprMod <- read.csv(file="ReprModBootCoef_plots.csv")
	Fec_model <- get(load("FecMod_NB.rda")) # load in top reproductive glm mixed model 
		# don't swap estimated values as we did above. 
		fef <- fixef(Fec_model)
		library(coefplot2)
		m.par.est_Fec <- fixef(Fec_model)
	
# FOR REPRODUCTION...
		ReprMod$Repr <- ReprMod$BootMean[1] + ReprMod$BootMean	
		ReprMod$Repr[1] <- ReprMod$BootMean[1]	
		ReprMod$Repr <- (mean(ReprMod$Repr)/ReprMod$Repr) - 1
		# lower
		ReprMod$ReprL <- ReprMod$X2.5.[1] + ReprMod$X2.5.
		ReprMod$ReprL[1] <- ReprMod$X2.5.[1]	
		ReprMod$ReprL <- mean(ReprMod$ReprL)/ReprMod$ReprL
		# upper
		ReprMod$ReprU <- ReprMod$X97.5.[1] + ReprMod$X97.5.
		ReprMod$ReprU[1] <- ReprMod$X97.5.[1]	
		ReprMod$ReprU <- mean(ReprMod$ReprU)/ReprMod$ReprU	
		ReprMod$name <- c("CALAPOOIA", "NA", "COAST", "HUNTER", "LOOK", "MOSBY", "ROCK", "THOMAS", "WILEY", "NA")
		ReprMod <- ReprMod[,c("name", "Repr", "ReprL", "ReprU")]
		binder <- merge(binder, ReprMod, by.x="name", by.y="name", all.x=TRUE, all.y=FALSE)

# FOR SURVIVAL...
		SurvMod$Surv <- SurvMod$BootMean[1] + SurvMod$BootMean	
		SurvMod$Surv[1] <- SurvMod$BootMean[1]	
		meany <- mean(SurvMod$Surv[c(1,3:9)])
		SurvMod$Surv <- (meany/SurvMod$Surv) - 1
		# lower
		SurvMod$SurvL <- SurvMod$X2.5.[1] + SurvMod$X2.5.
		SurvMod$SurvL[1] <- SurvMod$X2.5.[1]	
		SurvMod$SurvL <- meany/SurvMod$SurvL
		# upper
		SurvMod$SurvU <- SurvMod$X97.5.[1] + SurvMod$X97.5.
		SurvMod$SurvU[1] <- SurvMod$X97.5.[1]	
		SurvMod$SurvU <- meany/SurvMod$SurvU	
		SurvMod$name <- c("CALAPOOIA", "NA", "COAST", "HUNTER", "LOOK", "MOSBY", "ROCK", "THOMAS", "WILEY", "NA")
		SurvMod <- SurvMod[,c("name", "Surv", "SurvL", "SurvU")]
		binder <- merge(binder, SurvMod, by.x="name", by.y="name", all.x=TRUE, all.y=FALSE)

# FOR GROWTH...
		# used mean size of z - 3.853 & predict values for sites 
		g <- matrix(NA, nrow=8, ncol=4); g <- data.frame(g)
		g[,1] <- binder$name
		colnames(g) <- c("name", "Grow", "GrowL", "GrowU")	
		for(i in 1:3){
		g[1,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) 
		g[2,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[3,(i+3)]) + (3.853*(GrowMod[10,(i+3)]))
		g[3,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[4,(i+3)]) + (3.853*(GrowMod[11,(i+3)]))
		g[4,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[5,(i+3)]) + (3.853*(GrowMod[12,(i+3)]))
		g[5,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[6,(i+3)]) + (3.853*(GrowMod[13,(i+3)]))
		g[6,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[7,(i+3)]) + (3.853*(GrowMod[14,(i+3)]))
		g[7,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[8,(i+3)]) + (3.853*(GrowMod[15,(i+3)]))
		g[8,i+1] <- GrowMod[1,(i+3)] + (3.853*GrowMod[2,(i+3)]) + (GrowMod[9,(i+3)]) + (3.853*(GrowMod[16,(i+3)]))
		g[,i+1] <- g[,i+1]/6.719192  # mean(g[,2])
		}
		g$Grow <- g$Grow - 1
		binder <- merge(binder, g, by.x="name", by.y="name", all.x=TRUE, all.y=FALSE)

# FOR FECUNDITY 		
		# used mean size of z - 3.853 & predict values for sites 
		g <- matrix(NA, nrow=8, ncol=4); g <- data.frame(g)
		g[,1] <- binder$name
		colnames(g) <- c("name", "Fec", "FecL", "FecU")	
		#for(i in 1:3){
		i=1
		g[1,i+1] <- fef[1] + (3.853*fef[2]) 
		g[2,i+1] <- fef[1] + (3.853*fef[2]) + fef[3] + (3.853*fef[10])
		g[3,i+1] <- fef[1] + (3.853*fef[2]) + fef[4] + (3.853*fef[11])
		g[4,i+1] <- fef[1] + (3.853*fef[2]) + fef[5] + (3.853*fef[12])
		g[5,i+1] <- fef[1] + (3.853*fef[2]) + fef[6] + (3.853*fef[13])
		g[6,i+1] <- fef[1] + (3.853*fef[2]) + fef[7] + (3.853*fef[14])
		g[7,i+1] <- fef[1] + (3.853*fef[2]) + fef[8] + (3.853*fef[15])
		g[8,i+1] <- fef[1] + (3.853*fef[2]) + fef[9] + (3.853*fef[16])
		g[,i+1] <- g[,i+1]/mean(g[,2])  # mean(g[,2])
		#}
		g$Fec <- g$Fec - 1
		binder <- merge(binder, g, by.x="name", by.y="name", all.x=TRUE, all.y=FALSE)
	
# mean ENM 
		for(i in 1:dim(binder)[1]){
		binder$EN[i] <- mean(c(binder$LRavg[i], binder$GAMavg[i], binder$RFavg[i], binder$BRTavg[i], binder$MAXavg[i]))
		binder$EN_L[i] <- min(c(binder$LRavg[i], binder$GAMavg[i], binder$RFavg[i], binder$BRTavg[i], binder$MAXavg[i]))
		binder$EN_U[i] <- max(c(binder$LRavg[i], binder$GAMavg[i], binder$RFavg[i], binder$BRTavg[i], binder$MAXavg[i]))
		}
		
		for(i in 1:dim(binder2)[1]){
		binder2$EN[i] <- mean(c(binder2$LRavg[i], binder2$GAMavg[i], binder2$RFavg[i], binder2$BRTavg[i], binder2$MAXavg[i]))
		binder2$EN_L[i] <- min(c(binder2$LRavg[i], binder2$GAMavg[i], binder2$RFavg[i], binder2$BRTavg[i], binder2$MAXavg[i]))
		binder2$EN_U[i] <- max(c(binder2$LRavg[i], binder2$GAMavg[i], binder2$RFavg[i], binder2$BRTavg[i], binder2$MAXavg[i]))
		}
		
		
# QUICK CHECK 	
		par(mfrow=c(2,2))
		
		plot(binder$Latitude, binder$Surv, pch=19, col="red", ylim=c(-1,1))	
		points(binder$Latitude, binder$Grow, pch=19, col="blue", ylim=c(-0.5,0.5))	
		points(binder$Latitude, binder$Repr, pch=19, col="green", ylim=c(-0.5,0.5))	
		points(binder$Latitude, (binder$lam-1), pch=21, bg="yellow", ylim=c(-0.5,0.5))	
		points(binder$Latitude, binder$Fec, pch=21, bg="purple", ylim=c(-0.5,0.5))	
		
		plot(binder$EN, binder$Surv, pch=19, col="red", ylim=c(-1,1))	
		points(binder$EN, binder$Grow, pch=19, col="blue", ylim=c(-0.5,0.5))	
		points(binder$EN, binder$Repr, pch=19, col="green", ylim=c(-0.5,0.5))	
		points(binder$EN, (binder$lam-1), pch=21, bg="yellow", ylim=c(-0.5,0.5))	
		points(binder$EN, binder$Fec, pch=21, bg="purple", ylim=c(-0.5,0.5))	
		
		plot(binder2$Latitude, binder2$lam, pch=19, col="red", ylim=c(-1,1))	
		plot(binder2$EN, binder2$lam, pch=19, col="red", ylim=c(-1,1))	

		
# Prepare for ELEVATION
	enm_el <- enm_el[,c("ID1", "ID2", "Latitude", "Longitude", "Elevation", "LRavg", "GAMavg", "RFavg", "BRTavg", "MAXavg")]
	enm_el <- enm_el[which(enm_el$ID1=='elev'), ]
		for(i in 1:dim(enm_el)[1]){
		enm_el$EN[i] <- mean(c(enm_el$LRavg[i], enm_el$GAMavg[i], enm_el$RFavg[i], enm_el$BRTavg[i], enm_el$MAXavg[i]))
		enm_el$EN_L[i] <- min(c(enm_el$LRavg[i], enm_el$GAMavg[i], enm_el$RFavg[i], enm_el$BRTavg[i], enm_el$MAXavg[i]))
		enm_el$EN_U[i] <- max(c(enm_el$LRavg[i], enm_el$GAMavg[i], enm_el$RFavg[i], enm_el$BRTavg[i], enm_el$MAXavg[i]))
		}
	el <- read.csv(file="Old_Amy_elev_trans.csv")
	enm_el <- merge(enm_el, el, by.x="ID2", by.y="SITE")
	
# NORMALIZE VARIABLES FOR RELATIVE FITNESS METRIC 
	stdGROW <- (enm_el$GROW - mean(enm_el$GROW))/sd(enm_el$GROW)
  	stdFLOW <- (enm_el$Flower - mean(enm_el$Flower))/sd(enm_el$Flower)
	enm_el$fit <- ((stdGROW + stdFLOW)/2)
	# add 95% CIs from SE
	flrSE_UL <- enm_el$Flower + enm_el$flrSE 
	flrSE_LL <- enm_el$Flower - enm_el$flrSE  
  	flrSE_LL <- (flrSE_LL - mean(enm_el$Flower))/sd(enm_el$Flower)
  	flrSE_UL <- (flrSE_UL - mean(enm_el$Flower))/sd(enm_el$Flower)
	grSE_UL <- enm_el$GROW + enm_el$grSE 
	grSE_LL <- enm_el$GROW - enm_el$grSE  
  	grSE_LL <- (grSE_LL - mean(enm_el$GROW))/sd(enm_el$GROW)
  	grSE_UL <- (grSE_UL - mean(enm_el$GROW))/sd(enm_el$GROW)
	enm_el$fitUL <- ((grSE_UL + flrSE_UL)/2)
	enm_el$fitLL <- ((grSE_LL + flrSE_LL)/2)

	
	enm_el2 <- enm_el2[,c("ID1", "ID2", "Latitude", "Longitude", "Elevation", "LRavg", "GAMavg", "RFavg", "BRTavg", "MAXavg")]
	enm_el2 <- enm_el2[which(enm_el2$ID1=='elev'), ]
		for(i in 1:dim(enm_el2)[1]){
		enm_el2$EN[i] <- mean(c(enm_el2$LRavg[i], enm_el2$GAMavg[i], enm_el2$RFavg[i], enm_el2$BRTavg[i], enm_el2$MAXavg[i]))
		enm_el2$EN_L[i] <- min(c(enm_el2$LRavg[i], enm_el2$GAMavg[i], enm_el2$RFavg[i], enm_el2$BRTavg[i], enm_el2$MAXavg[i]))
		enm_el2$EN_U[i] <- max(c(enm_el2$LRavg[i], enm_el2$GAMavg[i], enm_el2$RFavg[i], enm_el2$BRTavg[i], enm_el2$MAXavg[i]))
		}
	el <- read.csv(file="Old_Amy_elev_trans.csv")
	enm_el2 <- merge(enm_el2, el, by.x="ID2", by.y="SITE")
	jusFit <- enm_el[,c("ID2", "fit", "fitUL", "fitLL")]
	enm_el2 <- merge(enm_el2, jusFit, by.x="ID2", by.y="ID2")

	
	
	
	#check 
	par(mfrow=c(1,3))
	plot(enm_el2$Elevation, enm_el2$fit, pch=19)
	points(enm_el2$Elevation, enm_el2$fitUL, pch=19, col="red")
	points(enm_el2$Elevation, enm_el2$fitLL, pch=19, col="blue")
  	plot(enm_el$EN, enm_el$fit, pch=19, main="long term")
	text(enm_el$EN, enm_el$fit, labels=paste0(enm_el$ID2, " ", enm_el$ELEV))
  	plot(enm_el2$EN, enm_el2$fit, pch=19, main="short term")
	points(enm_el2$EN, enm_el2$fitUL, pch=19, col="red")
	points(enm_el2$EN, enm_el2$fitLL, pch=19, col="blue")
  	text(enm_el2$EN, enm_el2$fit, labels=paste0(enm_el2$ID2, " ", enm_el2$Elevation))

#--------------------------------------------------------------------------------
# MAKE PLOTS
#--------------------------------------------------------------------------------
	
############################################################
############################################################
# Make a within & beyond barplot with CI's 
	setwd(path.fig)

	binder <- binder[order(binder$Latitude),]
	library(Hmisc)
	library(gplots)
	BarCols = c("coral2", "coral2", "coral2", "coral2", "blue4", "blue4", "blue4", "blue4") 
	BarColsB <- c("coral2", "coral2", "blue4", "blue4")
	## Add an alpha value to a colour
	BarCols2 <- as.character(add.alpha(BarCols, alpha=0.6)) 
	BarColsB2 <- as.character(add.alpha(BarColsB, alpha=0.6))

	
# VIOLIN PLOT 	~ site lambdas
	my_vioplot(siteLambdas$LOOK[complete.cases(siteLambdas$LOOK)],
			siteLambdas$ROCK[complete.cases(siteLambdas$ROCK)],
				siteLambdas$COAST[complete.cases(siteLambdas$COAST)],
				siteLambdas$MOSBY[complete.cases(siteLambdas$MOSBY)],
				siteLambdas$CALAPOOIA[complete.cases(siteLambdas$CALAPOOIA)],
				siteLambdas$WILEY[complete.cases(siteLambdas$WILEY)],
				siteLambdas$THOMAS[complete.cases(siteLambdas$THOMAS)],
				siteLambdas$HUNTER[complete.cases(siteLambdas$HUNTER)],
				col=BarCols, names=binder$name)
				abline(h=1, col="darkgrey", lty=2)

# binder - NRL (long term climate & no moisture)
# binder2 - NRL (short term climate + moisture)
	binder2 <- binder2[order(binder2$Latitude),]
	binder[,c(1:3)]
	binder2[,c(1:3)]
	# use this value to include moisture
	binder$lamMoist <- binder2$lam		

dev.off()
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# PLOT LAMBDAS (latitude vs. fitness)

setwd(path.fig)	
pdf(file = "03.4_IPM_finals_moisture_lat_el_lam2.pdf", width=(8.5 - (1.25 + 0.75)), height=(5 - (0.75 + 0.75)), family="Times")
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/


# PLOT LAMBDAS (latitude vs. fitness)
	par(mfrow=c(1,2))
		plot(binder$Latitude, binder$lam, ylab =expression("Population growth rate,  "*lambda), pch=21, cex=1,
		bg=BarCols, ylim=c(0, max(binder$lam*1.1)), 
		xlab=expression(paste("Latitude ", ""^o*"N", sep="")), main="Northern Range Limit", cex.main=1)
		abline(h=1, col="pink", lty=2); abline(v=43.85, lwd=3, col="pink")
	segments(x0 = binder$Latitude, y0 = binder$lam, x1 = binder$Latitude, y1 = binder$LL, col=BarCols)
	segments(x0 = binder$Latitude, y0 = binder$lam, x1 = binder$Latitude, y1 = binder$UL, col=BarCols)
		points(binder$Latitude, binder$lam,pch=21, cex=1,bg=BarCols)
#		points(binder$Latitude, binder$lamMoist, pch=24, cex=0.8,bg=BarCols)

		#legcol <- c(BarCols[1],BarCols[6], "black")
	#legend("bottomright", legend=c("within","beyond", "Moisture"), pch=c(21, 21, 24), 
#			cex=c(2,2,0.8), col="black", pt.bg=legcol, bg="white", box.col="white")
		z <- lm(binder$Latitude ~ binder$lam)		
		if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
			  print(paste0("lam", " & Latitude - ns")) }
#		z <- lm(binder$Latitude ~ binder$lamMoist)		
#		if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
#			  print(paste0("lam", " & Latitude - ns")) }
	
	par(mai=c(1.6,0.6,0.464,1.2))
	
# PLOT fitness of (elevation vs. fitness)
	plot(enm_el$Elevation, enm_el$fit, ylab ="Relative fitness", pch=21, cex=1,
		ylim=c(min(enm_el$fitLL)*0.95, max(enm_el$fitUL)*1.05),
		bg=c("coral2", "coral2", "blue4", "blue4"), 
		xlab="Elevation (m)", main="Angert and Schemske (2005)", cex.main=1)
		#abline(h=0, col="gray", lty=2)
		abline(v=1450, lwd=3, col="pink")
		segments(x0 = enm_el$Elevation, y0 = enm_el$fit, x1 = enm_el$Elevation, y1 = enm_el$fitLL, col=c("coral2", "coral2", "blue4", "blue4"))
		segments(x0 = enm_el$Elevation, y0 = enm_el$fit, x1 = enm_el$Elevation, y1 = enm_el$fitUL, col=c("coral2", "coral2", "blue4", "blue4"))
		points(enm_el$Elevation, enm_el$fit, ylab ="", pch=21, cex=1, bg=c("coral2", "coral2", "blue4", "blue4"))

	#plot.new()
	#plot.new()

	par(xpd=TRUE)
	legcol <- c(BarCols[1],BarCols[6])
	legend("bottom", legend=c("within","beyond"), inset=-1.05, pch=c(21, 21), 
			pt.cex=c(1,1),col="black", pt.bg=legcol, bg="white", box.col="black")
		
setwd(path.fig)
dev.off()



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
setwd(path.fig)	
pdf(file = "03.5_IPM_finals_vital_rates_ST.pdf", width=(8.5 - (1.25 + 0.75)), height=(11 - (0.75 + 0.75)), family="Times")
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/


# PLOTS NICHE MODEL PREDICTIONS 
		binder$LRavgST <-binder2$LRavg 
		binder$GAMavgST <-binder2$GAMavg 
		binder$RFavgST <-binder2$RFavg 
		binder$BRTavgST <-binder2$BRTavg 
		binder$MAXavgST <-binder2$MAXavg 
		binder$ENavgST <-binder2$EN 
		binder$EN_LST <-binder2$EN_L 
		binder$EN_UST <-binder2$EN_U
		
		mods <- c("LRavgST", "GAMavgST", "RFavgST", "BRTavgST", "MAXavgST")
		modNam <- c("GLM", "GAM", "RF", "BRT", "MAX")

	for(i in 1:length(mods)){
		plot(binder[,c(mods[i])], binder$Surv, ylab ="Relative Survivorship", pch=21, cex=1, 
			bg=BarCols,
			xlab= "Predicted Suitability", main=paste0(modNam[i]))
			z <- lm(binder$Surv ~ binder[,c(mods[i])])		
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0(mods[i], " & Latitude - ns")) }
			mtext(bquote(italic(R)^2 == .(format(summary(z)$adj.r.squared, digits = 2))), side=3, line=-1, adj = 1, cex=0.7)
			
		plot(binder[,c(mods[i])], binder$Grow, ylab ="Relative Growth", pch=21, cex=1,
		bg=BarCols,
			xlab= "Predicted Suitability", main=paste0(modNam[i]))
			z <- lm(binder$Grow ~ binder[,c(mods[i])])	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0(mods[i], " & Latitude - ns")) }
								mtext(bquote(italic(R)^2 == .(format(summary(z)$adj.r.squared, digits = 2))), side=3, line=-1, adj = 1, cex=0.7)
				
		plot(binder[,c(mods[i])], binder$Repr, ylab ="Probability of Flowering", pch=21, cex=1,
		bg=BarCols,
			xlab= "Predicted Suitability", main=paste0(modNam[i]))
			z <- lm(binder$Repr ~ binder[,c(mods[i])])	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0(mods[i], " & Latitude - ns")) }
									mtext(bquote(italic(R)^2 == .(format(summary(z)$adj.r.squared, digits = 2))), side=3, line=-1, adj = 1, cex=0.7)
			
		plot(binder[,c(mods[i])], binder$Fec, ylab ="Relative Fecundity", pch=21, cex=1,
		bg=BarCols, 
			xlab= "Predicted Suitability", main=paste0(modNam[i]))
			z <- lm(binder$Fec ~ binder[,c(mods[i])])	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0(mods[i], " & Latitude - ns")) }
						mtext(bquote(italic(R)^2 == .(format(summary(z)$adj.r.squared, digits = 2))), side=3, line=-1, adj = 1, cex=0.7)


				}	

	dev.off()
	
	
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# ENM VERSU FITNESS FOR LATITUDE AND ELEV
# ENSEMBLE

setwd(path.fig)	
	pdf(file = "03.6_IPM_finals_ENM_lambda_ST.pdf", width=(8.5 - (1.25 + 0.75)), height=(6.0), family="Times")
	par(ps=10, mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	
	# PLOT LAMBDAS 	- ENSEMBLE & NRL
		par(mfrow=c(2,2))
		plot(binder$EN, binder$lam, ylab =expression("Population growth rate,  "*lambda), pch=21, cex=0.8,
		bg=BarCols, ylim=c(0, max(binder$lam*1.1)), xlim=c(0, 1),
		xlab="Predicted Suitability (30 year average)", main="Northern Range Limit", cex.main=1)
		#abline(h=1, col="pink", lty=2)
			# CORRELATIONS SIGNIFIGANT 
#				z <- lm(binder$lamMoist ~ binder$EN)	
#					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=1, col="darkgrey") } else {
#					  print(paste0(mods[i], " & Elevation - ns")) }	
				z<- lm(binder$lam ~ binder$EN)	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=1, col="darkgrey") } else {
					  print(paste0(mods[i], " & Elevation - ns")) }
				segments(x0 = binder$EN, y0 = binder$lam, x1 = binder$EN, y1 = binder$LL, col=BarCols2)
				segments(x0 = binder$EN, y0 = binder$lam, x1 = binder$EN, y1 = binder$UL, col=BarCols2)
				segments(x0 = binder$EN, y0 = binder$lam, x1 = binder$EN_L, y1 = binder$lam, col=BarCols2)
				segments(x0 = binder$EN, y0 = binder$lam, x1 = binder$EN_U, y1 = binder$lam, col=BarCols2)
				points(binder$EN, binder$lam,pch=21, cex=1,bg=BarCols)
#				points(binder$EN, binder$lamMoist,pch=24, cex=0.8,bg=BarCols)
			
	# ENSEMBLE & ELEVATION
		plot(enm_el$EN, enm_el$fit, ylab ="Relative fitness", pch=21, cex=1,
		bg=c("coral2", "coral2", "blue4", "blue4"), xlim=c(0,1),  
		xlab="Predicted Suitability (30 year average)", main="Elevation (Angert and Schemske 2005)", cex.main=1)
			z <- lm(enm_el$fit ~ enm_el$EN)	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0("Ensemble", " & Elevation - ns")) }
			#segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN, y1 = enm_el$LL, col=BarCols)
				#segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN, y1 = enm_el$UL, col=BarCols)
				segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN_L, y1 = enm_el$fit, col=BarColsB2)
				segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN_U, y1 = enm_el$fit, col=BarColsB2)
				segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN, y1 = enm_el$fitLL, col=BarColsB2)
				segments(x0 = enm_el$EN, y0 = enm_el$fit, x1 = enm_el$EN, y1 = enm_el$fitUL, col=BarColsB2)	
				points(enm_el$EN, enm_el$fit, pch=21, cex=1,bg=c("coral2", "coral2", "blue4", "blue4"))
				
	# add short term climate predictions to elevation 
	enm_el$ENavgST <- enm_el2$EN
	enm_el$EN_LST <- enm_el2$EN_L
	enm_el$EN_UST <- enm_el2$EN_U
	
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 	SHORT TERM CLIMATE AVERAGE			
	# LATITUDE
		plot(binder$ENavgST, binder$lam, ylab =expression("Population growth rate,  "*lambda), pch=21, cex=0.8,
		bg=BarCols, ylim=c(0, max(binder$lam*1.1)), xlim=c(0, 1),
		xlab="Predicted Suitability (2014 - 2015)", main="Northern Range Limit", cex.main=1)
		#abline(h=1, col="pink", lty=2)
			# CORRELATIONS SIGNIFIGANT 
#				z <- lm(binder$lamMoist ~ binder$ENavgST)	
#					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=1, col="darkgrey") } else {
#					  print(paste0(mods[i], " & Elevation - ns")) }	
				z<- lm(binder$lam ~ binder$ENavgST)	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=1, col="darkgrey") } else {
					  print(paste0(mods[i], " & Elevation - ns")) }
				segments(x0 = binder$ENavgST, y0 = binder$lam, x1 = binder$ENavgST, y1 = binder$LL, col=BarCols2)
				segments(x0 = binder$ENavgST, y0 = binder$lam, x1 = binder$ENavgST, y1 = binder$UL, col=BarCols2)
				segments(x0 = binder$ENavgST, y0 = binder$lam, x1 = binder$EN_LST, y1 = binder$lam, col=BarCols2)
				segments(x0 = binder$ENavgST, y0 = binder$lam, x1 = binder$EN_UST, y1 = binder$lam, col=BarCols2)
				points(binder$ENavgST, binder$lam,pch=21, cex=1,bg=BarCols)
#				points(binder$ENavgST, binder$lamMoist,pch=24, cex=0.8,bg=BarCols)
			
	# ENSEMBLE & ELEVATION
		plot(enm_el$ENavgST, enm_el$fit, ylab ="Relative fitness", pch=21, cex=1,
		bg=c("coral2", "coral2", "blue4", "blue4"), xlim=c(0,1),
		xlab="Predicted Suitability (2001 - 2003)", main="Elevation (Angert and Schemske 2005)", cex.main=1)
			z <- lm(enm_el$fit ~ enm_el$ENavgST)	
					if(anova(z)[1, "Pr(>F)"] < 0.05){ abline(z, lty=2, lwd=2, col="darkgrey") } else {
					  print(paste0("Ensemble", " & Elevation - ns")) }
				segments(x0 = enm_el$ENavgST, y0 = enm_el$fit, x1 = enm_el$EN_LST, y1 = enm_el$fit, col=BarColsB2)
				segments(x0 = enm_el$ENavgST, y0 = enm_el$fit, x1 = enm_el$EN_UST, y1 = enm_el$fit, col=BarColsB2)
				segments(x0 = enm_el$ENavgST, y0 = enm_el$fit, x1 = enm_el$ENavgST, y1 = enm_el$fitLL, col=BarColsB2)
				segments(x0 = enm_el$ENavgST, y0 = enm_el$fit, x1 = enm_el$ENavgST, y1 = enm_el$fitUL, col=BarColsB2)	
				points(enm_el$ENavgST, enm_el$fit, pch=21, cex=1,bg=c("coral2", "coral2", "blue4", "blue4"))
					
		dev.off()
			
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>				
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


