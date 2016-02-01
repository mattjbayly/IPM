#============================================================================================#
# ANALYSIS OF WILD CARDINALIS TRANSPLANT DATA 
#============================================================================================#
# May 26, 2015
# SCRIPT 3: Site level fecundity/recruitment model
# file name: 03_RecruitmentModelComp.R

library(lme4)
library(nlme)
library(glmmML)
require(plyr)
library(arm)
require(eeptools)
require(stringr)
				
#============================================================================================#
# part A: Set directories, load libraries & import data
# Check basic data structure & reformat layout for use in basic IPM analysis
#============================================================================================#

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


#============================================================================================#
# part B: run models and data
#============================================================================================#
# start pdf for script 
setwd(path.fig)
	pdf(file="03_Recruitment_Explor2.pdf", width=11, height=8.5)



	# site level fecundity & recruitment
		setwd(path.obj)
		rec <- read.csv(file="siteFec.csv")
		# basic recruitment model with log transformed recruits & log transformed fruits
		#Rec_model <- get(load("Rec_model_site.rda")) # old model from previous script
		##################################

# ADD SEED COUNT DATA TO DATA FRAME
		
	# load in seed count per fruits (site level averages) 
		setwd(path.dat)
		seedPerF <- read.csv(file="fec2.seed.per.fruit.2010.2011.2012.csv")
		head(seedPerF)
	
	# Missing seed count data for two sites  
		# Little North Fork Middle Fork - named: "Little North Fork Middle Fork Feather"
		# South Fork Tuolumne - use values from North Fork Tuolumne
		seedPerF$site <- revalue(seedPerF$site, c("Little North Fork Middle Fork Feather"="Little North Fork Middle Fork"))
		copy <- seedPerF[which(seedPerF$site=="North Fork Tuolumne"),]
		copy$site <- "South Fork Tuolumne"
		seedPerF <- rbind(seedPerF, copy)
		
		justGrandMean <- seedPerF[,c("site","newgrandmean")]
		just2010 <- seedPerF[,c("site","mean2010")]; just2010$yrSpecific <- paste0(just2010$site, "2010")
		just2011 <- seedPerF[,c("site","mean2011")]; just2011$yrSpecific <- paste0(just2011$site, "2011")

	# merge seed counts to site level fec/rec - grand mean 
		rec2 <- merge(rec, justGrandMean, by.x="SiteID", by.y="site", all.x=TRUE)
	
	# year specific, but remember that year of fruits is the previous year not
	# the current year (e.g Year 2011 means fruit count from 2010 & recruit count from 2011
		rec2$Specific <- paste0(rec2$SiteID, (rec2$Year-1))
		recSp <- merge(rec2, just2010, by.x="Specific", by.y="yrSpecific", all.x=TRUE)
		recSp <- merge(recSp, just2011, by.x="Specific", by.y="yrSpecific", all.x=TRUE)
		# drop out bad rows 
		drops <- c("site.x","site.y", "Specific")
		recSp <- recSp[,!(names(recSp) %in% drops)]
	# merge year specific columns together
		SeedYrSp <- apply(recSp[, c("mean2010","mean2011")],1,function(x) x[!is.na(x)])
		rec <- cbind(recSp, SeedYrSp)
		drops <- c("mean2010","mean2011")
		rec <- rec[,!(names(rec) %in% drops)]
		names(rec)[names(rec)=="newgrandmean"] <- "SeedGrandMean"
		head(rec, 3)
		colnames(rec)
	# add two more variables called normalized seed counts
		
######################################################
# TRY SOME BASIC MODELS 
	# forced to go through the origin
	# can't have an intercept in this model
	# see summary at bottom
	Rec_model2 <- lm(Lrecruits ~ Lfruits*SeedGrandMean*Region -1,  data=rec)
	Rec_model3 <- lm(Lrecruits ~ Lfruits*SeedYrSp*Region -1,  data=rec)
	#summary(Rec_model2)
	#summary(Rec_model3)
	
	
	# basic plots
	palette(c("green", "blue", "red"))
	par(mfrow=c(2,3))
	plot(rec$fruits,rec$recruits,xlab="Fruits per site (t-1)",ylab="Recruits in site at t", main="raw", col=as.factor(rec$Region))
	# basic model - no seeds considered
		plot(rec$Lfruits,rec$Lrecruits,xlab="log Fruits per site (t-1)",ylab="log Recruits in site at t", main="ignore seed count", col=as.factor(rec$Region))
		lm1 <- lm(Lrecruits ~ Lfruits - 1, data=rec)
		abline(lm1); text(3,3,paste0("adj Rsqr " , round(summary(lm1)$adj.r.squared, 2)), cex=1.2)
	# mod * seed count (grand mean)
		rec$FrSeGrand <- (log(rec$fruits+1))*(log(rec$SeedGrandMean+1))
		plot(rec$FrSeGrand,rec$Lrecruits,xlab="log Fruits*seeds per site (t-1)",ylab="log Recruits in site at t", main="Fruit*seed grand mean", col=as.factor(rec$Region))
		lm2 <- lm(Lrecruits ~ FrSeGrand - 1, data=rec)
		abline(lm2); text(3,3,paste0("adj Rsqr " , round(summary(lm2)$adj.r.squared, 2)), cex=1.2)
	# mod * seed count (year specific) 	
		rec$FrSeYr <- (log(rec$fruits+1))*(log(rec$SeedYrSp+1))
		plot(rec$FrSeYr,rec$Lrecruits,xlab="log Fruits*seeds per site (t-1)",ylab="log Recruits in site at t", main="Fruit*seed yearly mean", col=as.factor(rec$Region))
		lm3 <- lm(Lrecruits ~ FrSeYr - 1, data=rec)
		abline(lm3); text(3,3,paste0("adj Rsqr " , round(summary(lm3)$adj.r.squared, 2)), cex=1.2)
	# collapse years - take average for both years 
		library(plyr)
		newRec <- ddply(rec, .(SiteID), summarize, FrSeGrand=mean(FrSeGrand), Lrecruits=mean(Lrecruits), SeedGrandMean=mean(SeedGrandMean), fruit=mean(fruits), recruits=mean(recruits))
		newRec$FrSeGrand <- (log(newRec$fruit+1))*(log(newRec$SeedGrandMean+1))
		dim(newRec) # only one observation per site
		plot(newRec$FrSeGrand,newRec$Lrecruits,xlab="log Fruits*seeds per site (t-1)",ylab="log Recruits in site at t", main="Collapse years", col=as.factor(rec$Region))
		lm4 <- lm(Lrecruits ~ FrSeGrand - 1, data=newRec)
		abline(lm4); text(3,3,paste0("adj Rsqr " , round(summary(lm4)$adj.r.squared, 2)), cex=1.2)
	# model lm4 is the winner here. 
	summary(lm4); AIC(lm2);AIC(lm3); AIC(lm4)

	# are there any North/South regional differences (drop North sites from dataset, then drop south sites
	# from dataset 
		recN <- rec[which(rec$Region=="N" | rec$Region=="C"),]
		newRecN <- ddply(recN, .(SiteID), summarize, FrSeGrand=mean(FrSeGrand), Lrecruits=mean(Lrecruits), SeedGrandMean=mean(SeedGrandMean), fruit=mean(fruits), recruits=mean(recruits))
		newRecN$FrSeGrand <- (log(newRecN$fruit+1))*(log(newRecN$SeedGrandMean+1))
		lm_NC <- lm(Lrecruits ~ FrSeGrand - 1, data=newRecN)
		# then run for south and central only 
		recS <- rec[which(rec$Region=="S" | rec$Region=="C"),]
		newRecS <- ddply(recS, .(SiteID), summarize, FrSeGrand=mean(FrSeGrand), Lrecruits=mean(Lrecruits), SeedGrandMean=mean(SeedGrandMean), fruit=mean(fruits), recruits=mean(recruits))
		newRecS$FrSeGrand <- (log(newRecS$fruit+1))*(log(newRecS$SeedGrandMean+1))
		lm_SC <- lm(Lrecruits ~ FrSeGrand - 1, data=newRecS)
		# compare models
		coffy <- c(coef(lm4), coef(lm_NC), coef(lm_SC))
		connyL <- c(confint(lm4)[1], confint(lm_NC)[1], confint(lm_SC)[1])
		connyU <- c(confint(lm4)[2], confint(lm_NC)[2], confint(lm_SC)[2])
		adjr2 <- c(round(summary(lm4)$adj.r.squared, 2),round(summary(lm_NC)$adj.r.squared, 2), round(summary(lm_SC)$adj.r.squared, 2))
		moddy <- c("all-lm4", "lm_NC", "lm_SC")
		(regional <- data.frame(modName=moddy, est=round(coffy, 4), LL=round(connyL, 4), UL=round(connyU, 4), AdjRsqur=adjr2, region=c("C", "N", "S")))
		# each of the north/central/south regions should use different estimates. 
		par(mfrow=c(1,1))
		plot(newRec$FrSeGrand,newRec$Lrecruits,xlab="log Fruits *seeds per fruit (t-1)",ylab="log Recruits in site at t", main="Regional differences", col=as.factor(rec$Region), ylim=c(0,5), xlim=c(0,60))
		abline(lm4, col="green")
		abline(lm_NC, col="blue")
		abline(lm_SC, col="red")

		
		
		
		
		
		
######################################################
# Rather than using the log transformed version of size what if we 
# just set the seed output to relative seed output per site in relation to the mean
# e.g. Seed output at Carlon = (Absolute seed output at Carlon)/(Average Seed output across all sites)
# but this might cause problems because values will be around 1.  
	
	# SUMMARY:
	# The data (as expected) is highly variable no model fits extremely well. But it seems 
	# like collapsing years and just looking at the relationship between the number of 
	# of fruits dropped at a site (over all years) and the number of recruits (over all years)
	# seems to give the best overall approximation.  The intercept is removed because 
	# can’t have a situation in the IPM where we simply add a site specific arbitrary number 
	# of new recruits each year. The estimate across the different models ranged from 
	# 0.06 – 0.076, but it is important to remember that these are on a weird log scale & need to be converted before entering the IPM. 
	# log(recruits #) = [log(fruit count+1)*log(seed/fruit +1)]*coef
	

	
#============================================================================================#
# part C: converting back & forth between transformed & untransformed data in the IPM 
#============================================================================================#
	# Methods from:
		# M. Rees, D.Z. Childs, and S.P. Ellner. Building Integral Projection Models: a User's Guide. Journal of Animal Ecology, 2014
		# Appendix S1: Probability densities vs. probabilities, and moving a kernel from one measurement scale to another
	
	# e.g. Probability of a seed surviving (assume here 1500 seeds/fruit)
	# coef = ~0.07, so for every FrSeGrand to get 100 recruits. 0.076*X = 100,  x=1313.232
	temp <- newRec
	temp$predicted_rec_NCS <- exp(predict(lm4, newdata=temp))
	temp$predicted_rec_NC <- exp(predict(lm_NC, newdata=temp))
	temp$predicted_rec_NS <- exp(predict(lm_SC, newdata=temp))

	# Format to work with 
	head(temp, 3)
	#EXAMPLE	recruit number = exp((log(SEEDSPERFRUIT+1)*log(FRUITPRODUCED+1))*(coef(MODEL)[[1]]))
	exp((log(1955.67+1)*log(135.7667+1))*(coef(lm4)[[1]]))	# expected number of recruits
	exp((log(1184.00+1)*log(2.0000+1))*(coef(lm4)[[1]]))	# expected number of recruits
	exp((log(1226.47+1)*log(489.0000+1))*(coef(lm4)[[1]]))	# expected number of recruits

#####################
# Fruit to recruit conversion ratio - for every fruit dropped there will be this many recruits next year: 
	mean((temp$predicted_rec_NCS)/((temp$fruit))) # in Center
	mean((temp$predicted_rec_NC)/((temp$fruit))) # in North
#####################
# for every seed there will be this many recruits produced the next year
	 mean((temp$predicted_rec_NCS)/((temp$SeedGrandMean)*(temp$fruit))) # if in the center, or 
	 mean((temp$predicted_rec_NC)/((temp$SeedGrandMean)*(temp$fruit))) # if in the north 
 

	####################
	# save model coefficients and seed counts for use in IPMs 
	setwd(path.obj)
 	(seedreg <- ddply(rec,~Region,summarise,meanSeedcount=mean(SeedGrandMean),Seedsd=sd(SeedGrandMean)))
	(mytable <- merge(seedreg, regional, by.x="Region", by.y="region"))
	mytable$Formula = "log(recruits + 1) = [log(seeds per fruit + 1)*log(number of fruits + 1)]*est"
	mytable
	write.csv(mytable, file="RecruitmentModCoefs.csv", row.names=FALSE)
 
#============================================================================================#
# part D: Other outstanding problems  
#============================================================================================#
# A. Intercept is still far >0 when we don't manually set it to zero
	# this means 
# B. What about all the young of the year recruits that 
	# incorrectly enter the model & recruit kernels in the incorrect year.
	# eg. seeds drop from plants in 2012 in July/August and seedlings germinate, but 
	# as far as the IPM is concerned these individuals are assumed to be from the 
	# fecundity of the previous year 2011. 
	# Optimal model would have a YOY (young of the year fecundity component). 
# C. Check model assumptions & violations ect. 

###########################################
# How well do these recruitment estimates above (estimated at the site level)
# predict plot level recruitment or region level recruitment. 
	setwd(path.obj)
	dir()
	
	plotTest <- read.csv(file="plotFecRec.csv")
	regionTest <- read.csv(file="regionFecRec.csv")

	# start with REGION first 
	(seedreg <- ddply(rec,~Region,summarise,meanSeed=mean(SeedGrandMean),sd=sd(SeedGrandMean)))
	require(stringr)
	regionTest$regNCS <- str_extract(regionTest$Region,"[[:upper:]]")
	regionTest <- merge(regionTest, seedreg, by.x="regNCS", by.y="Region", all.x=TRUE, all.y=FALSE)
	regionTest$FrSeGrand <- (log(regionTest$fruits+1))*(log(regionTest$meanSeed+1))
	regionTest$Lrecruits <- log(regionTest$recruits +1)
	plot(regionTest$FrSeGrand, regionTest$Lrecruits, xlab="log seed*fruit", ylab="log recruits", 
		col=as.factor(regionTest$regNCS), pch=16, main="Regional check of site level model")
		abline(lm4, col="green")
		abline(lm_NC, col="blue")
		abline(lm_SC, col="red")
		# doesn't look great but probably best possible from dataset 
		
	# then check PLOT level, just set average seed/fruit to  1622 (grand mean) for now	
	plotTest$FrSeGrand <- (log(plotTest$fruits + 1))*(log(1622+1))
	plotTest$Lrecruits <- log(plotTest$recruits + 1)
	plot(plotTest$FrSeGrand, plotTest$Lrecruits, xlab="log seed*fruit", ylab="log recruits", 
		pch=16, main="Plot level check of site level model")
		abline(lm4, col="green")
		abline(lm_NC, col="blue")
		abline(lm_SC, col="red")
		# very noisy, but ok
		
# close plotting device 
setwd(path.fig)
dev.off()
