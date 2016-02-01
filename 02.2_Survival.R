###########################################################
# SURVIVAL FUNCTIONs FOR 2013 - 2015 CARDINALIS TRANSPLANT DATA
	# This script is run remotely from "02.1_basicIPM_Transplant.R"

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table on content:
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# A. Exploratory plots across levels & factors from basic glm
			# Are there any major outliers ect. Was the log transformation of size sufficient? 
			
	# B. Determine appropriate random effect structure (single random effect of sites vs. plots within sites ect.) 
			# for glmer with PlotID & SiteID as random effects. 
			
	# C. Determine appropriate fixed structure (region, year, size, sites, interactions?)
			# for survival glmer with a Akaike-weights framework. 
			
	# D. Make post-hoc diagnostic plots 
			# are basic assumptions violated/any other problems? 
			
	# E. Bootstrap to get approximations for estimates & standard errors to use in final model. 
	
	# F. Save final results for incorporation into IPM. 

	# G. Make Response Plots.

#
##
###
####
#####
######
#######
######
#####
####
###
##
#

# SETUP 
	# load libraries
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
	dfclean <- dfclean[order(dfclean$z),]	# sort smallest -> largest (z) 

###############################################################################
	# Need to here (before running the rest of scripts & subsetting data ect)
	# make bootstrap replicates from the data. These bootstrap replicates will be used again 
	# in further scripts. 
	# Due to random effect structure will have to rely on bootstrapped estimates for parameters of 
	# interest to account for variation across levels. 
	#############################################################################	
	# Sorry gang... please install & load all these libraries too :(
	library(lme4); require(ggplot2); require(GGally); require(reshape2); require(lme4); require(compiler)
	require(parallel); require(boot)		
	# bootstrapping function
						sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
						cid <- unique(dat[, clustervar[1]])
						ncid <- length(cid)
						recid <- sample(cid, size = ncid * reps, replace = TRUE)
						if (replace) {
						rid <- lapply(seq_along(recid), function(i) {
						cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
						size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
						})
						} else {
						rid <- lapply(seq_along(recid), function(i) {
						cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
						})
						}
						dat <- as.data.frame(do.call(rbind, rid))
						dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
						labels = FALSE))
						dat$NewID <- factor(dat$NewID)
						return(dat)
						}
	# Resample data to a new dfs  
		set.seed(20)
		tmp <- sampler(dfclean, "SiteID", reps = 5000) # run at 5000 for now boot strap replaicate datasets, this means will get about ~ 800,000 observations
		bigdata <- cbind(tmp, dfclean[tmp$RowID, ])
		dim(dfclean)
		dim(bigdata)
		setwd(path.obj)
		write.csv(bigdata, file="BootReplicateDataSets.csv", row.names=FALSE)
		####################################################################################
		# Above is the 5000 random samples from the original dataframe. Next we will sample 
		# at the site level and then at the plot level to maintain the data structure. Final site level
		# confidence intervals will be compared using the following methods. 
			# 1. 5000 bootstrap replicates from the entire dataframe, sampled at random with replacement. 
			# 2. 5000 bootstrap replicates AT THE SITE LEVEL for each dataframe (# plants in each site is held constant for the replicates). 
			# 3. 5000 bootstrap replicates AT THE PLOT LEVEL for each dataframe (# plants in each plot is held constant for the replicates). 
#...... Sites............... 
		Tempdf <- dfclean
		Tempdf$RowID <- rownames(Tempdf)
		SampSites <- Tempdf[,c("SiteID", "RowID")]
		#SampSites$RowID <- rownames(dfclean)
		Usites <- unique(SampSites$SiteID)
		
		Mframe <- data.frame()	
			for(i in 1:length(Usites)){
				samp <- SampSites[which(SampSites$SiteID==Usites[i]),]
				mySamp <- sample(samp$RowID, ((dim(samp)[1])*5000), replace = TRUE)
				Replicate <- rep(1:5000, each=dim(samp)[1])
				NewID <- Replicate
				toBind <- data.frame(NewID=NewID, Replicate=Replicate, RowID= mySamp)
				Mframe <- rbind(Mframe, toBind)
			}
			siteBoot <- merge(Mframe, Tempdf, by.x="RowID", by.y="RowID", all.x=TRUE, all.y=FALSE)
			siteBoot <- siteBoot[,c("NewID", "RowID", "Replicate", "Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_13", "NewPlot_14", "Year")]
					
		setwd(path.obj)
		write.csv(siteBoot, file="BootReplicateDataSets_Sites.csv", row.names=FALSE)

#...... Plots............... 
		Tempdf <- dfclean
		Tempdf$RowID <- rownames(Tempdf)
		SampPlots <- Tempdf[,c("PlotID", "RowID")]
		#SampPlots$RowID <- rownames(dfclean)
		UPlots <- unique(SampPlots$PlotID)
		
		Mframe <- data.frame()	
			for(i in 1:length(UPlots)){
				samp <- SampPlots[which(SampPlots$PlotID==UPlots[i]),]
				mySamp <- sample(samp$RowID, ((dim(samp)[1])*5000), replace = TRUE)
				Replicate <- rep(1:5000, each=dim(samp)[1])
				NewID <- Replicate
				toBind <- data.frame(NewID=NewID, Replicate=Replicate, RowID= mySamp)
				Mframe <- rbind(Mframe, toBind)
			}
			PlotBoot <- merge(Mframe, Tempdf, by.x="RowID", by.y="RowID", all.x=TRUE, all.y=FALSE)
			PlotBoot <- PlotBoot[,c("NewID", "RowID", "Replicate", "Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_13", "NewPlot_14", "Year")]
					
		setwd(path.obj)
		write.csv(PlotBoot, file="BootReplicateDataSets_Plots.csv", row.names=FALSE)
			
#---------------------------------------------------------------------------------
	# need to exclude the 2015-2016 flower/fecundity observations from dataset
	# in this script
	dfclean <- dfclean[!is.na(dfclean$Surv),]	# exclude NA's from Flr & Fec 
#---------------------------------------------------------------------------------
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part A: Exploratory plots across levels & factors from basic GLM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Are there any major outliers ect. Was the log transformation of size sufficient (or too much)? 
	# Will run exploratory plots for each site, year & region combination and examine visually. 
	# Combine plots into a large pdf. A bit tedious, but informative & useful before proceeding 
	# with the remainder of the analysis. 
	#############################################################################	
			dev.off(); setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
			setwd(path.fig)#
			pdf(file="02.2.1_Survival_Explor.pdf", width=11, height=8.5)
			z <- dfclean$z # stand along object
			
				# Plotting frame 
					surv.plot.data <- within(dfclean, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = dfclean)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="Global model")
			points(z, jitter(dfclean$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = dfclean, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			
		# Run additional plots (optional code stored in R Code - functions subfolder)
		# additional plots in another script to save space on this sheet & make for efficent coding 
			setwd(path.funct)
			source("02.2_extra_plots_Surv.R")
			# check for figure in figures folder - looks OK 
		setwd(path.fig)
		dev.off() # close pdf & open to view



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part B: Determine appropriate random effect structure
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Many random effect structures are possible for a mixed effect model 
	# Especially when there are both plot and site level nesting (region/site/plot/individual). 
	# It might even be possible to have a random slope component in addition 
	# to a random intercept. 
	# Unfortunately defendable statistical methods to choose an optimal 
	# random effect structure for glm objects are still under development. 
	# discussion: http://www.researchgate.net/post/How_can_I_optimize_the_random_effect_structure_in_a_GLMM
	# good discussion of issue here too: http://glmm.wikidot.com/faq
	# Standard REML procedure not possible here for glm objects/ logistic regression. Should rely on own 
	# intuition, experimental setup and a solid hypothesis for each alternative 
	# random effect structure. 
	
		# Possible Alternatives (listed from most believable to least believable): 
				# (1|PlotID) random plot level intercept
				# (z|PlotID) random plot level intercept and random plot level slope
				# (1|SiteID/PlotID) random site level intercept and random site by plot level intercept
				# (1|SiteID/PlotID) random site level & slope intercept and random site by plot level intercept & slope
				# site + (z|PlotID) - fixed effect of site, plus random plot level int & slope
			# Not considered: 
				#F - (0+z|SiteID) = (-1+z|SiteID)	random slope of z within SiteID: no variation in intercept
				#G - (1|SiteID) + (0+z|SiteID)	uncorrelated random intercept and random slope within SiteID
	
	#############################################################################	
	# To avoid over-fitting will use cross-validation with successively smaller sub-samples of
	# the data (training & testing sets) to compare the complexity of alternative 
	# random effect structures. If training and testing sets are cut randomly from the 
	# full data frame then it might be possible to detect over fitting. 
		library(lme4)
		library(arm)
		library(sjPlot)
		
		
		holdOutPortions <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.95) # Portion of data to be used in each trial
		myFrame <- data.frame()
		for(i in 1:length(holdOutPortions)){
			# Function to split data into training & testing set (e.g. 50% training 50% testing, when holdOutPortions = 0.5)
				splitdf <- function(dataframe, seed=NULL) {
				if (!is.null(seed)) set.seed(seed)
				index <- 1:nrow(dataframe)
				trainindex <- sample(index, trunc(length(index)*holdOutPortions[i]))
				trainset <- dataframe[trainindex, ]
				testset <- dataframe[-trainindex, ]
				list(trainset=trainset,testset=testset)
			}
		
		splits <- splitdf(dfclean, seed=808) # splits dataframe 
		# str(splits) # training & testing set
		# Observations in each data frame
		lapply(splits,nrow)
		#view the first few columns in each data frame
		lapply(splits,head)
		# save the training and testing sets as data frames
		training <- splits$trainset
		testing <- splits$testset
		
		# To avoid having unique PlotID levels to training 
		# and testing sets add just one observation from each plot to each dataset
			df2 <- lapply(split(dfclean, dfclean$PlotID), function(subdf) subdf[sample(1:nrow(subdf), 1, replace=TRUE),])		
			df2 <- do.call('rbind', df2)
			head(df2, 3); dim(df2)
			training <- rbind(training, df2)
			testing <- rbind(testing, df2)

	# Run glmer with all potential fixed effect terms (no three way interaction)
		testA_plot <- glmer(Surv ~ z + (1|PlotID), na.action=na.omit, data=training, family = binomial)
		testA_site <- glmer(Surv ~ z + (1|SiteID), na.action=na.omit, data=training, family = binomial)
		testA_both <- glmer(Surv ~ z + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=training, family = binomial)
		testA_neither <- glm(Surv ~ z, na.action=na.omit, data=training, family = binomial)

		
	# Predict models to testing data
		ModPlotPreds <- predict(testA_plot, newdata = testing, type="response")
		ModSitePreds <- predict(testA_site, newdata = testing, type="response")
		ModBothPreds <- predict(testA_both, newdata = testing, type="response")
		ModNeithPreds <- predict(testA_neither, newdata = testing, type="response")
		
	library(pROC)
	testing$ModPlotPreds <- ModPlotPreds
	testing$ModSitePreds <- ModSitePreds
	testing$ModBothPreds <- ModBothPreds
	testing$ModNeithPreds <- ModNeithPreds
		
	modPlot <- roc(Surv ~ ModPlotPreds, data = testing); PlotA <- modPlot$auc[1]
	modSite <- roc(Surv ~ ModSitePreds, data = testing); SiteA <- modSite$auc[1]
	modBoth <- roc(Surv ~ ModBothPreds, data = testing); BothA <- modBoth$auc[1]
	modNeith <- roc(Surv ~ ModNeithPreds, data = testing); NeithA <- modNeith$auc[1]
	holdout <- holdOutPortions[i]
	N <- dim(training)[1]*holdOutPortions[i]
	
	toAdd <- data.frame(holdout, N, PlotA, SiteA, BothA, NeithA) # add results from each run to 
	
	myFrame	<- rbind(myFrame, toAdd)
	} # 
	
	# save final evaluation figure
	setwd(path.fig)
	pdf(file="02.2.2_Sur Random Effect Structure.pdf", width=11, height=8.5)
	colnames(myFrame) <- c("trainPor", "N", "PlotA", "SiteA", "BothA", "NeithA") 
	plot(myFrame$trainPor, myFrame$PlotA, xlab="Portion of data Training:Testing ratio", 
		ylab="AUC", col=1, ylim=c(min(myFrame[,3:6]) - 0.04, max(myFrame[,3:6]) + 0.04), pch=19, cex=1.5, main="Survivorship Random Effect Structure")
	lines(myFrame$trainPor, myFrame$PlotA, col=1)
	points(myFrame$trainPor, myFrame$SiteA, col=2, pch=19); lines(myFrame$trainPor, myFrame$SiteA, col=2)
	points(myFrame$trainPor, myFrame$BothA, col=3, pch=19); lines(myFrame$trainPor, myFrame$BothA, col=3)
	points(myFrame$trainPor, myFrame$NeithA, col="darkgrey", pch=19); lines(myFrame$trainPor, myFrame$NeithA, col="darkgrey")
	legend('topright', legend = c("(1|PlotID)", "(1|SiteID)", "(1|PlotID) + (1|SiteID) + (1|SiteID:PlotID)", "glm"), col = c(1:3, "darkgrey"), cex = 0.8, pch = 19)
	text(myFrame$trainPor, min(myFrame[,3:5]) - 0.02, labels = round(myFrame$N, 0))
	
	# plot of coefficients
		# https://strengejacke.wordpress.com/sjplot-r-package/
		
		# plot fixed effects
		#par(mfrow=c(1,3))
		#sjp.glmer(testA_plot, type = "fe")
		#sjp.glmer(testA_site, type = "fe")
		#sjp.glmer(testA_both, type = "fe")
	
	dev.off()
	
	# Having plot alone seems to work better than having site alone, but results seems counter 
	# intuitive to experimental design & field observations. 
	# Will have to reassess and see how sensitivity final results are to estimates from each random 
	# effect structure. 
	# Unlike the demography dataset we will use plot as the main random effect in our models.
	# Having only a site level intercept did not result in dramatic improvements from the basic glm (non-mixed model). 
	# There was little to no within plot variation because of the experimental design that standardized
	# conditions within plots. 
	# One concern is the low level of replication among some plots. 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part C: Determine appropriate fixed structure (region, year, size, sites, interactions?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# For survival glmer with a Akaike-weights framework. 
	# Largest possible full model:
		#"Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region"
	# Don't want year to go into model, so should potentially average predictions for 
	# model with & without year or evaluate difference whether year is included or excluded. 
	# Procedure followed here from: 
		# Chapter 21: GLMM applied on the spatial distribution of Koalas in a fragmented landscape
			# section 24.4.2 (pg. 483) Zuur 2009 - "Mixed effect models and extension in ecology with R"
			# code from: https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	#############################################################################	
	
	# Possible models(yes could be coded more efficiently.. sorry)  
	#mod2 <- glmer(Surv ~ z + SiteID + Region + z*SiteID + z*Region + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	#mod3 <- glmer(Surv ~ z + SiteID + Region + z*SiteID + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	#mod4 <- glmer(Surv ~ z + SiteID + Region + z*Region + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	#mod5 <- glmer(Surv ~ z + SiteID + Region + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod6 <- glmer(Surv ~ z + Region + z*Region + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod7 <- glmer(Surv ~ z + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod8 <- glmer(Surv ~ z + Region + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod9 <- glmer(Surv ~ z + SiteID + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod10 <- glmer(Surv ~ z + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)

	# package for getting effects & AIC values 
	require(MuMIn)
	
	#################
	# MODEL SELECTION 
	#################
		# use the mod.sel function to conduct model selection
		# and put output into object out.put
		out_put <- mod.sel(mod6,mod7,mod8,mod9,mod10)
		modList <- c("mod6","mod7","mod8","mod9","mod10")
		out_put
		
	# What about models with a random slope and intercept component: 
		mod6B <- glmer(Surv ~ z + Region + z*Region + (z|PlotID), na.action=na.omit, data=dfclean, family = binomial)
		mod7B <- glmer(Surv ~ z + SiteID + z*SiteID + (z|PlotID), na.action=na.omit, data=dfclean, family = binomial)
		mod8B <- glmer(Surv ~ z + Region + (z|PlotID), na.action=na.omit, data=dfclean, family = binomial)
		mod9B <- glmer(Surv ~ z + SiteID + (z|PlotID), na.action=na.omit, data=dfclean, family = binomial)
		mod10B <- glmer(Surv ~ z + (z|PlotID), na.action=na.omit, data=dfclean, family = binomial)
		out_putB <- mod.sel(mod6B,mod7B,mod8B,mod9B,mod10B)
		modListB <- c("mod6B","mod7B","mod8B","mod9B","mod10B")
		out_putB # looks the same out_put
	
	# Compare all at once - just to see 
	out_putC <- mod.sel(mod6,mod7,mod8,mod9,mod10,mod6B,mod7B,mod8B,mod9B,mod10B)
	out_putC

	# TOP MODELS
		# OPTIONAL - select models 95% cumulative weight criteria 
		# IMPORTANT: Weights have been renormalized!! 
		subset(out_put, cumsum(out_put$weight) <= .95) 
		
	# MAKE CLEAN AICc TABLE & SAVE
		# coerce the object out_put into a data frame 
		# elements 6-10 in out_put have what we want 
		sel.table <- as.data.frame(out_put)[7:11] 			
		# a little clean-up, lets round things a bit 
		sel.table[,2:3]<- round(sel.table[,2:3],2) 
		sel.table[,4:5]<- round(sel.table[,4:5],3) 
		# that’s better 
		sel.table 
		## lets be sure to put the model names in a column 
		sel.table$Model<-rownames(sel.table) 
		# replace Model name with formulas little tricky so be careful 
		for(i in 1:nrow(sel.table)) sel.table$Model[i]<- as.character(formula(paste(sel.table$Model[i])))[3] 
		# let's see what is in there 
		sel.table 
		#little reordering of columns 
		sel.table<-sel.table[,c(6,1,2,3,4,5)] 
		sel.table
		# write to a file, here a comma separated values format 
		# make sure your working directory is properly specified 
		setwd(path.obj)
		write.csv(sel.table,"SurvModAICTable.csv", row.names = T) 
			
	# VARIABLE IMPORTANCE (BASED ON AIC)	
		importance(out_put)
		importance(out_putC)

	# MODEL AVERAGING (over top models) ~ can't really do this properly with glm's wont work 
			MA.ests <- model.avg(out_put, revised.var = TRUE) 
			MA.ests$avg.model 
			#Here are the beta tilda bar MA estimates 
			MA.ests$coef.shrinkage
			# you can also obtain importance weights for individual params 
			MA.ests$importance 
			#create model averaged estimates for parameters in confidence set of 
			# models only using subset command 
			MA.ests<-model.avg(out_put, subset= delta < 5, revised.var = TRUE) 
			MA.ests$avg.model	
			# lets clean up a bit and write the table to a file 
			MA.est.table<-round(MA.ests$avg.model[,c(1:2,4:5)],6) 
			MA.est.table
			## write to a CSV file 
			setwd(path.obj)
			#write.csv(MA.est.table, "SurvModAvgEst.csv") 
			# Unfortunately estimating response values from these average 
			# parameters won't work for logistic regression.  
			
	# MAKE PREDICTIONS	
			# extract parameters and weights from confidence model set 
			# using get.models function 
			pred.parms<-get.models(out_put, subset= delta < 5) 
			# predict values using each model, here were just using the 
			# the example dataset, you could use a new dataset 
			model.preds = sapply(pred.parms, predict(type="response"), newdata = dfclean) 			
			# One word of warning, do not estimate model averaged parameters
			# for mixed models! You can, however, model average the predictions of GLMM. 	
			# All of the functions above can be used with objects created by
			# fitting linear models with the glm function. All of the above
			# applies to these glm objects. However, there is one important
			# difference. GLM’s such as Poisson regression
			# and Logistic regression can often fail to meet statistical
			# assumptions due to extra variance. This is often defined 
			# as over-dispersion (not an issue for normal linear regression!)
			# and requires the use of quasi-AIC for model selection.
		
		# Check for overdispersion 
			library(RVAideMemoire)
			overdisp.glmer(mod9) # seems ok, close to one

		# winning model:
			summary(mod9)
		
		library(glmmML) # run sample model with glmmML - do coefficents change?
		mod9_glmmML <-glmmML(formula=Surv ~ z + SiteID,
             cluster=PlotID,data=dfclean,family=binomial)
		summary(mod9_glmmML)
	

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part D: Make post-hoc diagnostic plots 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Are basic assumptions violated/any other problems? 
	# Procedure followed here from: 
		# Chapter 21: GLMM applied on the spatial distribution of Koalas in a fragmented landscape
			# Zuur 2009 - "Mixed effect models and extension in ecology with R"
		
		# 1. From the fitted model, calculate the residuals (ri). 
		# 2. Order the ri, giving r(i).
		# 3. Simulate M data sets from the fitted model.
		# 4. Fit the model to the M data sets.
		# 5. Compute the residuals ri for the models fitted to the M data 
		# 		sets and order them to get r(i)
		# 6. Calculate the medians of the ordered residuals from the M replicates. (Landwehr
		#		et al. (1984) use a slight modification here where they interpolate within the
		#		distribution of the simulated residuals to avoid plotting negative against positive residuals.)
		# 7. Plot the median simulated ordered (interpolated) residuals against the ordered
		#		residuals from the original model fit.
		# 8. Calculate confidence intervals for the simulated ordered (interpolated) residuals
		# 		from the M replicates.
		#9. Plot the median simulated ordered (interpolated) residuals against the upper and
		#		lower confidence intervals.
	
	#############################################################################	
	# load in model adequecy fuctions from Zuur. 
	setwd(path.funct)
	source("ZuurGLMMAdequecyFxns.R") # downloaded from: http://www.highstat.com/book2.htm
	
	# QUANTILE-QUANTILE PLOTS (note that the functions "fitted.glmmML" and "res.glmmML" are requiured here, loaded in line above)	
		#calculate the fitted values for the model
		Model_Best <- mod9_glmmML # best model from steps above
		Fitted <- fitted.glmmML(model=Model_Best, data=dfclean)
		#calculate the ordered residuals for the model
		Resids<-sort(res.glmmML(model=Model_Best,data=dfclean))
		Resids_0<-Resids[Resids<0]
		Resids_1<-Resids[Resids>=0]
		#specify the number of replicates
		Reps<-500
		#simulate Reps data sets
		Sims<-matrix(rbinom(n=length(Fitted)*Reps,size=1,p=Fitted),nrow=length(Fitted),ncol=Reps)
		#fit the model to each simulated data set
		Models<-apply(Sims,MARGIN=2,FUN=function(X){return(list(X,glmmML(formula=Surv ~ z + Region + z * Region,cluster=SiteID,data=dfclean,family=binomial)))})
		#calculate the ordered (interpolated) simulated residuals and the point-wise 95% confidence intervals
		Resids_Sim<-matrix(unlist(lapply(Models,FUN=function(X){TempData<-dfclean;TempData[,"Surv"]<-X[[1]];return(sort(res.glmmML(X[[2]],TempData)))})),ncol=Reps,nrow=nrow(dfclean))
		Resids_Sim_0<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X<0],ppoints(Resids_0,a=1))}),ncol=Reps,nrow=length(Resids_0))
		Resids_Sim_1<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X>=0],ppoints(Resids_1,a=1))}),ncol=Reps,nrow=length(Resids_1))
		Resids_Sim_0_Median<-apply(Resids_Sim_0,MARGIN=1,FUN=median)
		Resids_Sim_1_Median<-apply(Resids_Sim_1,MARGIN=1,FUN=median)
		Resids_Sim_0_Lower<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.025)})
		Resids_Sim_0_Upper<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.975)})
		Resids_Sim_1_Lower<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.025)})
		Resids_Sim_1_Upper<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.975)})
		#plot the qauntile-quantile plot with 95% confidence intervals and 1:1 line
		# save plots
		dev.off()
		setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
			pdf(file="02.2.3_Survival_model_adequacy.pdf", width=11, height=8.5)
		plot(Resids_Sim_0_Median,Resids_0,xlim=c(-1,1),ylim=c(-1,1),
			xlab="simulated quantiles",ylab="fitted quantiles", cex=0.5, col="#0000FF0A",
			main="Quantile-quantile plot with 95% pointwise confidence bounds")
		points(Resids_Sim_1_Median,Resids_1, cex=0.5, pch=16, col="#0000FF0A")
		lines(Resids_Sim_0_Median,Resids_Sim_0_Lower)
		lines(Resids_Sim_0_Median,Resids_Sim_0_Upper)
		lines(Resids_Sim_1_Median,Resids_Sim_1_Lower)
		lines(Resids_Sim_1_Median,Resids_Sim_1_Upper)
		abline(0,1,lty=3)

	# PARTIAL RESIDUAL PLOTS (note that the functions "fitted.glmmML" and "res.glmmML" are required here)
		#calculate the partial residuals for each covariate
		# see formula on pg. 489 of Zuur book cited above
		Part_Res1<-((res.glmmML(Model_Best,dfclean)-fitted.glmmML(Model_Best,dfclean))/
		((fitted.glmmML(Model_Best,dfclean)*(1-fitted.glmmML(Model_Best,dfclean)))))+
		(Model_Best$coefficients[2]*dfclean[,"z"])	
		#plot the partial residuals and the smoothed plot
		plot(dfclean[,"z"],Part_Res1,ylab="partial residuals", cex=0.5, pch=16, 
			col= "#0000FF4D", main="Partial residual plots for top model, red line is loess curve")
		lines(seq(min(dfclean[,"z"]),max(dfclean[,"z"]),length.out=100),
			predict(loess(formula=Part_Res1~z,data=dfclean),
			newdata=data.frame(z=seq(min(dfclean[,"z"]),
			max(dfclean[,"z"]),length.out=100))), lwd=2, col="red")	
		dev.off()	
		
		# no major violations, as expected some of the super giants die the next  
		# year while the model predicts that they should have survived. Some of the XXS
		# plants also survive despite predicted mortality. Relationship still seems 
		# very clearly to be linear, if size is plotted on the log scale. 
		
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part E: Bootstrap to get approximations for estimates & standard errors to use in final model.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# tutorial from: R Data Analysis Examples: Mixed Effects Logistic Regression
	# Source: 
		#R Data Analysis Examples: Mixed Effects Logistic Regression.  UCLA: Statistical Consulting Group. 
		#from http://www.ats.ucla.edu/stat/r/dae/melogit.htm  (accessed March 30, 2015).
		
	# Due to random effect structure will have to rely on bootstrapped estimates for parameters of 
	# interest to account for variation across levels. 
	#############################################################################	
	# Sorry gang... please install & load all these libraries too :(
	library(lme4)
	require(ggplot2)
	require(GGally)
	require(reshape2)
	require(lme4)
	require(compiler)
	require(parallel)
	require(boot)	
		
	# difference across models so far 
		fixef(mod9)
		coef(mod.surv.glm)
		# Standard error for models, on the logit scale(?)
		se <- sqrt(diag(vcov(mod9))); (tab <- cbind(Est = fixef(mod9), LL = fixef(mod9) - 1.96 * se, UL = fixef(mod9) + 1.96 * se))
		# basic glm below 
		se <- sqrt(diag(vcov(mod.surv.glm))); (tab <- cbind(Est = coef(mod.surv.glm), LL = coef(mod.surv.glm) - 1.96 * se, UL = coef(mod.surv.glm) + 1.96 * se))
	# bootstrapping function
						sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
						cid <- unique(dat[, clustervar[1]])
						ncid <- length(cid)
						recid <- sample(cid, size = ncid * reps, replace = TRUE)
						if (replace) {
						rid <- lapply(seq_along(recid), function(i) {
						cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
						size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
						})
						} else {
						rid <- lapply(seq_along(recid), function(i) {
						cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
						})
						}
						dat <- as.data.frame(do.call(rbind, rid))
						dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
						labels = FALSE))
						dat$NewID <- factor(dat$NewID)
						return(dat)
						}
	# Resample data to a new df 100X larger  
		#set.seed(20)
		#	tmp <- sampler(dfclean, "PlotID", reps = 500) # run at 500 for now, but need to bump this up to 1000 or 10,000 once things get all set up
		#	bigdata <- cbind(tmp, dfclean[tmp$RowID, ])
		# LOAD SAME REPLICATE DATASET TO ESTIMATE CI'S FOR LAMBDA ESTIMATES  
		setwd(path.obj)
		#write.csv(bigdata, file="BootReplicateDataSets.csv", row.names=FALSE)
		bigdata <- read.csv(file="BootReplicateDataSets.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		# remove NA values 
		bigdata <- subset(bigdata, !is.na(Surv))
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		dim(bigdata)
	
	
#bigdata <- bigdata[which(bigdata$Replicate == "1" | bigdata$Replicate == "2" | bigdata$Replicate == "3"), ]
	
	# Steps: refit model to sample data
		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
		# 3. Export data & lme4 on each cluster.
		# 4. Write a fxn to fit the model and return estimates. 
		# * the call glmer() is wrapped in try to avoid stopping the process if models do not converge. 
		f <- fixef(mod9) # store original parameters 
		r <- getME(mod9, "theta") 
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Surv ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
					silent = TRUE)
					if (class(object) == "try-error")
					return(object)
					c(fixef(object), getME(object, "theta"))
				}
		# now move on to actually carrying out the bootstrapping, use 'parLapplyLB'  which loops through every replicate
		# giving them out to each node of the cluster to estimate models. 
		# results from all nodes are aggregated back into a single list & stored in object 'res'. 
		# R will then shut down local clusters. 
		start <- proc.time()
		res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
		end <- proc.time()
		# shut down the cluster
		stopCluster(cl)
	# Summarize results - calculate # of models that successfully converged.
		# Will do this be checking whether a particular result is numeric or not. 
		# calculate proportion of models that successfully converged
		success <- sapply(res, is.numeric)
		mean(success)
	# Convert list of bootstrap results into a matrix & 2.5th & 97.5th percentiles for each 
		# parameter for a table of results. Include original estimates and standard errors. 	
		# combine successful results
		bigres <- do.call(cbind, res[success])
		# calculate 2.5th and 97.5th percentiles for 95% CI
		(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))	
	# All results
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		# round and print
		round(finaltable, 3)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part F: Save final results for incorporation into IPM.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# results are fairly comparable, but possible to see a bit of bias. 
		# save results to (save final model & save table of boot strapped CI's 
		setwd(path.obj)
		write.csv(finaltable,"SurvModBootCoef.csv", row.names = T) 
		save(mod9, file = "SurvMod6.rda")
		save(res, file = "SurvModBootReplicates500.rda")

		
	##################################################################################			
	# FOR SITES - Steps: refit model to sample data
		bigdata <- read.csv(file="BootReplicateDataSets_Sites.csv") 
		bigdata$Year <- as.factor(bigdata$Year)
		# remove NA values 
		bigdata <- subset(bigdata, !is.na(Surv))
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		dim(bigdata)
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Surv ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
					silent = TRUE)
					if (class(object) == "try-error")
					return(object)
					c(fixef(object), getME(object, "theta"))
				}
		start <- proc.time()
		res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
		end <- proc.time()
		stopCluster(cl)
		success <- sapply(res, is.numeric)
		mean(success)
		bigres <- do.call(cbind, res[success])
		(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))	
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		round(finaltable, 3)		
		setwd(path.obj)
		write.csv(finaltable,"SurvModBootCoef_sites.csv", row.names = T) 
		save(res, file = "SurvModBootReplicates500_sites.rda")

	##################################################################################			
	# FOR PLOTS - Steps: refit model to sample data
		bigdata <- read.csv(file="BootReplicateDataSets_Plots.csv") 
		bigdata$Year <- as.factor(bigdata$Year)
		# remove NA values 
		bigdata <- subset(bigdata, !is.na(Surv))
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		dim(bigdata)
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Surv ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
					silent = TRUE)
					if (class(object) == "try-error")
					return(object)
					c(fixef(object), getME(object, "theta"))
				}
		start <- proc.time()
		res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
		end <- proc.time()
		stopCluster(cl)
		success <- sapply(res, is.numeric)
		mean(success)
		bigres <- do.call(cbind, res[success])
		(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))	
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		round(finaltable, 3)		
		setwd(path.obj)
		write.csv(finaltable,"SurvModBootCoef_plots.csv", row.names = T) 
		save(res, file = "SurvModBootReplicates500_plots.rda")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part G: Make Response Plots.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Make plots of results 
		# Graph out results to evaluate & check for other possible problems.  
		# conditional partial plots are fine, but run into the problem with MEMs of 
		# group level vs. average. Therefore we really want the average marginal probability.
		# will get the average marginal probability by plotting the average of each group - 
		# so will have to calculate the conditional probabilities for each group and then average them.
		summary(dfclean$z); summary(dfclean$z1) # values range from 0 - 10. 
		jvalues <- with(dfclean, seq(from = min(z), to = max(z), length.out = 100)) # sequence of values to plot.
		tmpdat <- dfclean # copy the object over as back up, predictions will be added to this data frame.  
		# calculate predicted probabilities and store in a list
			pp <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				predict(mod9, newdata = tmpdat, type = "response")
			})
		# Now we can display all the predicted probabilities & plot them 
		# average marginal predicted probability across a few different sizes z 
		round((sapply(pp[c(0.5, 1, 2, 4, 6, 8, 9, 10)], mean)*100), 1) # in percent
		# get the means with lower and upper quartiles
			plotdat <- t(sapply(pp, function(x) {
				c(M = mean(x), quantile(x, c(0.25, 0.75)))
			}))
		# add in z values and convert to data frame
		plotdat <- as.data.frame(cbind(plotdat, jvalues))
		# better names and show the first few rows
		colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "z")
		head(plotdat)
		setwd(path.fig)
		pdf(file="02.2.4_Survival_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_line() + ylim(c(0, 1))	
		# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))
		# refit model with SiteID & re-plot
		# convert SiteID down to just nine levels (the sites) 
		library(plyr)
		# calculate predicted probabilities and store in a list
		# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(dfclean$SiteID), function(stage) {
					  tmpdat$SiteID[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod9, newdata = tmpdat, type = "response")
					  })
					})
		# get means and quartiles for all jvalues for each level of SiteID
					plotdat2 <- lapply(biprobs, function(X) {
					  temp <- t(sapply(X, function(x) {
						c(M=mean(x), quantile(x, c(.25, .75)))
					  }))
					  temp <- as.data.frame(cbind(temp, jvalues))
					  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "z")
					  return(temp)
					})

		# collapse to one data frame
		plotdat2 <- do.call(rbind, plotdat2)
######## add SiteID
		plotdat2$SiteID <- factor(rep(levels(dfclean$SiteID), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for SITES
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SiteID), alpha = .15) +
					  geom_line(aes(colour = SiteID), size = 2) +
					  scale_fill_manual(values=c("green", "blue", "red", "yellow", "pink", "purple", "orange", "lightgrey")) +
					scale_colour_manual(values=c("green", "blue", "red", "yellow", "pink", "purple", "orange", "lightgrey"))+
					  ylim(c(0, 1))
					  
######### add Region
		# calculate predicted probabilities and store in a list
		rm(biprobs); rm(plotdat2)	
					
					biprobs <- lapply(levels(dfclean$Region), function(stage) {
					  tmpdat$Region[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod9, newdata = tmpdat, type = "response")
					  })
					})
			# get means and quartiles for all jvalues for each level of SiteID
					plotdat2 <- lapply(biprobs, function(X) {
					  temp <- t(sapply(X, function(x) {
						c(M=mean(x), quantile(x, c(.25, .75)))
					  }))
					  temp <- as.data.frame(cbind(temp, jvalues))
					  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "z")
					  return(temp)
					})	
		plotdat2 <- do.call(rbind, plotdat2)
		plotdat2$Region <- factor(rep(levels(dfclean$Region), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for region
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  scale_fill_manual(values=c("blue", "red")) +
					scale_colour_manual(values=c("blue", "red"))+
					  ylim(c(0, 1))
		

		dev.off()
	
	#############################################################################	
	
	

					
