###########################################################
# GROWTH FUNCTIONS FOR 2013-2015 CARDINALIS TRANSPLANT DATA
	# This script is run remotely from "02.1_basicIPM_Transplant.R"
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table on content:
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Part A: Exploratory plots across levels & factors from basic LM
			# Are there any major outliers ect. Was the log transformation of size sufficient? 
			
## Part B: Determine appropriate random effect structure & check model assumptions
			# with PlotID & SiteID as random effects. 

## Part C: Determine appropriate fixed structure (region, year, size, sites, interactions?)

## Part D: Bootstrap to get approximations for estimates & standard errors to use in final model.
			# are basic assumptions violated/any other problems? 

## Part E: Save final results for incorporation into IPM.

## Part F: Make Response Plots.



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

	# need to exclude the 2015-2016 flower/fecundity observations from dataset
	# in this script
	dim(dfclean)[1] # all plants left in Fall of 2014
	dfclean <- dfclean[!is.na(dfclean$z1),]	# exclude NA's from Flr & Fec 
	dim(dfclean) # number of plants that survived to final census July 201
	dfclean <- dfclean[order(dfclean$z),]	# sort smallest -> largest (z) 
	# ok we are now ready to get going



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part A: Exploratory plots across levels & factors from basic 'lm'
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Are there any major outliers ect. Was the log transformation of size sufficient? 
	# Will run exploratory plots for each site, year & region combination and examine visually. 
	# Combine plots into a large pdf. A bit tedious, but informative & useful before proceeding 
	# with the remainder of the analysis. 
	#############################################################################	
	dev.off()
	setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
	pdf(file="02.3.1_Growth_Explor.pdf", width=11, height=8.5)
		z <- dfclean$z # stand along object 
		# Plotting frame 
		# Growth (given survival), conditional on survival 
		grow.plot.data <- subset(dfclean, !is.na(z1))
		# Growth regression, simple linear - will do mixed model  with regions later
		mod.grow <- lm(z1 ~ z, data = grow.plot.data)
		palette("default")      # reset r color palette back to the default
		plot(z1 ~ z,
				 data = grow.plot.data,
				 xlim = plot.range, ylim = plot.range, pch = 21, cex = 1,
				 col=as.factor(grow.plot.data$SiteID),
				 bg=as.factor(grow.plot.data$PlotID),
				 xlab = expression("Initial size, "*italic(z)),
				 ylab = expression("Final size, "*italic(z)*"'"),
				 main="Global model (stem lengths log transformed)")
		# fitted line from growth regression 
		abline(mod.grow, col="red")
		abline(0, 1, lty=2) # line of zero growth (not positive or negative), useful for negative growth at large size classes. 
		add_panel_label(ltype="b")	
	# Run additional plots (optional code stored in R Code - functions subfolder)
	# additional plots in another script to save space on this sheet & make for efficent coding 
	setwd(path.funct)
	require(RCurl)
	source("02.3.2.extraplotsGrowth.R")
	# check for figure in figures folder - looks OK 
	# dev.off()
	# From exploratory work: On the log scale of total stem length, relationships
	# between t & t + 1 are defiantly linear. Variance structure might be an issue, 
	# but evidence of increasing variability was only really evident for some odd sites 
	# and was only supported by the odd outlier at a given site. 
	# Will explore below. But the growth model here is likely to remain pretty simple. 
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part B: Determine appropriate random effect structure & check model assumptions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Many random effect structures are possible for a mixed effect model 
	# Especially when there are both plot and site level nesting. 
	# It might even be possible to have a random slope component in addition 
	# to a random intercept.
	# Will use REML to evaluate alternative random effect structures. 
	# Source: Zuur et al 2009 Mixed Effect Models and Extension in Ecology with R
	# 	Chapter 5: Mixed Effects Modelling for Nested Data
			# Section: 5.8.2 The Good Approach (pages 127 - 142)
	
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
		# Start with as many fixed parameters as practically possible. 
			grow.plot.data <- subset(dfclean, !is.na(z1))
			grow.plot.data$ID <- as.factor(grow.plot.data$ID) # individual ID just to see what happens in an overly complex model  
			library(lme4)
			library(nlme) # for comparing random effect structures 
			Bbase <- gls(z1 ~ z*Region, method = "REML", data = grow.plot.data) # no nested factors
			Bplot <- lme(z1 ~ z*Region, data = grow.plot.data, random = ~1 | PlotID, method = "REML")
			Bsite <- lme(z1 ~ z*Region, data = grow.plot.data, random = ~1 | SiteID, method = "REML")
			Bboth <- lme(z1 ~ z*Region, data = grow.plot.data, random = ~1 | SiteID/PlotID, method = "REML")
			B_Stupidly_complex <- lme(z1 ~ z*Region, data = grow.plot.data, random = ~z | SiteID/PlotID, method = "REML")
			#B_Really_Stupidly_complex <- lmer(z1 ~ z*Region*Year + (z|SiteID:PlotID) + (1|ID), data = grow.plot.data)
			AIC(Bbase, Bplot, Bsite, Bboth)
			#AIC(B_Really_Stupidly_complex)
			anova(Bplot, Bsite)
			anova(Bplot, Bboth)
			anova(Bplot, Bbase)
			anova(Bsite, Bbase)
			
	#	UNEQUAL VARIANCE?
	#	We might expect that as size increases the variance 
	# 	in the individual growth rate might decrease & become more predictable. 
	#	If this is the case then we might want to specify a changing variance structure with 
	# 	size. This is unlikely for the cardinalis data as seen from the exploration plots.
	#	The log transformation of size might have already taken away the unequal variance
	#	problem for us, but we can do quick visual checks here & then specify an alternative 
	#	variance structure if needed. 
	
	# Check visually for unequal across variables: 
	setwd(path.fig); pdf(file="02.3.2_Growth_RandStructAssum.pdf", width=11, height=8.5)
	possible_rand_structs <- c("Bbase", "Bplot", "Bsite", "Bboth")
	for(i in 1:length(possible_rand_structs)){	
		this_mod <- get(possible_rand_structs[i])
		E2 <- resid(this_mod, type = "normalized")
		F2 <- fitted(this_mod)
		op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
		MyYlab <- "Residuals"
		plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab, main=paste0("Model: ", possible_rand_structs[i])) # weak/no evidence for unequal variance
		boxplot(E2 ~ SiteID, data = grow.plot.data, main = "SiteID", ylab = MyYlab) # weak/no evidence for unequal variance
		boxplot(E2 ~ Region, data =  grow.plot.data, main = "Region", ylab = MyYlab) # center region tends to be more variable
		plot(x =  grow.plot.data$z, y = E2, ylab = MyYlab, main = "Initial size", xlab = "log(initial size)") # weak/no evidence for unequal variance
		par(op)
	}
	# will revisit this section after adjusting the fixed effect structure,
	# but no strong evidence for alternative variance structure at this point. 
	# will chose a random effect structure of random site intercept for now
	# since it makes the most intuitive sense, fixed well with experimental design 
	# matched field observations and has enough data in each level (some plots occasional had only
	# a couple plants in them & replication of plots within site was also low ~3-5).
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part C: Determine appropriate fixed structure (region, year, size, sites, interactions?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Following:  Zuur et al. 2009: Mixed Effects Models and Extensions in Ecology with R, Springer)
		# Section: 5.10.7 Steps 7 and 8 of the Protocol: The Optimal Fixed Structure
		# pg. 135 
		summary(Bsite)
		anova(Bsite)
	# Drop superfluous 3-way interaction 
	# Switch from REML -> ML to compare fixed effects. 
		Form <- formula(z1 ~ z + SiteID + SiteID*z)
		M1.Full <- lme(Form, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		M1.A <- lme(z1 ~ z + SiteID, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		M1.B <- lme(z1 ~ z, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		anova(M1.Full, M1.A)
		anova(M1.Full, M1.B)
		anova(M1.A, M1.B)

		# should include site specific slope!
	# drop terms from above update model & try to drop additional terms
		Form2 <- formula(z1 ~ z + Region + Region:z)
		M2.Full <- lme(Form2, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		M2.B <- lme(z1 ~ z + Region, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		M2.C <- lme(z1 ~ z, random =~ 1 | PlotID, method = "ML", data = grow.plot.data)		
		anova(M2.Full, M2.B)
		anova(M2.Full, M2.C)
		anova(M2.B, M2.C)
		# region doesn't really seem to be doing anything here...
		anova(M1.Full, M2.C)

		
	# 5.10.8 Step 9 of the Protocol: Refit with REML and Validate the Model
		M3 <- lme(z1 ~ z + SiteID + SiteID*z, random =~1|PlotID,
				method = "REML", data = grow.plot.data)
		summary(M3)
		methods(class = "lme")
		denom <- as.numeric(VarCorr(M3)[1,1]) + as.numeric(VarCorr(M3)[2,1])
		# Variance due to SiteID = 0.09, PlotID = 0.12
		as.numeric(VarCorr(M3)[1,1])/denom
	
	# 5.10.10 pg. 139 
		# is the interactive effect of year & region really inconsequential? 
		# fit with a smoother spline to visualize
		# also good to check to make sure the effect is indeed linear
		library(lattice)
		xyplot(z1 ~ z|Region*Year,
			data = grow.plot.data, ylab = "Residuals",
			xlab = "Initial size (log)",
			panel = function(x,y){
			panel.grid(h = -1, v = 2)
			panel.points(x, y, col = 1)
			panel.loess(x, y, span = 0.5, col = 1,lwd=2)})
		# its pretty darn linear on the log scale		
		dev.off()
		
		# fit a gamm just to make extra sure 
		library(mgcv)
		M6 <- gamm(z1 ~ SiteID + s(z),
			random = list(PlotID = ~ 1), data = grow.plot.data)
		anova(M6$gam)
		plot(M6$gam) # looks pretty linear 
		#plot(M6$lme)
	
##################################################################################
	# The MuMIn Package tutorial: 
	# Look at estimates with effect sizes & AIC scores for competing models 
	# package for getting effects & AIC values 
	# Code from tutorial available at  https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	
	require(MuMIn)
	#################
	# MODEL SELECTION 
		# use the mod.sel function to conduct model selection
		# and put output into object out.put
		# rename models from previous section to be a bit more useful 
			
		M1.Full <- M1.Full; M1A.int <- M1.A; M1B.intslope <- M1.B; M2C.justsize <- M2.C
		out_put <- mod.sel(M1.Full, M1A.int, M1B.intslope, M2C.justsize)
		modList <- c("M1.Full", "M1A.int", "M1B.intslope", "M2C.justsize")
		out_put

	# TOP MODELS
		# OPTIONAL - select models 95% cumulative weight criteria 
		# IMPORTANT: Weights have been renormalized!! 
		subset(out_put, cumsum(out_put$weight) <= .95) # just the full model really
		
	# MAKE CLEAN AICc TABLE & SAVE
		# coerce the object out_put into a data frame 
		# elements 6-10 in out_put have what we want 
		sel.table <- as.data.frame(out_put)		
		# a little clean-up, lets round things a bit 
		sel.table[,c(1,3,6:9)]<- round(sel.table[,c(1,3,6:9)],2) 
		# that’s better 
		sel.table 
		## lets be sure to put the model names in a column 
		sel.table$Model<-rownames(sel.table); sel.table
		# replace Model name with formulas little tricky so be careful 
		for(i in 1:nrow(sel.table)) 
			sel.table$Model[i]<- as.character(formula(paste(sel.table$Model[i])))[3] 
		# let's see what is in there 
		sel.table 
		#little reordering of columns 
		sel.table<-sel.table[,c(10, 5, 6, 7, 8, 9)]
		sel.table
		# write to a file, here a comma separated values format 
		# make sure your working directory is properly specified 
		setwd(path.obj)
		write.csv(sel.table,"GrowModAICTable.csv", row.names = T) 
			
	# VARIABLE IMPORTANCE (BASED ON AIC)	
		importance(out_put)
	
	# MODEL AVERAGING (over top models) 
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
			#write.csv(MA.est.table, "GrowModAvgEst.csv") 
			
	# MAKE PREDICTIONS	
			# extract parameters and weights from confidence model set 
			# using get.models function 
			pred.parms<-get.models(out_put) 
			# predict values using each model, here were just using the 
			# the example dataset, you could use a new dataset 
			model.preds = sapply(pred.parms, predict, newdata = dfclean) 			

	# Weight predictions by model akaki weights
		# We can weigh the predictions (by the AIC weights) for each of the top models
		# & average across
		# we also are using matrix multiplication %*% 
		mod.ave.preds<-model.preds %*% Weights(out_put) 
		head(mod.ave.preds) # just to see

	# Make partial plots - won't really work with categorical predictors 
		# size ranges from the observed minimum to maximum 
		z_range=c(min(grow.plot.data$z):max(grow.plot.data$z)) 
		# create plotdata data frame with mean values 
		#plotdata<-as.data.frame(lapply(lapply(grow.plot.data[5:8],mean),rep,length(z_range)))) 
		#plotdata<-cbind(z_range,plotdata) 
		# now predict density for the plot data with each model 
		#model.preds = sapply(pred.parms, predict, newdata = plotdata) 
		# weight the prediction from each model by its AIC weight 
		# and sum (matrix multiplication) 
		#mod.ave4plot<-model.preds %*% Weights(out.put) 
		# plot the model averaged predicted densities vs elevation 
		#plot(mod.ave4plot~ elev, type = 'l', xlab="Elevation (m)", ylab="Model averaged predicted density") 

		

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part D: Bootstrap to get approximations for estimates & standard errors to use in final model.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# exclude year from final model 
	library(lme4)
	library(lmerTest)
	# refit with lmer
	M2.A <- lmer(z1 ~ z + SiteID + SiteID:z + (1|PlotID), data = grow.plot.data)		

	# bootMer from lme4 is similar to bootCase from car (but poorly documented)
	# Make sure to set re.form=NA to output only predictions using the fixed effects.
	bootfit2 <- bootMer(M2.A, FUN=function(x)predict(x, grow.plot.data, re.form=NA),
							nsim=999)

	# We have to find that the $t element contains the predictions output by bootMer.
	head(apply(bootfit2$t, 2, sd))
	bb <- bootfit2 

####################################################################################	
# load bootstrapping function 
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

# resample our data with 1000 replicates													
		set.seed(20)
		#tmp <- sampler(grow.plot.data, "SiteID", reps = 500) # shold be in the 1000's
		#bigdata <- cbind(tmp, grow.plot.data[tmp$RowID, ])
		# LOAD SAME REPLICATE DATASET TO ESTIMATE CI'S FOR LAMBDA ESTIMATES  
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		#toMerge <- bigdata[,c("Replicate", "ID")]
		#toMerge$ID <- as.factor(toMerge$ID)
		#toMerge$Replicate <- as.factor(toMerge$Replicate)
		#bigdata <- merge(toMerge, grow.plot.data, by="ID", all.x=TRUE)
		# remove NA values 
		bigdata <- subset(bigdata, !is.na(z1))
		dim(bigdata)
	
		
	require(ggplot2)
	require(GGally)
	require(reshape2)
	require(lme4)
	require(compiler)
	require(parallel)
	require(boot)
	library(lme4)
	# Next we refit the model on the resampled data. 
	#First we store the estimates from our original model, 
	#which we will use as start values for the bootstrap models.
	# Then we make a local cluster with 4 nodes (the number of processors 
	#on our machine; set to the number of processors you have on yours). 
	#Next, we export the data and load the lme4 package on the cluster. 
	#Finally, we write a function to fit the model and return the estimates.
	# The call to lmer() is wrapped in try because not all models may 
	#converge on the resampled data. This catches the error and returns 
	#it, rather than stopping processing.

	f <- fixef(M2.A)
	r <- getME(M2.A, "theta")

	cl <- makeCluster(4)
	clusterExport(cl, c("bigdata", "f", "r"))
	clusterEvalQ(cl, require(lme4))	
		
	myboot <- function(i) {
		object <- try(lmer(z1 ~ z + SiteID + SiteID:z + (1|PlotID), na.action=na.omit, data = bigdata,
			subset = Replicate == i), silent = TRUE)
		if (class(object) == "try-error")
			return(object)
		c(fixef(object), getME(object, "theta"))
	}	
	start <- proc.time()
	res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
	end <- proc.time()
	
	# shut down the cluster
	stopCluster(cl)	
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
		# Standard error for models
		se <- sqrt(diag(vcov(M2.A)))
		(tab <- cbind(Est = fixef(M2.A), LL = fixef(M2.A) - 1.96 * se, UL = fixef(M2.A) + 1.96 * se))
		# basic lm below 
		se <- sqrt(diag(vcov(mod.grow))); (tab <- cbind(Est = coef(mod.grow), LL = coef(mod.grow) - 1.96 * se, UL = coef(mod.grow) + 1.96 * se))
		
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		# round and print
		round(finaltable, 3)

		
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part E: Save final results for incorporation into IPM.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		# results are fairly comparable, but possible to see a bit of bias. 
		# save results to (save final model & save table of boot strapped CI's 
		setwd(path.obj)
		write.csv(finaltable,"GrowModBootCoef.csv", row.names = T) 
		save(M2.A, file = "GrowMod.rda")
		save(res, file = "GrowModBootReplicates500.rda")


	
	##################################################################################			
	# FOR SITES - Steps: refit model to sample data
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets_Sites.csv") 
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		bigdata <- subset(bigdata, !is.na(z1))
		dim(bigdata)
		f <- fixef(M2.A)
		r <- getME(M2.A, "theta")
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))			
			myboot <- function(i) {
				object <- try(lmer(z1 ~ z + SiteID + SiteID:z + (1|PlotID), na.action=na.omit, data = bigdata,
					subset = Replicate == i), silent = TRUE)
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
		se <- sqrt(diag(vcov(M2.A)))
		(tab <- cbind(Est = fixef(M2.A), LL = fixef(M2.A) - 1.96 * se, UL = fixef(M2.A) + 1.96 * se))
		se <- sqrt(diag(vcov(mod.grow))); (tab <- cbind(Est = coef(mod.grow), LL = coef(mod.grow) - 1.96 * se, UL = coef(mod.grow) + 1.96 * se))	
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		round(finaltable, 3)
		write.csv(finaltable,"GrowModBootCoef_sites.csv", row.names = T) 
		save(res, file = "GrowModBootReplicates500_sites.rda")


	##################################################################################			
	# FOR PLOTS - Steps: refit model to sample data
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets_Plots.csv") 
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		bigdata <- subset(bigdata, !is.na(z1))
		dim(bigdata)
		f <- fixef(M2.A)
		r <- getME(M2.A, "theta")
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))			
			myboot <- function(i) {
				object <- try(lmer(z1 ~ z + SiteID + SiteID:z + (1|PlotID), na.action=na.omit, data = bigdata,
					subset = Replicate == i), silent = TRUE)
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
		se <- sqrt(diag(vcov(M2.A)))
		(tab <- cbind(Est = fixef(M2.A), LL = fixef(M2.A) - 1.96 * se, UL = fixef(M2.A) + 1.96 * se))
		se <- sqrt(diag(vcov(mod.grow))); (tab <- cbind(Est = coef(mod.grow), LL = coef(mod.grow) - 1.96 * se, UL = coef(mod.grow) + 1.96 * se))	
		finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci)
		round(finaltable, 3)
		write.csv(finaltable,"GrowModBootCoef_plots.csv", row.names = T) 
		save(res, file = "GrowModBootReplicates500_plots.rda")


	

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part F: Make Response Plots.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Make plots of results 
		# Graph out results to evaluate & check for other possible problems.  
		# conditional partial plots are fine, but run into the problem with MEMs of 
		# group level vs. average. Therefore we really want the average marginal probability.
		# will get the average marginal probability by plotting the average of each group - 
		# so will have to calculate the conditional probabilities for each group and then average them.
		summary(grow.plot.data$z); summary(grow.plot.data$z1) # values range from 0 - 10. 
		jvalues <- with(grow.plot.data, seq(from = min(z), to = max(z), length.out = 100)) # sequence of values to plot.
		tmpdat <- grow.plot.data # copy the object over as back up, predictions will be added to this data frame.  
		# calculate predicted probabilities and store in a list
			pp <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				predict(M2.A, newdata = tmpdat)
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
		colnames(plotdat) <- c("PredictedSize", "Lower", "Upper", "z")
		head(plotdat)
		setwd(path.fig)
		pdf(file="02.3.3_Growth_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedSize)) + geom_line() + ylim(c(0, 8.053))	
		#We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedSize)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 1) + ylim(c(0, 8.053))
		# refit model with year & region & re-plot
		# convert regions down to just three levels WITHIN/BEYOND
		library(plyr)
		# calculate predicted sizes and store in a list
					biprobs <- lapply(levels(grow.plot.data$Region), function(stage) {
					  tmpdat$Region[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(M2.A, newdata = tmpdat)
					  })
					})
		# get means and quartiles for all jvalues for each level of Region
					plotdat2 <- lapply(biprobs, function(X) {
					  temp <- t(sapply(X, function(x) {
						c(M=mean(x), quantile(x, c(.25, .75)))
					  }))
					  temp <- as.data.frame(cbind(temp, jvalues))
					  colnames(temp) <- c("PredictedSize", "Lower", "Upper", "z")
					  return(temp)
					})

		# collapse to one data frame
		plotdat2 <- do.call(rbind, plotdat2)
		# add Region
		plotdat2$Region <- factor(rep(levels(grow.plot.data$Region), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for N, C, S
		ggplot(plotdat2, aes(x = z, y = PredictedSize)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  scale_fill_manual(values=c("green", "blue", "red")) +
					scale_colour_manual(values=c("green", "blue", "red"))+
					  ylim(c(0, 8.053)) + xlim(c(0, 8.053)) 

		# can view distribution of predicted sizes here:
	
	#############################################################################	
	# convert regions down to just site levels 
		library(plyr)
		# calculate predicted sizes and store in a list
					biprobs <- lapply(levels(grow.plot.data$SiteID), function(stage) {
					  tmpdat$SiteID[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(M2.A, newdata = tmpdat)
					  })
					})
		# get means and quartiles for all jvalues for each level of SiteID
					plotdat2 <- lapply(biprobs, function(X) {
					  temp <- t(sapply(X, function(x) {
						c(M=mean(x), quantile(x, c(.25, .75)))
					  }))
					  temp <- as.data.frame(cbind(temp, jvalues))
					  colnames(temp) <- c("PredictedSize", "Lower", "Upper", "z")
					  return(temp)
					})

		# collapse to one data frame
		plotdat2 <- do.call(rbind, plotdat2)
		# add SiteID
		plotdat2$SiteID <- factor(rep(levels(grow.plot.data$SiteID), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for SiteID
			# Will use Wes Anderson color palettes because he is simply the best
			library(RColorBrewer)
			#display.brewer.all()
			n = length(levels(plotdat2$SiteID))
			brewer.pal(n, "Set1")
		
		ggplot(plotdat2, aes(x = z, y = PredictedSize)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SiteID), alpha = .15) +
					  geom_line(aes(colour = SiteID), size = 2) +
					  scale_fill_manual(values=brewer.pal(n, "Set1")) +
					scale_colour_manual(values=brewer.pal(n, "Set1"))+
					  ylim(c(0, 8.053)) + xlim(c(0, 8.053)) 

		# can view distribution of predicted sizes here:
		
		dev.off()
	





	
