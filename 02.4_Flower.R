###########################################################
# FLOWERING FUNCTION FOR 2013 - 2015 TRANSPLANT DATA

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table on content:
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# A. Exploratory plots across levels & factors from basic GLM
			# Are there any major outliers ect. Was the log transformation of size sufficent? 
			
	# B. Determine appropriate random effect structure (single random effect of sites vs. plots within sites ect.) 
			# for glmer with PlotID & SiteID as random effects. 
			
	# C. Determine appropriate fixed structure (region, year, size, sites, interactions?)
			# for flowering glmer with a Akaike-weights framework. 
			
	# D. Make post-hoc diagnostic plots 
			# are basic assumptions violated/any other problems? 
			
	# E. Bootstrap to get approximations for estimates & standard errors to use in final model. 
	
	# F. Save final results for incorporation into IPM. 

	# G. Make Response Plots.

	
	
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

	# need to include all flowering data in this script 2014/2015 & 2015/2016
	dim(dfclean)[1] # all plants in dataset
	dfclean <- dfclean[!is.na(dfclean$Repr),]	# exclude NA's from Flr & Fec 
	dim(dfclean)[1] # all plants left with flowering data
	# ok we are now ready to get going
	
	
	
	
	
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
##############################
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1A - Load vital rate regression models & coefficients from scripts 04.2.1 - 04.2.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## sort by size and print a sample to the screen
	dfclean <- dfclean[order(dfclean$z),]					# sort smallest -> largest (z) 
	

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part A: Exploratory plots across levels & factors from basic GLM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Are there any major outliers ect. Was the log transformation of size sufficient? 
	# Will run exploratory plots for each site, year & region combination and examine visually. 
	# Combine plots into a large pdf. A bit tedious, but informative & useful before proceeding 
	# with the remainder of the analysis. 
	#############################################################################	
			dev.off()
			setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
			pdf(file="02.4.1_Flower_Explor.pdf", width=11, height=8.5)
			Repr.plot.data <- subset(dfclean, !is.na(Repr))
			z <- Repr.plot.data$z # stand alone object 
				# Plotting frame 
					Repr.plot.data <- within(Repr.plot.data, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global flower model
			mod.Repr.glm <- glm(Repr ~ z, family = binomial, data = Repr.plot.data)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.data)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of Flowering", main="Global model")
			points(z, jitter(Repr.plot.data$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = Repr.plot.data, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
	
		# Run additional plots (optional code stored in R Code - functions subfolder)
		# additional plots in another script to save space on this sheet & make for efficent coding 
			setwd(path.funct)
			source("02.4.extraplotsRepr.R")
			# check for figure in figures folder - looks OK 
			

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part B: Determine appropriate random effect structure
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Many random effect structures are possible for a mixed effect model 
	# Especially when there are both plot and site level nesting. 
	# It might even be possible to have a random slope component in addition 
	# to a random intercept. 
	# Unfortunately defendable statistical methods to choose an optimal 
	# random effect structure for glm objects are still under development. 
	# discussion: http://www.researchgate.net/post/How_can_I_optimize_the_random_effect_structure_in_a_GLMM
	# good discussion of issue here too: http://glmm.wikidot.com/faq
	# REML procedure not possible here for glm objects. Should rely on own 
	# intuition, experimental setup and a solid hypothesis for each alternative 
	# random effect structure. Some approaches are covered here. 
	
		# Possible Alternatives (listed from most believable to least believable): 
				#A - (1|PlotID)	
				#B - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
				#C - (z|PlotID) = (1+z|PlotID)	random slope of z within PlotID with correlated intercept for site
				#D - (z|site/PlotID) = (z|site)+(z|site:PlotID) random slope of z within both site and plot. intercept varying among sites and among PlotIDs within sites (nested random effects)
				#E - site+(1|site:PlotID)	fixed effect of sites plus random variation in intercept among PlotIDs within sites, unlikely & difficult to work with. Too many site levels. 				
			# Not considered: 
				#F - (0+z|SiteID) = (-1+z|SiteID)	random slope of z within SiteID: no variation in intercept
				#G - (1|SiteID) + (0+z|SiteID)	uncorrelated random intercept and random slope within SiteID
	#############################################################################	
	# To avoid over-fitting will use cross-validation with successively smaller sub-samples of
	# the data (training & testing sets) to compare the complexity of alternative 
	# random effect structures. If training and testing sets are cut randomly from the 
	# full data frame 
		library(lme4)
		library(arm)
		library(sjPlot)
		
		
		holdOutPortions <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.95)
		myFrame <- data.frame()
		for(i in 1:length(holdOutPortions)){
			# Function to split data into training & testing set 
				splitdf <- function(dataframe, seed=NULL) {
				if (!is.null(seed)) set.seed(seed)
				index <- 1:nrow(dataframe)
				trainindex <- sample(index, trunc(length(index)*holdOutPortions[i]))
				trainset <- dataframe[trainindex, ]
				testset <- dataframe[-trainindex, ]
				list(trainset=trainset,testset=testset)
			}
		
		splits <- splitdf(Repr.plot.data, seed=808) # splits dataframe 
		# str(splits) # training & testing set
		# Observations in each data frame
		lapply(splits,nrow)
		#view the first few columns in each data frame
		lapply(splits,head)
		# save the training and testing sets as data frames
		training <- splits$trainset
		testing <- splits$testset
		
		# To avoid having unique PlotID levels in only one of training 
		# and testing sets add one observation from each plot (lowest level) to each dataset
			df2 <- lapply(split(Repr.plot.data, Repr.plot.data$PlotID), function(subdf) subdf[sample(1:nrow(subdf), 1, replace=TRUE),])		
			df2 <- do.call('rbind', df2)
			head(df2, 3); dim(df2)
			training <- rbind(training, df2)
			testing <- rbind(testing, df2)

	# Run glmer with all potential fixed effect terms (no three way interaction)
		testA_plot <- glmer(Repr ~ z + Year + Region + z*Region + (1|PlotID), na.action=na.omit, data=training, family = binomial)
		testA_site <- glmer(Repr ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=training, family = binomial)
		testA_both <- glmer(Repr ~ z + Year + Region + z*Region + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=training, family = binomial)
		testA_neither <- glm(Repr ~ z + Year + Region + z*Region, na.action=na.omit, data=training, family = binomial)

		
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

			
	modPlot <- roc(Repr ~ ModPlotPreds, data = testing); PlotA <- modPlot$auc[1]
	modSite <- roc(Repr ~ ModSitePreds, data = testing); SiteA <- modSite$auc[1]
	modBoth <- roc(Repr ~ ModBothPreds, data = testing); BothA <- modBoth$auc[1]
	modNeith <- roc(Repr ~ ModNeithPreds, data = testing); NeithA <- modNeith$auc[1]
	holdout <- holdOutPortions[i]
	N <- dim(training)[1]*holdOutPortions[i]
	
	toAdd <- data.frame(holdout, N, PlotA, SiteA, BothA, NeithA)
	
	myFrame	<- rbind(myFrame, toAdd)
	}
	
	# save final evaluation figure
	setwd(path.fig)
	pdf(file="02.4.2_Flower Random Effect Structure.pdf", width=11, height=8.5)
	colnames(myFrame) <- c("trainPor", "N", "PlotA", "SiteA", "BothA", "NeithA") 
	plot(myFrame$trainPor, myFrame$PlotA, xlab="Portion of data Training:Testing ratio", 
		ylab="AUC", col=1, ylim=c(min(myFrame[,3:6]) - 0.04, max(myFrame[,3:6]) + 0.04), pch=19, cex=1.5, main="Flowering Random Effect Structure")
	lines(myFrame$trainPor, myFrame$PlotA, col=1)
	points(myFrame$trainPor, myFrame$SiteA, col=2, pch=19); lines(myFrame$trainPor, myFrame$SiteA, col=2)
	points(myFrame$trainPor, myFrame$BothA, col=3, pch=19); lines(myFrame$trainPor, myFrame$BothA, col=3)
	points(myFrame$trainPor, myFrame$NeithA, col="darkgrey", pch=19); lines(myFrame$trainPor, myFrame$NeithA, col="darkgrey")
	legend('topright', legend = c("(1|PlotID)", "(1|SiteID)", "(1|PlotID) + (1|SiteID) + (1|SiteID:PlotID)", "glm"), col = c(1:3, "darkgrey"), cex = 0.8, pch = 19)
	text(myFrame$trainPor, min(myFrame[,3:5]) - 0.02, labels = round(myFrame$N, 0))
	
	# plot of coefficients
		# https://strengejacke.wordpress.com/sjplot-r-package/
		
		# plot fixed effects
		par(mfrow=c(1,3))
		#sjp.glmer(testA_plot, type = "fe")
		#sjp.glmer(testA_site, type = "fe")
		#sjp.glmer(testA_both, type = "fe")
	
	dev.off()
	
	# Having plot alone seems to work better than having site alone, but these results don't match experimental design or make intuative sence 
	# also very low sample sizes for a large number of plots. 
	# Will have to reassess and see how sensitivity final results are to estimates from each random 
	# effect structure. Site should remain as a random effect because that was more
	# geared to how the study was designed. Intra-plot variability should be high 
	# because of how transects stretched over knolls, depressions, banks and streams. 
	# also replication within plots is occasionally very low (some only have few plant observations). 
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part C: Determine appropriate fixed structure (region, year, size, sites, interactions?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# For flowering glmer with a Akaike-weights framework. 
	# Largest possible full model:
		#"Repr ~ z + SiteID + Region + z*SiteID + z*Region"
	# Procedure followed here from: 
		# Chapter 21: GLMM applied on the spatial distribution of Koalas in a fragmented landscape
			# section 24.4.2 (pg. 483) Zuur 2009 - "Mixed effect models and extension in ecology with R"
			# code from: https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	#############################################################################	
	
	# Possible models(yes could be coded more efficiently)  
	#mod1 <- glmer(Repr ~ z + SiteID + Region + z*SiteID + z*Region + z*SiteID*Region + (1|PlotID) + (1|PlotID) + (1|PlotID:PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod2 <- glmer(Repr ~ z + I(z^2) + SiteID + Region + z*SiteID + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod3 <- glmer(Repr ~ z + I(z^2) + SiteID + Region + z*SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod4 <- glmer(Repr ~ z + I(z^2) + SiteID + Region + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod5 <- glmer(Repr ~ z + I(z^2) + SiteID + Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod6 <- glmer(Repr ~ z + I(z^2) + Region + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod7 <- glmer(Repr ~ z + I(z^2) + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod8 <- glmer(Repr ~ z + I(z^2) + Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod9 <- glmer(Repr ~ z + I(z^2) + SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod10 <- glmer(Repr ~ z + I(z^2) + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod2linear <- glmer(Repr ~ z + SiteID + Region + z*SiteID + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod3linear <- glmer(Repr ~ z + SiteID + Region + z*SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod4linear <- glmer(Repr ~ z + SiteID + Region + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	#mod5linear <- glmer(Repr ~ z + SiteID + Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod6linear <- glmer(Repr ~ z + Region + z*Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod7linear <- glmer(Repr ~ z + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod8linear <- glmer(Repr ~ z + Region + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod9linear <- glmer(Repr ~ z + SiteID + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)
	mod10linear <- glmer(Repr ~ z + (1|PlotID), na.action=na.omit, data=Repr.plot.data, family = binomial)

	# package for getting effects & AIC values 
	require(MuMIn)
	
	#################
	# MODEL SELECTION 
	#################
	# which models worked?

		out_put <- mod.sel(mod6,mod7,mod8,mod9,mod10,mod6linear,mod7linear,mod8linear,mod9linear,mod10linear)
		modList <- c("mod6","mod7","mod8","mod9","mod10", "mod6linear","mod7linear","mod8linear","mod9linear","mod10linear")
		out_put
	# Don't include quadratic term
	
		# use the mod.sel function to conduct model selection
		# and put output into object out.put
		out_put <- mod.sel(mod6linear,mod7linear,mod8linear,mod9linear,mod10linear)
		modList <- c("mod6linear","mod7linear","mod8linear","mod9linear","mod10linear")
		out_put

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
		write.csv(sel.table,"ReprModAICTable.csv", row.names = T) 
			
	# VARIABLE IMPORTANCE (BASED ON AIC)	
		importance(out_put)
		# very strong year effect, much greater than region. 
	
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
			#write.csv(MA.est.table, "ReprModAvgEst.csv") 
			# Unfortunately estimating response values from these average 
			# parameters won't work for logistic regression.  			
			
	# MAKE PREDICTIONS	
			# extract parameters and weights from confidence model set 
			# using get.models function 
		# pred.parms<-get.models(out_put, subset= delta < 5) 
			# predict values using each model, here were just using the 
			# the example dataset, you could use a new dataset 
		# model.preds = sapply(pred.parms, predict(type="response"), newdata = Repr.plot.data) 			
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
			overdisp.glmer(mod9linear) # seems ok, a bit under-dispersed 
			overdisp.glmer(mod8linear) # seems ok, a bit under-dispersed 


		# winning model:
			summary(mod9linear)
		
		library(glmmML) # run sample model with glmmML - do coefficents change?
		mod9linear_glmmML <-glmmML(formula=Repr ~ z + I(z^2) + Year + z*Year + Region,
             cluster=SiteID,data=Repr.plot.data,family=binomial)
		summary(mod9linear_glmmML)
		fixef(mod9linear)

	
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
	# load in model adequacy functions from Zuur. 
	setwd(path.funct)
	source("ZuurGLMMAdequecyFxns.R") # downloaded from: http://www.highstat.com/book2.htm
	
	# QUANTILE-QUANTILE PLOTS (note that the functions "fitted.glmmML" and "res.glmmML" are requiured here, loaded in line above)	
		#calculate the fitted values for the model
		Model_Best <- mod9linear_glmmML # best model from steps above
		Fitted <- fitted.glmmML(model=Model_Best, data=Repr.plot.data)
		#calculate the ordered residuals for the model
		Resids<-sort(res.glmmML(model=Model_Best,data=Repr.plot.data))
		Resids_0<-Resids[Resids<0]
		Resids_1<-Resids[Resids>=0]
		#specify the number of replicates
		Reps<-500
		#simulate Reps data sets
		Sims<-matrix(rbinom(n=length(Fitted)*Reps,size=1,p=Fitted),nrow=length(Fitted),ncol=Reps)
		#fit the model to each simulated data set
		Models<-apply(Sims,MARGIN=2,FUN=function(X){return(list(X,glmmML(formula=Repr ~ z + SiteID,cluster=PlotID,data=Repr.plot.data,family=binomial)))})
		#calculate the ordered (interpolated) simulated residuals and the point-wise 95% confidence intervals
		Resids_Sim<-matrix(unlist(lapply(Models,FUN=function(X){TempData<-Repr.plot.data;TempData[,"Repr"]<-X[[1]];return(sort(res.glmmML(X[[2]],TempData)))})),ncol=Reps,nrow=nrow(Repr.plot.data))
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
		dev.off()
		setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
		pdf(file="02.4.3_Flower_model_adequacy.pdf", width=11, height=8.5)
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
		Part_Res1<-((res.glmmML(Model_Best,Repr.plot.data)-fitted.glmmML(Model_Best,Repr.plot.data))/
		((fitted.glmmML(Model_Best,Repr.plot.data)*(1-fitted.glmmML(Model_Best,Repr.plot.data)))))+
		(Model_Best$coefficients[2]*Repr.plot.data[,"z"])	
		#plot the partial residuals and the smoothed plot
		plot(Repr.plot.data[,"z"],Part_Res1,ylab="partial residuals", cex=0.5, pch=16, 
			col= "#0000FF4D", main="Partial residual plots for top model, red line is loess curve")
		lines(seq(min(Repr.plot.data[,"z"]),max(Repr.plot.data[,"z"]),length.out=100),
			predict(loess(formula=Part_Res1~z,data=Repr.plot.data),
			newdata=data.frame(z=seq(min(Repr.plot.data[,"z"]),
			max(Repr.plot.data[,"z"]),length.out=100))), lwd=2, col="red")	
		
		dev.off()	
		
# Comment: "Things seem OK from the quantile-quantile plot, but 
			# the partial residual plots raises quite a bit of 
			# concern. Notice that the loess fit to residuals is not 
			# strait at all, but actually curves down sharply at the end. 
			# This suggests that a linear assumption in the glmm of the initial size
			# z might not be the most appropriate model for the data. 
			# Should a quadratic term for size be included in the model?
		
		# Refit with quadratic term: 
		mod9linear_glmmML_NoQuadratic <-glmmML(formula=Repr ~ z + SiteID, cluster=PlotID,data=Repr.plot.data,family=binomial)
		mod9linear_glmmML_Quadratic <-glmmML(formula=Repr ~ z + I(z^2) + SiteID, cluster=PlotID,data=Repr.plot.data,family=binomial)
		summary(mod9linear_glmmML_Quadratic)
		AIC(mod9linear_glmmML_NoQuadratic) # lower AIC from model with quadratic term ... 
		# but maybe this is just a regional effect 
		Quadratic2NoRegion <-glmmML(formula=Repr ~ z + I(z^2) + SiteID, cluster=PlotID,data=Repr.plot.data,family=binomial)
		AIC(Quadratic2NoRegion)
		# but aic model selection above doesn't jibe.. ?

	
	
	
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
		fixef(mod9linear)
		coef(mod.Repr.glm)
		# Standard error for models, on the logit scale(?)
		se <- sqrt(diag(vcov(mod9linear))); (tab <- cbind(Est = fixef(mod9linear), LL = fixef(mod9linear) - 1.96 * se, UL = fixef(mod9linear) + 1.96 * se))
		# basic glm below 
		se <- sqrt(diag(vcov(mod.Repr.glm))); (tab <- cbind(Est = coef(mod.Repr.glm), LL = coef(mod.Repr.glm) - 1.96 * se, UL = coef(mod.Repr.glm) + 1.96 * se))
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
	# Resample data to a new df 500X larger  
		#set.seed(20)
		#tmp <- sampler(Repr.plot.data, "SiteID", reps = 500) # run at 500 for now, but need to bump this up to 1000 or 10,000 once things get all set up
		#bigdata <- cbind(tmp, Repr.plot.data[tmp$RowID, ])
		# LOAD SAME REPLICATE DATASET TO ESTIMATE CI'S FOR LAMBDA ESTIMATES  
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		bigdata <- bigdata[!is.na(bigdata$Repr),]	# exclude NA's from Flr & Fec 

	
	# Steps: refit model to sample data
		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
		# 3. Export data & lme4 on each cluster.
		# 4. Write a fxn to fit the model and return estimates. 
		# * the call glmer() is wrapped in try to avoid stopping the process if models do not converge. 
		f <- fixef(mod9linear) # store original parameters 
		r <- getME(mod9linear, "theta") 
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Repr ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
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
		write.csv(finaltable,"ReprModBootCoef.csv", row.names = T) 
		save(mod9linear, file = "Reprmod9L.rda")
		
		# save boot replicate datasets
		# write.csv(bigdata,"BootReplicateDataSets_Repr.csv", row.names = T) 
		# save boot replicate coefficents
		save(res, file = "ReprModBootReplicates500.rda")

			
	##################################################################################			
	# FOR SITES - Steps: refit model to sample data		
		bigdata <- read.csv(file="BootReplicateDataSets_Sites.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		bigdata <- bigdata[!is.na(bigdata$Repr),]	# exclude NA's from Flr & Fec 
		f <- fixef(mod9linear) # store original parameters 
		r <- getME(mod9linear, "theta") 
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Repr ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
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
		write.csv(finaltable,"ReprModBootCoef_sites.csv", row.names = T) 
		save(res, file = "ReprModBootReplicates500_sites.rda")

		
	##################################################################################			
	# FOR PLOTS - Steps: refit model to sample data		
		bigdata <- read.csv(file="BootReplicateDataSets_Plots.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		bigdata$Replicate <- as.factor(bigdata$Replicate)
		bigdata <- bigdata[!is.na(bigdata$Repr),]	# exclude NA's from Flr & Fec 
		f <- fixef(mod9linear) # store original parameters 
		r <- getME(mod9linear, "theta") 
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Repr ~ z + SiteID + (1|PlotID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial),
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
		write.csv(finaltable,"ReprModBootCoef_plots.csv", row.names = T) 
		save(res, file = "ReprModBootReplicates500_plots.rda")

	

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part G: Make Response Plots.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Make plots of results 
		# Graph out results to evaluate & check for other possible problems.  
		# conditional partial plots are fine, but run into the problem with MEMs of 
		# group level vs. average. Therefore we really want the average marginal probability.
		# will get the average marginal probability by plotting the average of each group - 
		# so will have to calculate the conditional probabilities for each group and then average them.
		summary(Repr.plot.data$z); summary(Repr.plot.data$z1) # values range from 0 - 10. 
		jvalues <- with(Repr.plot.data, seq(from = min(z), to = max(z), length.out = 100)) # sequence of values to plot.
		tmpdat <- Repr.plot.data # copy the object over as back up, predictions will be added to this data frame.  
		# calculate predicted probabilities and store in a list
			pp <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				predict(mod9linear, newdata = tmpdat, type = "response")
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
		pdf(file="02.4.4_Flower_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_line() + ylim(c(0, 1))	
		# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))
		# refit model with year & Year & re-plot
		# convert Years down to just three levels 1011,1112,1213 and make year a facot variable 
		library(plyr)
		# calculate predicted probabilities and store in a list

#############################################
# RESPONSE PLOTS FOR EACH YEAR 		
	# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(Repr.plot.data$Year), function(stage) {
					  tmpdat$Year[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod9linear, newdata = tmpdat, type = "response")
					  })
					})
		# get means and quartiles for all jvalues for each level of Year
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
		# add Year
		plotdat2$Year <- factor(rep(levels(Repr.plot.data$Year), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for N, C, S
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Year), alpha = .15) +
					  geom_line(aes(colour = Year), size = 2) +
					  scale_fill_manual(values=c("green", "blue", "red")) +
					scale_colour_manual(values=c("green", "blue", "red"))+
					  ylim(c(0, 1))

##############################################
# do the same thing for REGION NOW 
# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(Repr.plot.data$Region), function(stage) {
					  tmpdat$Region[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod9linear, newdata = tmpdat, type = "response")
					  })
					})
		# get means and quartiles for all jvalues for each level of Region
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
		# add Region
		plotdat2$Region <- factor(rep(levels(Repr.plot.data$Region), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for N, C, S
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  scale_fill_manual(values=c("blue", "red")) +
					scale_colour_manual(values=c("blue", "red"))+
					  ylim(c(0, 1))
		
	
	#############################################################################	
	
##############################################
# do the same thing for SITE NOW 
	
# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(Repr.plot.data$SiteID), function(stage) {
					  tmpdat$SiteID[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod9linear, newdata = tmpdat, type = "response")
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
		# add SiteID
		plotdat2$SiteID <- factor(rep(levels(Repr.plot.data$SiteID), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for each site
			library(RColorBrewer)
			#display.brewer.all()
			n = length(levels(plotdat2$SiteID))
			brewer.pal(n, "Set1")	
		
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SiteID), alpha = .15) +
					  geom_line(aes(colour = SiteID), size = 2) +
					  scale_fill_manual(values=brewer.pal(n, "Set1")) +
					scale_colour_manual(values=brewer.pal(n, "Set1"))+
					  ylim(c(0, 1))
		
		dev.off()
	
	#############################################################################	
		





	
