###########################################################
# SURVIVAL FUNCTION FOR 2010 - 2012 DEMOGRAPHY DATA
	# This script is run remotely from "04.1_basicIPMdemography.R"

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table on content:
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# A. Exploratory plots across levels & factors from basic GLM
			# Are there any major outliers ect. Was the log transformation of size sufficent? 
			
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


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part A: Exploratory plots across levels & factors from basic GLM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Are there any major outliers ect. Was the log transformation of size sufficient? 
	# Will run exploratory plots for each site, year & region combination and examine visually. 
	# Combine plots into a large pdf. A bit tedious, but informative & useful before proceeding 
	# with the remainder of the analysis. 
	#############################################################################	
			dev.off(); setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
			pdf(file="04.2.1_Survival_Explor.pdf", width=11, height=8.5)
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
			source("demo_04.2.1.extraplots.R")
			# check for figure in figures folder - looks OK 
			

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part B: Determine appropriate random effect structure
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Many random effect structures are possible for a mixed effect model 
	# Especially when there are both plot and site level nesting. 
	# It might even be possible to have a random slope component in addition 
	# to a random intercept. 
	# Unfortunately defend-able statistical methods to choose an optimal 
	# random effect structure for glm objects are still under development. 
	# discussion: http://www.researchgate.net/post/How_can_I_optimize_the_random_effect_structure_in_a_GLMM
	# good discussion of issue here too: http://glmm.wikidot.com/faq
	# REML procedure not possible here for glm objects. Should rely on own 
	# intuition, experimental setup and a solid hypothesis for each alternative 
	# random effect structure. 
	
		# Possible Alternatives (listed from most believable to least believable): 
				#A - (1|SiteID)	random SiteID intercept
				#B - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
				#C - (z|SiteID) = (1+z|SiteID)	random slope of z within SiteID with correlated intercept for site
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
			# Function to split data into training & testing set (50% training 50% testing)
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
		# and testing sets add one observation from each plot to each dataset
			df2 <- lapply(split(dfclean, dfclean$PlotID), function(subdf) subdf[sample(1:nrow(subdf), 1, replace=TRUE),])		
			df2 <- do.call('rbind', df2)
			head(df2, 3); dim(df2)
			training <- rbind(training, df2)
			testing <- rbind(testing, df2)

	# Run glmer with all potential fixed effect terms (no three way interaction)
		testA_plot <- glmer(Surv ~ z + Year + Region + z*Region + (1|PlotID), na.action=na.omit, data=training, family = binomial)
		testA_site <- glmer(Surv ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=training, family = binomial)
		testA_both <- glmer(Surv ~ z + Year + Region + z*Region + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=training, family = binomial)
		testA_neither <- glm(Surv ~ z + Year + Region + z*Region, na.action=na.omit, data=training, family = binomial)

		
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
	
	toAdd <- data.frame(holdout, N, PlotA, SiteA, BothA, NeithA)
	
	myFrame	<- rbind(myFrame, toAdd)
	}
	
	# save final evaluation figure
	setwd(path.fig)
	pdf(file="04.2.1_Sur Random Effect Structure.pdf", width=11, height=8.5)
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
		par(mfrow=c(1,3))
		sjp.glmer(testA_plot, type = "fe")
		sjp.glmer(testA_site, type = "fe")
		sjp.glmer(testA_both, type = "fe")
	
	dev.off()
	
	# Having plot alone seems to work better than having site alone, but results seems counter 
	# intuitive to experimental design & field observations. 
	# Will have to reassess and see how sensitivity final results are to estimates from each random 
	# effect structure. Site should remain as a random effect because that was more
	# geared to how the study was designed. Intra-plot variability should be high 
	# because of how transects stretched over knolls, depressions, banks and streams. 
	# also replication within plots is occasionally very low (some only have few plant observations). 
	


	
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
	
	# Possible models(yes could be coded more efficiently)  
	#mod1 <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod2 <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod3 <- glmer(Surv ~ z + Year + Region + z*Year + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod4 <- glmer(Surv ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod5 <- glmer(Surv ~ z + Year + Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod6 <- glmer(Surv ~ z + Region + z*Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod7 <- glmer(Surv ~ z + Year + z*Year + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod8 <- glmer(Surv ~ z + Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod9 <- glmer(Surv ~ z + Year + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
	mod10 <- glmer(Surv ~ z + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)

	# package for getting effects & AIC values 
	require(MuMIn)
	
	#################
	# MODEL SELECTION 
	#################
		# use the mod.sel function to conduct model selection
		# and put output into object out.put
		out_put <- mod.sel(mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10)
		modList <- c("mod2","mod3","mod4","mod5","mod6","mod7","mod8","mod9","mod10")
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
		write.csv(sel.table,"SurvModAICTable.csv", row.names = T) 
			
	# VARIABLE IMPORTANCE (BASED ON AIC)	
		importance(out_put)
	
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
			overdisp.glmer(mod4) # seems ok 

		# winning model:
			summary(mod6)
		
		library(glmmML) # run sample model with glmmML - do coefficents change?
		mod6_glmmML <-glmmML(formula=Surv ~ z + Region + z * Region,
             cluster=SiteID,data=dfclean,family=binomial)
		summary(mod6_glmmML)
	
	
	
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
		Model_Best <- mod6_glmmML # best model from steps above
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
		dev.off(); setwd(path.fig)# close existing plotting frame, switch directories and open a new file for saving plots  
			pdf(file="04.2.1_Survival_model_adequacy.pdf", width=11, height=8.5)
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
		# very clearly to be linear, if size is ploted on the log scale. 
	
	
	
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
		fixef(mod6)
		coef(mod.surv.glm)
		# Standard error for models, on the logit scale(?)
		se <- sqrt(diag(vcov(mod6))); (tab <- cbind(Est = fixef(mod6), LL = fixef(mod6) - 1.96 * se, UL = fixef(mod6) + 1.96 * se))
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
		set.seed(20)
			#tmp <- sampler(dfclean, "SiteID", reps = 500) # run at 500 for now, but need to bump this up to 1000 or 10,000 once things get all set up
			#bigdata <- cbind(tmp, dfclean[tmp$RowID, ])
		# LOAD SAME REPLICATE DATASET TO ESTIMATE CI'S FOR LAMBDA ESTIMATES  
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets.csv")
		bigdata$Year <- as.factor(bigdata$Year)
		toMerge <- bigdata[,c("Replicate", "ID")]
		toMerge$ID <- as.factor(toMerge$ID)
		toMerge$Replicate <- as.factor(toMerge$Replicate)
		bigdata <- merge(toMerge, dfclean, by="ID", all.x=TRUE)
		# remove NA values 
		bigdata <- subset(bigdata, !is.na(Surv))
		dim(bigdata)
	
	
#bigdata <- bigdata[which(bigdata$Replicate == "1" | bigdata$Replicate == "2" | bigdata$Replicate == "3"), ]
	
	# Steps: refit model to sample data
		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
		# 3. Export data & lme4 on each cluster.
		# 4. Write a fxn to fit the model and return estimates. 
		# * the call glmer() is wrapped in try to avoid stopping the process if models do not converge. 
		f <- fixef(mod6) # store original parameters 
		r <- getME(mod6, "theta") 
		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
				myboot <- function(i) {
					object <- try(glmer(Surv ~ z + Region + z*Region + (1|SiteID), na.action=na.omit, data = bigdata, subset = Replicate == i, family = binomial,
					start = list(fixef = f, theta = r)), silent = TRUE)
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
		save(mod6, file = "SurvMod6.rda")
		save(res, file = "SurvModBootReplicates500.rda")




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
				predict(mod6, newdata = tmpdat, type = "response")
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
		pdf(file="04.2.1_Survival_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_line() + ylim(c(0, 1))	
		# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))
		# refit model with year & region & re-plot
		# convert regions down to just three levels N, C & S and make year a facot variable 
		library(plyr)
		# calculate predicted probabilities and store in a list
		# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(dfclean$Region), function(stage) {
					  tmpdat$Region[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(mod6, newdata = tmpdat, type = "response")
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
		plotdat2$Region <- factor(rep(levels(dfclean$Region), each = length(jvalues)))
		# show first few rows
		head(plotdat2)	
		# set color scale for N, C, S
		ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  scale_fill_manual(values=c("green", "blue", "red")) +
					scale_colour_manual(values=c("green", "blue", "red"))+
					  ylim(c(0, 1))

		# It looks like this distribution is skewed a bit across each of the levels
		# can view distribution of predicted probabilities here:
		# center
		ggplot(data.frame(Probs = biprobs[[1]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75)) +
		ggtitle("Center - distribution of predicted values")
		# north
		ggplot(data.frame(Probs = biprobs[[2]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75)) + 
		ggtitle("North - distribution of predicted values")
		# south
		ggplot(data.frame(Probs = biprobs[[3]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75)) + 				
		ggtitle("South - distribution of predicted values")
		dev.off()
	
	#############################################################################	
	
	





	
