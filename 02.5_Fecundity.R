###########################################################
# FECUNDITY FUNCTIONs 2013-2015 for transplant data 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table on content:
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# A. Exploratory plots across levels & factors from basic GLM/LM
			# Are there any major outliers ect. Was the log transformation of size sufficient? 
			
	# B. Determine appropriate random effect structure (single random effect of sites vs. plots within sites ect.) 
			# for glmer with PlotID & SiteID as random effects. 
			
	# C. Determine appropriate fixed structure (region, year, size, sites, interactions?)
			# for fecundity glmer with a Akaike-weights framework. 
			
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
	dfclean$z <- log(dfclean$z); summary(dfclean$z)
	dfclean$z1 <- log(dfclean$z1); summary(dfclean$z1)
	# round down fecundity data to actual fruit counts:
	dfclean$Fec <- round(dfclean$Fec, digits = 0) # round 
	dfclean$Fec <- as.integer(dfclean$Fec) # make integer for Poisson regression later.
	plot.range <- c(0, 9) # range of size data to be plotted (~ 0 - 9 log(total stem length))

	# need to include all flowering data in this script 2014/2015 & 2015/2016
	dim(dfclean)[1] # all plants in dataset
	dfclean <- dfclean[!is.na(dfclean$Fec),]	# exclude NA's from Flr & Fec 
	dim(dfclean)[1] # all plants left with flowering data
	# ok we are now ready to get going
## need to sort data by size if plotting line of predicted values  
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
			pdf(file="02.5.1_Fec_Explor.pdf", width=11, height=8.5)
			Fec.plot.data <- subset(dfclean, !is.na(Fec))
			dim(Fec.plot.data) # much less observations than other vital rates
			z <- Fec.plot.data$z # stand alone object 
				# Plotting frame 
					Fec.plot.data <- within(Fec.plot.data, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecundity model
			mod.Fec.glm <- glm(Fec ~ z, family = poisson, data = Fec.plot.data)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.data)
			Fec.ps # (actual values for plotting vs. predicted) 
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 2, xlim = plot.range, xlab  =expression("Size , "*italic(z)), ylab = "Number of fruit", main="Global model")
			points(z, Fec.plot.data$Fec, col="lightgrey", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = Fec.plot.data, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
	
		# Run additional plots (optional code stored in R Code - functions subfolder)
		# additional plots in another script to save space on this sheet & make for efficent coding 
			setwd(path.funct)
			source("02.5.extraplotsFec.R")
			# check for figure in figures folder - looks OK 
			# some sites certainty have very few observation ~ <10 plants made fruits 
			# definitely can't have plot level intercept random effect
			
			

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part B: Determine appropriate random effect structure
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Many random effect structures are possible for a mixed effect model 
	# Especially when there are both plot and site level nesting. 
	# It might even be possible to have a random slope component in addition 
	# to a random intercept. 
	# Unfortunately defendable statistical methods to choose an optimal 
	# random effect structure for glmm objects are still under development. 
	# discussion: http://www.researchgate.net/post/How_can_I_optimize_the_random_effect_structure_in_a_GLMM
	# good discussion of issue here too: http://glmm.wikidot.com/faq
	# REML procedure not possible here for glm objects. Should rely on own 
	# intuition, experimental setup and a solid hypothesis for each alternative 
	# random effect structure. Some approaches are covered here - using AIC & reduction in deviance 
	# Cannot use cross validation like in the survival & flower model (won't work here)
	
		# Possible Alternatives (listed from most believable to least believable): 
				#A - (1|PlotID)	random PlotID intercept
				#B - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
				#C - (z|PlotID) = (1+z|PlotID)	random slope of z within SiteID with correlated intercept for site
				#D - (z|site/PlotID) = (z|site)+(z|site:PlotID) random slope of z within both site and plot. intercept varying among sites and among PlotIDs within sites (nested random effects)
				#E - site+(1|site:PlotID)	fixed effect of sites plus random variation in intercept among PlotIDs within sites, unlikely & difficult to work with. Too many site levels. 				
			# Not considered: 
				#F - (0+z|SiteID) = (-1+z|SiteID)	random slope of z within SiteID: no variation in intercept
				#G - (1|SiteID) + (0+z|SiteID)	uncorrelated random intercept and random slope within SiteID

		library(lme4)
		library(arm)
		library(sjPlot)
		library(MuMIn)
		
		
		
		holdOutPortions <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.95)
		myFrame <- data.frame()
		myFrame2 <- data.frame()

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
		
		splits <- splitdf(Fec.plot.data, seed=808) # splits dataframe 
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
			df2 <- lapply(split(Fec.plot.data, Fec.plot.data$PlotID), function(subdf) subdf[sample(1:nrow(subdf), 1, replace=TRUE),])		
			df2 <- do.call('rbind', df2)
			head(df2, 3); dim(df2)
			training <- rbind(training, df2)
			testing <- rbind(testing, df2)
			# need to drop NA's out of testing dataset.
			testing <- testing[complete.cases(testing[,c("Fec", "z", "SiteID", "PlotID")]),]
			
	# Run glmer with all potential fixed effect terms (no three way interaction)
		testA_plot <- glmer(Fec ~ z + (1|PlotID), na.action=na.omit, data=training, family = poisson)
		testA_site <- glmer(Fec ~ z + (1|SiteID), na.action=na.omit, data=training, family = poisson)
		testA_both <- glmer(Fec ~ z + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=training, family = poisson)
		testA_neither <- glm(Fec ~ z, na.action=na.omit, data=training, family = poisson)

	# Predict models to testing data
		ModPlotPreds <- predict(testA_plot, newdata = testing, type="response")
		ModSitePreds <- predict(testA_site, newdata = testing, type="response")
		ModBothPreds <- predict(testA_both, newdata = testing, type="response")
		ModNeithPreds <- predict(testA_neither, newdata = testing, type="response")

	testing$ModPlotPreds <- ModPlotPreds
	testing$ModSitePreds <- ModSitePreds
	testing$ModBothPreds <- ModBothPreds
	testing$ModNeithPreds <- ModNeithPreds

	# AIC FRAME		
	PlotA <- round(AIC(testA_plot), 0)
	BothA <- round(AIC(testA_both), 0)
	SiteA <- round(AIC(testA_site), 0)
	NeithA <- round(AIC(testA_neither), 0)
	holdout <- holdOutPortions[i]
	N <- dim(training)[1]*holdOutPortions[i]
	toAdd <- data.frame(holdout, N, PlotA, SiteA, BothA, NeithA)
	myFrame	<- rbind(myFrame, toAdd)
	
	# pseudo-r2 frame 
	PlotA <- as.numeric(r.squaredGLMM(testA_plot)[1])
	BothA <- as.numeric(r.squaredGLMM(testA_both)[1])
	SiteA <- as.numeric(r.squaredGLMM(testA_site)[1])
	NeithA <- as.numeric(r.squaredGLMM(testA_neither)[1])
	toAdd <- data.frame(holdout, N, PlotA, SiteA, BothA, NeithA)
	myFrame2 <- rbind(myFrame2, toAdd)

	}
	
	# save final evaluation figure
	setwd(path.fig)
	pdf(file="02.5.2_Fecundity Random Effect Structure.pdf", width=11, height=8.5)
	colnames(myFrame) <- c("trainPor", "N", "PlotA", "SiteA", "BothA", "NeithA") 
	plot(myFrame$trainPor, myFrame$PlotA, xlab="Portion of data Training:Testing ratio", 
		ylab="AIC score", col=1, ylim=c(min(myFrame[,3:5]) - 0.04, max(myFrame[,3:5]) + 0.04), pch=19, cex=1.5, main="Fecundity Random Effect Structure")
	lines(myFrame$trainPor, myFrame$PlotA, col=1)
	points(myFrame$trainPor, myFrame$SiteA, col=2, pch=19); lines(myFrame$trainPor, myFrame$SiteA, col=2)
	points(myFrame$trainPor, myFrame$BothA, col=3, pch=19); lines(myFrame$trainPor, myFrame$BothA, col=3)
	points(myFrame$trainPor, myFrame$NeithA, col="darkgrey", pch=19); lines(myFrame$trainPor, myFrame$NeithA, col="darkgrey")
	legend('topright', legend = c("(1|PlotID)", "(1|SiteID)", "(1|PlotID) + (1|SiteID) + (1|SiteID:PlotID)", "glm"), col = c(1:3, "darkgrey"), cex = 0.8, pch = 19)
	text(myFrame$trainPor, min(myFrame[,3:5]) - 0.02, labels = round(myFrame$N, 0))
	
	colnames(myFrame2) <- c("trainPor", "N", "PlotA", "SiteA", "BothA", "NeithA") 
	plot(myFrame2$trainPor, myFrame2$PlotA, xlab="Portion of data Training:Testing ratio", 
		ylab="pseudo R-squared", col=1, ylim=c(min(myFrame2[,3:5]) - 0.04, max(myFrame2[,3:5]) + 0.04), pch=19, cex=1.5, main="Fecundity Random Effect Structure")
	lines(myFrame2$trainPor, myFrame2$PlotA, col=1)
	points(myFrame2$trainPor, myFrame2$SiteA, col=2, pch=19); lines(myFrame2$trainPor, myFrame2$SiteA, col=2)
	points(myFrame2$trainPor, myFrame2$BothA, col=3, pch=19); lines(myFrame2$trainPor, myFrame2$BothA, col=3)
	points(myFrame2$trainPor, myFrame2$NeithA, col="darkgrey", pch=19); lines(myFrame2$trainPor, myFrame2$NeithA, col="darkgrey")
	legend('topright', legend = c("(1|PlotID)", "(1|SiteID)", "(1|PlotID) + (1|SiteID) + (1|SiteID:PlotID)", "glm"), col = c(1:3, "darkgrey"), cex = 0.8, pch = 19)
	text(myFrame2$trainPor, min(myFrame2[,3:5]) - 0.02, labels = round(myFrame2$N, 0))
	
	dev.off()
	
	# Having plot alone seems to work better than having site alone.
	# Will have to reassess and see how sensitivity final results are to estimates from each random 
	# effect structure. Plot should remain as a random effect because that was more
	# geared to how the study was designed. Intra-plot variability should be low 
	# because of how plots were established in nearly uniform areas.
	# but unfortunately replication within plots is occasionally very low (some only have few plant observations). 
	
	testA_plot <- glmer(Fec ~ z + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_site <- glmer(Fec ~ z + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_both <- glmer(Fec ~ z + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_neither <- glm(Fec ~ z, na.action=na.omit, data=Fec.plot.data, family = poisson)
	
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part C: Determine appropriate fixed structure (region, year, size, sites, interactions?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# For fecundity glmer with a Akaike-weights framework. 
	# Largest possible full model:
		#"Fec ~ z + Year + Region + z*Year + z*Region + z*Year*Region"
	# Procedure followed here from: 
		# Chapter 9: GLM and GAM for count data & (example in chapter 16)
			# Zuur 2009 - "Mixed effect models and extension in ecology with R"
			# code from: https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	#############################################################################	
	
	# try at first with just a simple glm - no mixed effects & no nested structure
		mod2 <- glm(Fec ~ z + Year + SiteID + z*Year + z*SiteID, data=Fec.plot.data, family = poisson)
		mod3 <- glm(Fec ~ z + Year + SiteID + z*Year, data=Fec.plot.data, family = poisson)
		mod4 <- glm(Fec ~ z + Year + SiteID + z*SiteID, data=Fec.plot.data, family = poisson)
		mod5 <- glm(Fec ~ z + Year + SiteID, data=Fec.plot.data, family = poisson)
		mod6 <- glm(Fec ~ z + SiteID + z*SiteID, data=Fec.plot.data, family = poisson)
		mod7 <- glm(Fec ~ z + Year + z*Year, data=Fec.plot.data, family = poisson)
		mod8 <- glm(Fec ~ z + SiteID, data=Fec.plot.data, family = poisson)
		mod9 <- glm(Fec ~ z + Year, data=Fec.plot.data, family = poisson)
		mod10 <- glm(Fec ~ z, data=Fec.plot.data, family = poisson)
	
	# compare AIC score of models with no random effect structure 
		out_put <- mod.sel(mod3,mod4,mod5,mod6,mod4,mod8,mod9,mod10)
		modList <- c("mod3","mod4","mod5","mod6","mod4","mod8","mod9","mod10")
		out_put # same order of models as below
	
	# check these non-mixed models for overdispersion 
		# - Overdispersion may be real or apparent here (apparent: outliers, nested structure, missing covariates
		# or non-linear effects & real overdispersion: exists when we cannot identify the causes
		# and the variation really is greater than the mean. 
		# - Looks like some major overdipersion issues here. Might want to switch 
		# to a negative binomial distribution.
		# try to add a quadratic term for size. Does this help?
		FullQuasi <- glm(Fec ~ z + I(z^2) + Year + SiteID + z*Year + z*SiteID, data=Fec.plot.data, family = quasipoisson)
			drop1(FullQuasi,test = "F")
				FullQuasi <- glm(Fec ~ z + Year + SiteID + z*SiteID, data=Fec.plot.data, family = quasipoisson)
					drop1(FullQuasi,test = "F")
							#drop1(FullQuasi,test = "F")
		# Hesitant to use this model with only interactions, also since we have seen the random effect structure above 
		# we should use mixed models instead
									
	# Visually check assumptions of top non-mixed effect model (just the glms above)
		# Should plot all variables here & fitted/residuals to check assumptions
		# also would be good to explore visually how overdispersion could be improved. 
		# View three types of residuals 
		EP <- resid(FullQuasi, type = "pearson") # Pearsons 
		ED <- resid(FullQuasi, type = "deviance") # Deviance residuals 
		mu <- predict(FullQuasi, type = "response") # Predicted values from model 
		E <- Fec.plot.data$z - mu # classic residuals 
		EP2 <- E / sqrt(as.numeric(summary(FullQuasi)[14]) * mu)
		setwd(path.fig)
		pdf(file="02.5.3_Fecundity_overdispersion.pdf", width=11, height=8.5)
		par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
		plot(FullQuasi, ask=FALSE)
		op <- par(mfrow = c(2, 2))
		# Response residuals (observed minus fitted values, also called ordinary residuals),
		plot(x = mu, y = E, main = "Response residuals", sub="Response residuals (observed minus fitted values, also called ordinary residuals)")
		# Pearson residuals, scaled Pearson residuals (the overdispersion is taken into account)
		plot(x = mu, y = EP, main = "Pearson residuals")
		plot(x = mu, y = EP2, main = "Pearson residuals scaled", sub="Pearson residuals, scaled Pearson residuals (the overdispersion is taken into account)")
		# the deviance residuals for the optimal quasi-Poisson mod
		plot(x = mu, y = ED, main = "Deviance residuals", sub="the deviance residuals for the optimal quasi-Poisson mod")
		par(op)
		# dev.off(), continue with negative binomial plots 
	
	# still major violation problems, will try negative binomial 
	# Try with negative binomial:
		library(MASS)
		NegBinom <- glm.nb(Fec ~ z + I(z^2) + Year + SiteID + z*Year + z*SiteID, data=Fec.plot.data, link = "log")
		summary(NegBinom, cor = FALSE) # watch=out cor must equal false for correct output of coefficents 
		# find optimal models or set of top models
		anova(NegBinom, test = "Chi")
		drop1(NegBinom, test = "Chi") # drop z:Year
			NegBinom <- glm.nb(Fec ~ z + I(z^2) + Year + SiteID + z*SiteID, data=Fec.plot.data, link = "log")
				drop1(NegBinom, test = "Chi") # drop z:SiteID
				# test with anova 
				NegBinomSimp <- glm.nb(Fec ~ z + SiteID, data=Fec.plot.data, link = "log")
				anova(NegBinom, NegBinomSimp)
		stepAIC(NegBinom)
		# from stepwise selection top = Fec ~ z + SiteID (seems reasonable & probably not super-overfitting data)
		# Since p-values are approximate here anything close to the magic 5% threshold without sound ecol explan
		# should probably be dropped. 
		NegBinom <- glm.nb(Fec ~ I(z^2), data=Fec.plot.data, link = "log")
		op <- par(mfrow = c(2, 2))
		 plot(NegBinom, ask=FALSE) # seems better than poisson
		mtext("Negative binomial", side=3, line=1, outer=TRUE, cex=0.8, font=2)
		par(op)
		#dev.off()
	
######################################################################################	
######################################################################################

# Now compare models with random effect structure glmer, nested site level intercepts. 		
	# Possible competing models(yes could be coded more efficiently)  
	# mod2 <- glmer(Fec ~ z + Year + SiteID + z*Year + z*SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod3 <- glmer(Fec ~ z + Year + SiteID + z*Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod4 <- glmer(Fec ~ z + Year + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod5 <- glmer(Fec ~ z + Year + SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod6 <- glmer(Fec ~ z + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod7 <- glmer(Fec ~ z + Year + z*Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod8 <- glmer(Fec ~ z + SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod9 <- glmer(Fec ~ z + Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod10 <- glmer(Fec ~ z + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod13 <- glmer(Fec ~ z + I(z^2) + Year + SiteID + z*Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod14 <- glmer(Fec ~ z + I(z^2) + Year + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod15 <- glmer(Fec ~ z + I(z^2) + Year + SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod16 <- glmer(Fec ~ z + I(z^2) + SiteID + z*SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod17 <- glmer(Fec ~ z + I(z^2) + Year + z*Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod18 <- glmer(Fec ~ z + I(z^2) + SiteID + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod19 <- glmer(Fec ~ z + I(z^2) + Year + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod20 <- glmer(Fec ~ z + I(z^2) + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)

	# package for getting effects & AIC values 
	require(MuMIn)
	
	#################
	# MODEL SELECTION 
	#################
		# use the mod.sel function to conduct model selection
		# and put output into object out.put
		# First check if use of quadratic term is necessary
		out_put <- mod.sel(mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)
		modList <- c("mod3","mod4","mod5","mod6","mod7","mod8","mod9","mod10","mod13","mod14","mod15","mod16","mod17","mod18","mod19","mod20")
		out_put
		#out_put <- mod.sel(mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10)
		#modList <- c("mod3","mod4","mod5","mod6","mod7","mod8","mod9","mod10")
		#out_put
		# Mod4 & mod14 seem to work  best here 
		# Fec ~ z + Year + SiteID + z*SiteID + (1|SiteID)
		
	# TOP MODELS
		# OPTIONAL - select models 95% cumulative weight criteria 
		# IMPORTANT: Weights have been renormalized!! 
		subset(out_put, cumsum(out_put$weight) <= .95) 
		
# MAKE CLEAN AICc TABLE & SAVE
		# coerce the object out_put into a data frame 
		sel.table <- as.data.frame(out_put)[8:12] 			
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
		# write to a file, here a comma separated values format 
		# make sure your working directory is properly specified 
		setwd(path.obj)
		write.csv(sel.table,"FecModAICTable.csv", row.names = T) 
			
	# VARIABLE IMPORTANCE (BASED ON AIC)	
		importance(out_put)
		# very strong year effect, much greater than region. 
	# winning model:
			summary(mod14)
			# similar choice whether we use glm or glmer, but z:year is switched for z:region 
	
	# Check for overdispersion in final GLMM
		# advice followed from: "How can I deal with overdispersion in GLMMs?"
		# http://glmm.wikidot.com/faq
		# WARNING: This is a very crude approximation of a dispersion parameter. 
			overdisp_fun <- function(model) {
							  ## number of variance parameters in 
							  ##   an n-by-n variance-covariance matrix
							  vpars <- function(m) {
								nrow(m)*(nrow(m)+1)/2
							  }
							  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
							  rdf <- nrow(model.frame(model))-model.df
							  rp <- residuals(model,type="pearson")
							  Pearson.chisq <- sum(rp^2)
							  prat <- Pearson.chisq/rdf
							  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
							  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
							}
		
		# this 'overdisp_fun' is essentially the same as what is available from the 
		# gof & ados packages 
		overdisp_fun(mod14) # looking at the ratio & still about 7... 
		overdisp_fun(mod4) # looking at the ratio & still about 6... 
		overdisp_fun(mod13) # looking at the ratio & still about 7... 
		#overdisp_fun(mod2)
		# major OD with models above... clearly zero-inflated & extremely high over dispersion
		# will attempt to resolve below  
		
	
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part D: ZIPGLMM - 'Fit a Zero Inflated General Linear Mixed Model'
# there will be errors in this section ~ just continue running the script. 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We won't actually use a ZIPGLMM here, since in the tranplant dataset pretty much 
	# all plants that flowered also produced fruits. Over dispersion isn't as much of an 
	# issue here as it was in the demography dataset. A simple Poission model or negative 
	# binomial model will work fine for the fecundity function. 	

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part D2: General poisson 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
# try with negative binomial & alternative types  
		library(glmmADMB)
		library(R2jags)
		library(reshape)
		library(coda)
		library(AICcmodavg)
		gfit1 <- glmmadmb(Fec ~z + z*SiteID + (1|PlotID), data=Fec.plot.data, zeroInflation=FALSE, family="poisson")
		#gfit2_ZI <- glmmadmb(Fec ~z + SiteID + (1|PlotID), data=Fec.plot.data, zeroInflation=TRUE, family="poisson")
		gfit3_NB <- glmmadmb(Fec ~z + z*SiteID + (1|PlotID), data=Fec.plot.data, zeroInflation=FALSE, family="nbinom1")
		#gfit4_NBZI <- glmmadmb(Fec ~z + SiteID + (1|PlotID), data=Fec.plot.data, zeroInflation=TRUE, family="nbino1")
		gfit5_NB2 <- glmmadmb(Fec ~z + z*SiteID + (1|PlotID), data=Fec.plot.data, zeroInflation=FALSE, family="nbinom1")
	
# Plot coefficients across competing models  
		palette("default")
		coefplot2(list(P=gfit1, NB1=gfit3_NB, NB2=gfit5_NB2), 
				varnames=NULL, intercept=TRUE, legend=TRUE)

# What does AIC suggest? 
		AIC(gfit1,gfit3_NB,gfit5_NB2)		
		overdisp_fun(gfit1) 
		overdisp_fun(gfit3_NB) 
		overdisp_fun(gfit5_NB2) 
			
# Look at the variance to mean ratio. Calculate the mean and variance for each group (PlotID)
				library(plyr)
				library(ggplot2)
				mvtab <- ddply(Fec.plot.data,
					.(SiteID:PlotID),
						summarise,
						callmean=mean(Fec),
						callvar=var(Fec))
				# plot the mean vs variance for each plot 
				q1 <- qplot(callmean,callvar,data=mvtab)
				# add on alternative model forms 
				print(q1+
					## linear (quasi-Poisson/NB1) fit
					geom_smooth(method="lm",formula=y~x-1)+
					## smooth (loess)
					geom_smooth(colour="red")+
					## semi-quadratic (NB2/LNP)
					geom_smooth(method="lm",formula=y~I(x^2)+offset(x)-1,colour="purple")+
					## Poisson (v=m)
					geom_abline(a=0,b=1,lty=2))

		# from the graph produced above a straight line implies a quasi-poisson or negative binomial 
		# of the "type 1" group with variance proportional to the mean is likely to be best. 
		# While a variance to mean ratio in the form Variance = mean + (constant)mean^2 implies 
		# a negative binomial TYPE 2 (the standard in R) or a lognormal-Poission fit is likely to 
		# work best. 			
						
		# For are data here the relationship looks fairly linear, but there is just one plot 
		# that seems to act as a major outlier and screw everything up. Its really just Wiley plot 10 
		# which we can forget about here because it was weird. Defiantly not a problem with ZI 
		# because all flowering plants produced at least one fruit. 
		
		# **** What should we do? ******
		# Visually NB1 seems ok (similar to other); AIC suggests NB1;
		# Variance to mean ration suggests NB1
		
# Bootstrap the data to conduct a parametric estimate of the dispersion paarameter (given mixed effects) 
	# rational behind this & code used is found in 
		# "Using observation-level random effects to model overdispersion in count data in ecology and evolution"
		# Harrison XA. (2014) Using observation-level random effects to model overdispersion in count data in ecology and evolution. PeerJ 2:e616 https://dx.doi.org/10.7717/peerj.616


	FecMod_NB <- glmmadmb(Fec ~z + z*SiteID + (1|PlotID),
			data=Fec.plot.data, zeroInflation=FALSE,
			family="nbinom1")
					
	summary(FecMod_NB)
	names(FecMod_NB)
	coef(FecMod_NB)
	setwd(path.obj)
	save(FecMod_NB, file="FecMod.rda")

	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part E: Bootstrap to get approximations for estimates & CI for lambda bootstrap estimates.
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
	
#-------------------------------------------------------------------------------
# SKIP LINES & DO NOT RUN BOOTSTRAPPED REPLICATES ~ NO TIME (84 DAYS)
#-------------------------------------------------------------------------------		
#	# bootstrapping function
#						sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
#						cid <- unique(dat[, clustervar[1]])
#						ncid <- length(cid)
#						recid <- sample(cid, size = ncid * reps, replace = TRUE)
#						if (replace) {
#						rid <- lapply(seq_along(recid), function(i) {
#						cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
#						size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
#						})
#						} else {
#						rid <- lapply(seq_along(recid), function(i) {
#						cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
#						})
#						}
#						dat <- as.data.frame(do.call(rbind, rid))
#						dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
#						labels = FALSE))
#						dat$NewID <- factor(dat$NewID)
#						return(dat)
#						}
#	# Resample data to a new df 100X larger  
#		set.seed(20)
#		# skip the next two lines because we want to use the same bootstrapped datasets for each vital rate so later we can get CI for lambda estimates  
#		#tmp <- sampler(Fec.plot.data, "SiteID", reps = 500) # run at 500 for now, but need to bump this up to 1000 or 10,000 once things get all set up
#		#bigdata <- cbind(tmp, Fec.plot.data[tmp$RowID, ])
#		setwd(path.obj)
#		bigdata <- read.csv(file="BootReplicateDataSets.csv")
#		bigdata$Year <- as.factor(bigdata$Year)
#		bigdata$Replicate <- as.factor(bigdata$Replicate)
#		bigdata <- bigdata[!is.na(bigdata$Fec),]	# exclude NA's from Flr & Fec 
#		
#		# are there enough observations at each site level?
#		table(dfclean$SiteID)	
#		table(bigdata$SiteID)	
#		mySub1 <- bigdata[bigdata$Replicate == 1, ]
#		mySub2 <- bigdata[bigdata$Replicate == 2, ]
#		mySub3 <- bigdata[bigdata$Replicate == 3, ]
#		table(mySub1$SiteID)	
#		table(mySub2$SiteID)	
#		table(mySub3$SiteID)	
#
#
#		# Steps: refit model to sample data
#		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
#		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
#		# 3. Export data & lme4 on each cluster.
#		# 4. Write a fxn to fit the model and return estimates. 
#		# * the call MCMCglmm() is wrapped in try to avoid stopping the process if models do not converge. 
#		
#		# Loop through bootstrap estimates 
#			mod3 <- Fec ~ z + z*SiteID + (1|PlotID) # formula to work with 	
#		
#		# Make empty dataframe to inset estimates into 
#			Names <- rownames(coeftab(FecMod_NB)[1])
#			Replicates <- unique(bigdata$Replicate)
#			C = matrix("NA", nrow=length(Replicates), ncol=length(Names)); 
#			colnames(C) <- Names
#			C <- data.frame(C) 
#			C[,1:length(Names)] <- as.numeric(as.character(C[,1:length(Names)]))
#	
#		# Now run big loop through replicates 
#		# WARNING ~ might take about 6 hrs, leave it to run overnight
#			cl <- makeCluster(4) # start clusters from snow library
#			clusterExport(cl, c("bigdata"))
#			clusterEvalQ(cl, Sys.getenv("HOST"))
#			clusterEvalQ(cl, library(glmmADMB))
#			for(i in 1:length(Replicates)){
#				mySub <- bigdata[bigdata$Replicate == i, ]
#				object <- try(glmmadmb(Fec ~ z + z*SiteID + (1|PlotID), data=mySub, easyFlag=TRUE, zeroInflation=FALSE, family="nbinom1"), silent=TRUE)
#				if(class(object) == "try-error"){
#					next
#				}
#				else {
#					C[i,c("X.Intercept.")] <-  coeftab(object)["(Intercept)",1]
#					C[i,c("z")] <-  coeftab(object)["z",1]
#					C[i,c("SiteIDCOAST")] <-  coeftab(object)["SiteIDCOAST",1]
#					C[i,c("SiteIDHUNTER")] <-  coeftab(object)["SiteIDHUNTER",1]
#					C[i,c("SiteIDLOOK")] <-  coeftab(object)["SiteIDLOOK",1]
#					C[i,c("SiteIDMOSBY")] <-  coeftab(object)["SiteIDMOSBY",1]
#					C[i,c("SiteIDROCK")] <-  coeftab(object)["SiteIDROCK",1]
#					C[i,c("SiteIDTHOMAS")] <-  coeftab(object)["SiteIDTHOMAS",1]
#					C[i,c("SiteIDWILEY")] <-  coeftab(object)["SiteIDWILEY",1]
#					C[i,c("z.SiteIDCOAST")] <-  coeftab(object)["z:SiteIDCOAST",1]
#					C[i,c("z.SiteIDHUNTER")] <-  coeftab(object)["z:SiteIDHUNTER",1]
#					C[i,c("z.SiteIDLOOK")] <-  coeftab(object)["z:SiteIDLOOK",1]
#					C[i,c("z.SiteIDMOSBY")] <-  coeftab(object)["z:SiteIDMOSBY",1]
#					C[i,c("z.SiteIDROCK")] <-  coeftab(object)["z:SiteIDROCK",1]
#					C[i,c("z.SiteIDTHOMAS")] <-  coeftab(object)["z:SiteIDTHOMAS",1]
#					C[i,c("z.SiteIDWILEY")] <-  coeftab(object)["z:SiteIDWILEY",1]
#					rm(object)
#					(print(i)); (Sys.time())
#				}
#			}
#		stopCluster(cl)
#
#		# adjust column names
#		Cdf <- data.frame(C)
#		# save boot replicate coefficents
#		setwd(path.obj)
#		save(Cdf, file = "FecModBootReplicates500.csv")


		
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part F: Save final results for incorporation into IPM.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#		# results are fairly comparable, but possible to see a bit of bias. 
#		# save results to (save final model & save table of boot strapped CI's 
#		setwd(path.obj)
#		write.csv(finaltable,"FecModBootCoef.csv", row.names = T) 
#		
#		# save boot replicate coefficents
#		save(res, file = "FecModBootReplicates500.rda")

		
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part G: Make Response Plots.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Make plots of results 
		# Graph out results to evaluate & check for other possible problems.  
		# conditional partial plots are fine, but run into the problem with MEMs of 
		# group level vs. average. Therefore we really want the average marginal probability.
		# will get the average marginal probability by plotting the average of each group - 
		# so will have to calculate the conditional probabilities for each group and then average them.
		summary(Fec.plot.data$z); summary(Fec.plot.data$z1) # values range from 0 - 10. 
		jvalues <- with(Fec.plot.data, seq(from = min(z), to = max(z), length.out = 100)) # sequence of values to plot.
		tmpdat <- Fec.plot.data # copy the object over as back up, predictions will be added to this data frame.  
		# calculate predicted probabilities and store in a list
		pp <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				# predict(mod4, newdata = tmpdat, type = "response")
				count_number= exp(coef(FecMod_NB)["(Intercept)"]) + exp((coef(FecMod_NB)["z"])*tmpdat$z)
			})
			
		LL <- lapply(jvalues, function(j) {
				tmpdat$z <- j
			count_numberL= exp(confint(FecMod_NB)[c("(Intercept)"), c("2.5 %")]) + exp((confint(FecMod_NB)[c("z"), c("2.5 %")])*tmpdat$z)
			})
		
		UL <- lapply(jvalues, function(j) {
				tmpdat$z <- j			
			count_numberU= exp(confint(FecMod_NB)[c("(Intercept)"), c("97.5 %")]) + exp((confint(FecMod_NB)[c("z"), c("97.5 %")])*tmpdat$z)
			})	
			
	
		# Now we can display all the predicted probabilities & plot them 
		# average marginal predicted probability across a few different sizes z 
		round((sapply(pp[c(0.5, 1, 2, 4, 6, 8, 9, 10)], mean)*100), 1) # in percent
		# get the means with lower and upper quartiles
			plotdat <- t(sapply(pp, function(x) {
				c(M = mean(x), quantile(x, c(0.25, 0.75)))
			}))
			plotdat[,2] <- t(sapply(LL, function(x) {
				LL = mean(x)
			}))
			plotdat[,3] <- t(sapply(UL, function(x) {
				UL = mean(x)
			}))
			
			
		# add in z values and convert to data frame
		plotdat <- as.data.frame(cbind(plotdat, jvalues))
		# better names and show the first few rows
		colnames(plotdat) <- c("PredictedFecundity", "Lower", "Upper", "z")
		plotdat$Site <- "CALAPOOIA"
		head(plotdat)
		setwd(path.fig)
		pdf(file="02.5.4_Fecundity_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedFecundity)) + geom_line()	
		# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedFecundity)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 2)

#############################################################################
#############################################################################		
#############################################################################		
# refit model with SiteID & re-plot
# extremely inefficient coding here sorry 
	# to loop through
	to_loop <- levels(dfclean$SiteID)[2:length(levels(dfclean$SiteID))] #all sites but the one used as a dummy variable. 
	
	for(i in 1:length(to_loop)){
		paste0("SiteID", (to_loop)[i])

		pp2 <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				# predict(mod4, newdata = tmpdat, type = "response")
				count_number= exp((coef(FecMod_NB)["(Intercept)"]) + (coef(FecMod_NB)[paste0("SiteID",(to_loop)[i])]) +
									(coef(FecMod_NB)["z"]*tmpdat$z) + ((coef(FecMod_NB)[paste0("z:SiteID", (to_loop)[i])])*tmpdat$z))					
			})
			assign(paste0("pp_", to_loop[i]), pp2) # assign values for each site 
				
		LL2 <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				# predict(mod4, newdata = tmpdat, type = "response")
				count_numberL= exp((confint(FecMod_NB)["(Intercept)","2.5 %"]) + (confint(FecMod_NB)[paste0("SiteID",(to_loop)[i]),"2.5 %"]) +
					(confint(FecMod_NB)["z","2.5 %"]*tmpdat$z) + ((confint(FecMod_NB)[paste0("z:SiteID", (to_loop)[i]), "2.5 %"])*tmpdat$z))					
			})
			assign(paste0("LL_", to_loop[i]), LL2) # assign values for each site 
			
		UL2 <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				# predict(mod4, newdata = tmpdat, type = "response")
				count_numberU= exp((confint(FecMod_NB)["(Intercept)","97.5 %"]) + (confint(FecMod_NB)[paste0("SiteID",(to_loop)[i]),"97.5 %"]) +
					(confint(FecMod_NB)["z","97.5 %"]*tmpdat$z) + ((confint(FecMod_NB)[paste0("z:SiteID", (to_loop)[i]), "97.5 %"])*tmpdat$z))					
			})
			assign(paste0("UL_", to_loop[i]), UL2) # assign values for each site 	
	
		# Now we can display all the predicted probabilities & plot them 
		# average marginal predicted probability across a few different sizes z 
		round((sapply(pp2[c(0.5, 1, 2, 4, 6, 8, 9, 10)], mean)*100), 1) # in percent
		# get the means with lower and upper quartiles
			plotdatN <- t(sapply(pp2, function(x) {
				c(M = mean(x), quantile(x, c(0.25, 0.75)))
			}))
			plotdatN[,2] <- t(sapply(LL2, function(x) {
				LL2 = mean(x)
			}))
			plotdatN[,3] <- t(sapply(UL2, function(x) {
				UL2 = mean(x)
			}))
			# add in z values and convert to data frame
			plotdatN <- as.data.frame(cbind(plotdatN, jvalues))
			# better names and show the first few rows
			plotdatN$Site <- (to_loop)[i]
			colnames(plotdatN) <- c("PredictedFecundity", "Lower", "Upper", "z", "Site")		
			# bind to existing dataframe 
			plotdat <- rbind(plotdat, plotdatN)
			}
			
		# show first few rows
		head(plotdat)	
		# set color scale for sites
		ggplot(plotdat, aes(x = z, y = PredictedFecundity)) +
					  geom_line(aes(colour = Site), size = 2) +
					  ylim(0, 550)

		# It looks like this distribution is skewed a bit across each of the levels
		# can view distribution of predicted probabilities here:

		dev.off()
	
	#############################################################################	
		
	