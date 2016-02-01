###########################################################
# FECUNDITY FUNCTION FOR 2010 - 2012 DEMOGRAPHY DATA
	# This script is run remotely from "04.1_basicIPMdemography.R"

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
# add on Fruit data from 2012 
setwd(path.dat)
addOn2012 <- read.csv(file="demo_addonFecRepr2012Data.csv")
# log transformation of basic size data 
addOn2012$z <- log(addOn2012$z + 1)
# convert regions down to just three levels N, C & S and make year a facot variable 
addOn2012$Region <- revalue(addOn2012$Region, c("C1"="C", "C2"="C", "C3"="C", "N1"="N", "N2"="N", "S1"="S", "S2"="S"))
addOn2012$Year <- as.factor(addOn2012$Year)
# also want to round of fruit counts for plants to whole numbers so that they work with poisson distribution 
addOn2012$Fec <- round(addOn2012$Fec, digits = 0)
addOn2012$Fec <- as.integer(addOn2012$Fec)
colnames(addOn2012); dim(addOn2012)
colnames(dfclean); dim(dfclean)
dim(dfclean)[1] + dim(addOn2012)[1] #
##################
# bind dataframes together to estimate vital rates coefficients for reproductive functions. 
dfclean <- rbind(dfclean, addOn2012)	
dim(dfclean) # should be: 9556
	


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
			pdf(file="04.2.4_Fec_Explor.pdf", width=11, height=8.5)
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
			source("demo_04.2.4.extraplotsFec.R")
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
				#A - (1|SiteID)	random SiteID intercept
				#B - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
				#C - (z|SiteID) = (1+z|SiteID)	random slope of z within SiteID with correlated intercept for site
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

	# Run glmer with all potential fixed effect terms (no three way interaction)
		testA_plot <- glmer(Fec ~ z + Year + Region + z*Region + (1|PlotID), na.action=na.omit, data=training, family = poisson)
		testA_site <- glmer(Fec ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=training, family = poisson)
		testA_both <- glmer(Fec ~ z + Year + Region + z*Region + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=training, family = poisson)
		testA_neither <- glm(Fec ~ z + Year + Region + z*Region, na.action=na.omit, data=training, family = poisson)

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
	pdf(file="04.2.4_Fecundity Random Effect Structure.pdf", width=11, height=8.5)
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
	
	# Having plot alone seems to work better than having site alone, but results seems counter 
	# -intuitive to experimental design & field observations. 
	# Will have to reassess and see how sensitivity final results are to estimates from each random 
	# effect structure. Site should remain as a random effect because that was more
	# geared to how the study was designed. Intra-plot variability should be high 
	# because of how transects stretched over knolls, depressions, banks and streams. 
	# also replication within plots is occasionally very low (some only have few plant observations). 
	
	testA_plot <- glmer(Fec ~ z + Year + Region + z*Region + (1|PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_site <- glmer(Fec ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_both <- glmer(Fec ~ z + Year + Region + z*Region + (1|PlotID) + (1|SiteID) + (1|SiteID:PlotID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	testA_neither <- glm(Fec ~ z + Year + Region + z*Region, na.action=na.omit, data=Fec.plot.data, family = poisson)
	
	
	
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
		mod2 <- glm(Fec ~ z + Year + Region + z*Year + z*Region, data=Fec.plot.data, family = poisson)
		mod3 <- glm(Fec ~ z + Year + Region + z*Year, data=Fec.plot.data, family = poisson)
		mod4 <- glm(Fec ~ z + Year + Region + z*Region, data=Fec.plot.data, family = poisson)
		mod5 <- glm(Fec ~ z + Year + Region, data=Fec.plot.data, family = poisson)
		mod6 <- glm(Fec ~ z + Region + z*Region, data=Fec.plot.data, family = poisson)
		mod7 <- glm(Fec ~ z + Year + z*Year, data=Fec.plot.data, family = poisson)
		mod8 <- glm(Fec ~ z + Region, data=Fec.plot.data, family = poisson)
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
		FullQuasi <- glm(Fec ~ z + I(z^2) + Year + Region + z*Year + z*Region, data=Fec.plot.data, family = quasipoisson)
			drop1(FullQuasi,test = "F")
				FullQuasi <- glm(Fec ~ z + I(z^2) + z*Year + z*Region, data=Fec.plot.data, family = quasipoisson)
					drop1(FullQuasi,test = "F")
						FullQuasi <- glm(Fec ~ I(z^2) + z*Year + z*Region, data=Fec.plot.data, family = quasipoisson)
							drop1(FullQuasi,test = "F")
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
		pdf(file="04.2.4_Fecundity_overdispersion.pdf", width=11, height=8.5)
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
		NegBinom <- glm.nb(Fec ~ z + I(z^2) + Year + Region + z*Year + z*Region, data=Fec.plot.data, link = "log")
		summary(NegBinom, cor = FALSE) # watch=out cor must equal false for correct output of coefficents 
		# find optimal models or set of top models
		anova(NegBinom, test = "Chi")
		drop1(NegBinom, test = "Chi") # drop z:Year
			NegBinom <- glm.nb(Fec ~ z + I(z^2) + Year + Region + z*Region, data=Fec.plot.data, link = "log")
				drop1(NegBinom, test = "Chi") # drop z:region
				# test with anova 
				NegBinomSimp <- glm.nb(Fec ~ z + Year + Region, data=Fec.plot.data, link = "log")
				anova(NegBinom, NegBinomSimp)
		stepAIC(NegBinom)
		# from stepwise selection top = Fec ~ z + I(z^2) + Year + Region + z * Region (seems reasonable & probably not super-overfitting data)
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
	# mod2 <- glmer(Fec ~ z + Year + Region + z*Year + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod3 <- glmer(Fec ~ z + Year + Region + z*Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod4 <- glmer(Fec ~ z + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod5 <- glmer(Fec ~ z + Year + Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod6 <- glmer(Fec ~ z + Region + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod7 <- glmer(Fec ~ z + Year + z*Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod8 <- glmer(Fec ~ z + Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod9 <- glmer(Fec ~ z + Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod10 <- glmer(Fec ~ z + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod13 <- glmer(Fec ~ z + I(z^2) + Year + Region + z*Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod14 <- glmer(Fec ~ z + I(z^2) + Year + Region + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod15 <- glmer(Fec ~ z + I(z^2) + Year + Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod16 <- glmer(Fec ~ z + I(z^2) + Region + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod17 <- glmer(Fec ~ z + I(z^2) + Year + z*Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod18 <- glmer(Fec ~ z + I(z^2) + Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod19 <- glmer(Fec ~ z + I(z^2) + Year + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
	mod20 <- glmer(Fec ~ z + I(z^2) + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)

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
		# model 2 seems to be the 'best' based on Akaike-weights alone, but it is highly suspect. 
		# mod2 has two interactive terms in it - too complex. This is because the data is probably tightly
		# tightly fitting a limited dataset. When we run mod2 through the bootstrapping procedure below
		# we find that confidence intervals for both the two way interactive terms become  
		# non-signifigant and overlap zero. 
		
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
			summary(mod13)
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
		overdisp_fun(mod14) # looking at the ratio & still about 41... 
		overdisp_fun(mod4) # looking at the ratio & still about 30... 
		overdisp_fun(mod13) # looking at the ratio & still about 30... 
		#overdisp_fun(mod2)
		# major with models above... clearly zero-inflated & extremely high over dispersion
		# will attempt to resolve below  
		
	
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part D: ZIPGLMM - 'Fit a Zero Inflated General Linear Mixed Model'
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
# (ZIP) CHECK FOR ZERO INFLATION & INCORPORATE ZI INTO POISSION (ZIP)
	# (ZIP, Zero Inflated Poisson + GLMM, General linear mixed effect model)
	# Owls example: a zero-inflated, generalized linear mixed model for count data
	# Authors: Ben Bolker, Mollie Brooks, Beth Gardner, Cleridy Lennert, Mihoko Minami (Oct 23, 2012)
	#https://groups.nceas.ucsb.edu/non-linear-modeling/projects/owls/WRITEUP/owls.pdf
	# Zero-inflated generalized linear mixed-effects models (ZIGLMMs) 
		hist(Fec.plot.data$Fec, breaks=500, xlim=c(0,500))
		# problem is likely that most flowering plants didn't produce fruits 
		# also a few produced thousands 
		hist(Fec.plot.data$Fec)
		head(sort(-Fec.plot.data$Fec))
		# since we are modelling size on the log scale its no wonder we were having 
		# so many problems with the models above. 
		table(Fec.plot.data$Fec == 0) # about 1/5th of the plants that flowered, didn't produce any fruits
		ppois(0, mean(Fec.plot.data$Fec))
	# Setup priors, formula ect. 
		library(MCMCglmm)
		offvec <- c(1,1,2,rep(1,5)) # set up an offset variable for later. Will be in position 3 of the parameter vector 
		# See comments in section 3.1 of Bolker's tutorial, awkward coding to specify a fixed effect in MCMCglmm 
		#"trait" is a component of the MCMCglmm formula, not a variable in out dataset, so don't worry		
		# there is a zeroinflated binary aspect to the model (binary 0/1, were any fruits produced. There is also a count aspect to the model (1,2,3,4... to n fruits produced by a plant) 
		fixef2 <- Fec ~ trait-1+ ## intercept terms for both count and binary compoenents 
			at.level(trait,1):((Year+Region)*z)
		# make a base version of the prior. 
		prior_overdisp <- list(R=list(V=diag(c(1,1)),nu=0.002,fix=2),
								G=list(list(V=diag(c(1,1e-6)),nu=0.002,fix=2)))
	# Fit the ZIPGLMM
		mt1 <- system.time(mfit1 <- MCMCglmm(fixef2, # formula from above 
											rcov = ~idh(trait):units, 
											random= ~idh(trait):SiteID, # nested structure, similar to random intercept of SiteID 
											prior= prior_overdisp, # setup from Bolker
											data= Fec.plot.data, 
											family="zipoisson",
											verbose=FALSE))
		# Bolker fxn to abbreviate variable names, not even sure if necessary 
		abbfun <- function(x) { 
			gsub("(Sol\\.)*(trait|at.level\\(trait, 1\\):)*", "",
			gsub("Year", "Year", x))
			}
		colnames(mfit1$Sol) <- abbfun(colnames(mfit1$Sol))
		# check results of model 
		mfit1$Residual$nrt<-2 # weird error, resolved from online forum http://comments.gmane.org/gmane.comp.lang.r.lme4.devel/11706
		summary(mfit1)
	# Visualize results
		install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")
		library(coefplot2)
		op <- par(mfrow=c(2,1))
		vn1 <- abbfun(rownames(coeftab(mfit1)))
		coefplot2(mfit1, intercept=TRUE, varnames=vn1) # basic plot (warning z is not scaled)
		coefplot2(mfit1, var.idx=c(1,3), ptype="vcov", main="")
		# support for interactive term of size*Year or size*Region is weak,
		# but marginal results are interesting and do fit with our field observations 
	# Trace Plots - Now, look at trace plots. We are hopping that trace plots look like messy 
		# noise with no obvious directional trends. 
		print(xyplot(mfit1$Sol, layout=c(3,3)))
		#zi_Fec is questionable, effective size is also rather lower than we would have hoped. 
		round(sort(effectiveSize(mfit1$Sol)))
		# also can check for convergence by 'gewke.diag' which gives similarity between first 10% & last 50% of the chain.
		geweke.diag(mfit1$Sol) # based on Z-statistic (is |x| < 1.96?)
		# again, interaction between region & size doesn't have support. 
	# Explore variance parameters 
		vv <- mfit1$VCV 
		# drop uninformative ZI random effects
		vv <- vv[,c("Fec.SiteID", "Fec.units")]
		print(xyplot(vv, layout=c(1,2)))
		effectiveSize(vv)
		print(densityplot(mfit1$Sol, layout=c(3,3)))
		print(densityplot(vv, layout=c(2,1)))
	# Compare with other models 
		basicPoisson <- glm(Fec ~ z + Year + Region + z*Year + z*Region, data=Fec.plot.data, family = poisson)
		basicPoisson_mixed <- glmer(Fec ~ z + Year + Region + z*Year + z*Region + (1|SiteID), na.action=na.omit, data=Fec.plot.data, family = poisson)
		coefplot2(list("glm"=basicPoisson, "glmer"=basicPoisson_mixed, "MCMCglmm"=mfit1),
			intercept=TRUE, legend=TRUE, legend.x="right")
		
	# try with negative binomial 
		library(glmmADMB)
		gt2 <- system.time(gfit2 <- glmmadmb(Fec ~z + Year + Region + z*Year + z*Region+
				(1|SiteID),
				data=Fec.plot.data,
				zeroInflation=TRUE,
				family="nbinom"))
				
		cg2tab <- coeftab(gfit2)[2:10,]
		cm1tab <- coeftab(mfit1)[3:10,]
				
		coefplot2(list(MCMCglmm=cm1tab,glmmADMB_NB=cg2tab),merge.names=FALSE,intercept=TRUE,
			legend=TRUE)
			
		AIC(basicPoisson,basicPoisson_mixed,gfit2)	
		
		
	dev.off()	
######################################################################################	
	# Find top model & predict data 
	# load in functions from Merow et al 2014 doi: 10.1111/ecog.00839
	setwd(path.funct); dir()
	source('MerowFunctions2014IPMSDM.R') # stepDIC, predict.MCMCglmm.cm, shrink.matrix
	# find optimal model with DIC values 
	mod2 <- Fec ~ trait-1 + at.level(trait,1):((Year+Region)*z)
	mod3 <- Fec ~ trait-1 +  at.level(trait,1):(Region) + at.level(trait,1):((Year)*z)
	mod4 <- Fec ~ trait-1 +  at.level(trait,1):(Year) + at.level(trait,1):((Region)*z)
	mod5 <- Fec ~ trait-1 +  at.level(trait,1):(Year) + at.level(trait,1):(z) + at.level(trait,1):(Region)
	mod6 <- Fec ~ trait-1 +  at.level(trait,1):((Region)*z)
	mod7 <- Fec ~ trait-1 +  at.level(trait,1):((Year)*z)
	mod8 <- Fec ~ trait-1 +  at.level(trait,1):(z) + at.level(trait,1):(Region)
	mod9 <- Fec ~ trait-1 +  at.level(trait,1):(z) + at.level(trait,1):(Year)
	mod10 <- Fec ~ trait-1 +  at.level(trait,1):(z)
	modList <- c(mod2, mod3,mod4,mod5,mod6, mod7, mod8, mod9, mod10)
	
	frammy <- data.frame()
	for(i in 1:length(modList)){
		current <- modList[i]
		mt1 <- system.time(mfit1 <- MCMCglmm(current[[1]], rcov = ~idh(trait):units, random= ~idh(trait):SiteID, prior= prior_overdisp, data= Fec.plot.data,  family="zipoisson", verbose=FALSE))
		myDIC <- DIC(mfit1)
		add_this <- data.frame(mod=paste0("mod",i+1), DIC=round(myDIC, 2))
		frammy <- rbind(frammy, add_this)
		}
	(FecMods <- frammy[order(frammy$DIC), ]) # not sure if akaki weights apply here? maybe just from coef plot 
	# mod 3 wins! same with glmer and glm, same structure as above!!!
	# redo model 

	mt1 <- system.time(mfit1 <- MCMCglmm(mod3, rcov = ~idh(trait):units, random= ~idh(trait):SiteID, prior= prior_overdisp, data= Fec.plot.data,  family="zipoisson", verbose=FALSE))
	# fix glitch
	mfit1$Residual$nrt<-2 #
	coeftab(mfit1); vn1 <- abbfun(rownames(coeftab(mfit1)))
	coefplot2(mfit1, intercept=TRUE, varnames=vn1) # basic plot (warning z is not scaled)
	Feccoefficents <- coeftab(mfit1)

	# save top model 
	setwd(path.obj)
	save(mfit1, file='FecMod.rda')

	# Now need to break down the dataset into the same bootstrapped datasets used in the other
	# scripts & run this same model for each of these to get coefficent estimates. 

	
	
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
		# skip the next two lines because we want to use the same bootstrapped datasets for each vital rate so later we can get CI for lambda estimates  
		# tmp <- sampler(Fec.plot.data, "SiteID", reps = 500) # run at 500 for now, but need to bump this up to 1000 or 10,000 once things get all set up
		# bigdata <- cbind(tmp, Fec.plot.data[tmp$RowID, ])
		setwd(path.obj)
		bigdata <- read.csv(file="BootReplicateDataSets.csv")
		bigdata$Year <- as.factor(bigdata$Year)
	# Steps: refit model to sample data
		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
		# 3. Export data & lme4 on each cluster.
		# 4. Write a fxn to fit the model and return estimates. 
		# * the call MCMCglmm() is wrapped in try to avoid stopping the process if models do not converge. 
		
		# Loop through bootstrap estimates 
			mod3 <- Fec ~ trait-1 +  at.level(trait,1):(Region) + at.level(trait,1):((Year)*z) # formula to work with 
		# Make empty dataframe to inset estimates into 
			Names <- rownames(coeftab(mfit1)[1])
			Replicates <- unique(bigdata$Replicate)
			C = matrix("NA", nrow=length(Replicates), ncol=length(Names)); 
	
		# Now run big loop through replicates 
		# WARNING ~ might take like 6 hrs, leave it to run overnight
			for(i in 1:length(Replicates)){
				ReplicateOfbigdata <- subset(bigdata, subset = Replicate == i)
				object <- try(mooty <- MCMCglmm(mod3, rcov = ~idh(trait):units, random= ~idh(trait):SiteID, 
					prior= prior_overdisp, data = ReplicateOfbigdata,  family="zipoisson", verbose=FALSE))	
				C[i,] <-  coeftab(object)[[1]]
			}

		# adjust column names
		C <- data.frame(C) 
		colnames(C) <- Names
		Cdf <- data.frame(C)
		# save boot replicate coefficents
		setwd(path.obj)
		save(Cdf, file = "FecModBootReplicates500.csv")

		
	

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
				count_linear = as.numeric(Feccoefficents["Sol.traitFec",][1]) + as.numeric(Feccoefficents["Sol.at.level(trait, 1):z",][1])*tmpdat$z # linear prediction for the count part of the model 
				zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][1])	# linear prediction for the zero inflated part of the model 
				zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
				count_number = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			
			})
			
		LL <- lapply(jvalues, function(j) {
				tmpdat$z <- j
			count_linear = as.numeric(Feccoefficents["Sol.traitFec",][2]) + as.numeric(Feccoefficents["Sol.at.level(trait, 1):z",][2])*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][2])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberL = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			})
		
		UL <- lapply(jvalues, function(j) {
				tmpdat$z <- j			
			count_linear = as.numeric(Feccoefficents["Sol.traitFec",][5]) + as.numeric(Feccoefficents["Sol.at.level(trait, 1):z",][5])*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][5])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberU = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
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
		head(plotdat)
		setwd(path.fig)
		pdf(file="04.2.4_Fecundity_ResponsePlots.pdf", width=11, height=8.5)
		# plot average marginal predicted probabilities
		ggplot(plotdat, aes(x = z, y = PredictedFecundity)) + geom_line()	
		# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
		ggplot(plotdat, aes(x = z, y = PredictedFecundity)) + geom_linerange(aes(ymin = Lower,
		ymax = Upper)) + geom_line(size = 2)

#############################################################################
#############################################################################		
#############################################################################		
# refit model with year & region & re-plot
# extremely inefficient coding here sorry 

		# calculate predicted probabilities and store in a list
		ppn <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				count_linear = -3.905168 + 1.027125*tmpdat$z # linear prediction for the count part of the model 
				zi_linear = -3.074622	# linear prediction for the zero inflated part of the model 
				zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
				count_number = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			
			})
			
		LLn <- lapply(jvalues, function(j) {
				tmpdat$z <- j
			count_linear = -5.272003 + 0.963153*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][2])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberL = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			})
		
		ULn <- lapply(jvalues, function(j) {
				tmpdat$z <- j			
			count_linear = -2.4871 + as.numeric(Feccoefficents["Sol.at.level(trait, 1):z",][5])*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][5])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberU = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			})	
			
	
		# Now we can display all the predicted probabilities & plot them 
		# average marginal predicted probability across a few different sizes z 
		round((sapply(ppn[c(0.5, 1, 2, 4, 6, 8, 9, 10)], mean)*100), 1) # in percent
		# get the means with lower and upper quartiles
			plotdatN <- t(sapply(ppn, function(x) {
				c(M = mean(x), quantile(x, c(0.25, 0.75)))
			}))
			plotdatN[,2] <- t(sapply(LLn, function(x) {
				LLn = mean(x)
			}))
			plotdatN[,3] <- t(sapply(ULn, function(x) {
				ULn = mean(x)
			}))
			# add in z values and convert to data frame
			plotdatN <- as.data.frame(cbind(plotdatN, jvalues))
			# better names and show the first few rows
			colnames(plotdatN) <- c("PredictedFecundity", "Lower", "Upper", "z")		

################################################################################			
	# calculate predicted probabilities and store in a list
		pps <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				count_linear = -3.60296 + 1.027125*tmpdat$z # linear prediction for the count part of the model 
				zi_linear = -3.074622	# linear prediction for the zero inflated part of the model 
				zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
				count_number = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			
			})
			
		LLs <- lapply(jvalues, function(j) {
				tmpdat$z <- j
			count_linear = -4.861641 + 0.963153*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][2])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberL = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			})
		
		ULs <- lapply(jvalues, function(j) {
				tmpdat$z <- j			
			count_linear = -2.3327 + as.numeric(Feccoefficents["Sol.at.level(trait, 1):z",][5])*tmpdat$z # linear prediction for the count part of the model 
			zi_linear = as.numeric(Feccoefficents["Sol.traitzi_Fec",][5])	# linear prediction for the zero inflated part of the model 
			zi_prob = exp(zi_linear)/(1+exp(zi_linear)) # logistic transformation for zero inflated part 
			count_numberU = exp(count_linear)*(1-zi_prob)# final transformation to incorporate count & zi part of model 	
			})	
			
	
		# Now we can display all the predicted probabilities & plot them 
		# average marginal predicted probability across a few different sizes z 
		round((sapply(pps[c(0.5, 1, 2, 4, 6, 8, 9, 10)], mean)*100), 1) # in percent
		# get the means with lower and upper quartiles
			plotdatS <- t(sapply(pps, function(x) {
				c(M = mean(x), quantile(x, c(0.25, 0.75)))
			}))
			plotdatS[,2] <- t(sapply(LLs, function(x) {
				LLs = mean(x)
			}))
			plotdatS[,3] <- t(sapply(ULs, function(x) {
				ULs = mean(x)
			}))
			# add in z values and convert to data frame
			plotdatS <- as.data.frame(cbind(plotdatS, jvalues))
			# better names and show the first few rows
			colnames(plotdatS) <- c("PredictedFecundity", "Lower", "Upper", "z")	
		
		plotdat$Region <- "South"
		plotdatS$Region <- "Center"
		plotdatN$Region <- "North"
		plotdat2 <- rbind(plotdatS, plotdatN, plotdat) 

		# show first few rows
		head(plotdat2)	
		# set color scale for N, C, S
		ggplot(plotdat2, aes(x = z, y = PredictedFecundity)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  scale_fill_manual(values=c("green", "blue", "red")) +
					scale_colour_manual(values=c("green", "blue", "red"))

		# It looks like this distribution is skewed a bit across each of the levels
		# can view distribution of predicted probabilities here:

		dev.off()
	
	#############################################################################	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
			
		


	



	





	
