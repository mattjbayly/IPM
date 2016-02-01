###########################################################
# SURVIVAL FUNCTION FOR 2010 - 2012 DEMOGRAPHY DATA

## 1.1 - SURVIVAL (binary indicator 'Surv' = 1 if survived)
	dev.off()# close plotting frame 
	z <- dfclean$z
	
surv.plot.data <- within(dfclean, {
    z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
    z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
	z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
	z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
    rm(z.quantiles)
})

# survival model
	mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = dfclean)
	summary(mod.surv.glm)

# calculate groupwise summary statistics (estimates across size range - 16 classes) 
		# (actual values for plotting vs. predicted) 
	surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
	surv.ps

	plot(Surv.mean ~ z.mean,
		 data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1),
		 xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving")

	# Could be fitted line from MEM
	lines(fitted(mod.surv.glm) ~ z, data = dfclean, col = "red")
	add_panel_label(ltype="a") # A.) B.) C.) ... 

	
###########################################################
# Explore with mixed effect model structure 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part 1: glmer exploration of nested random effects (plot nested within site nested within region)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# tutorial from: R Data Analysis Examples: Mixed Effects Logistic Regression
	#http://www.ats.ucla.edu/stat/r/dae/melogit.htm 
	library(lme4)
	require(ggplot2)
	require(GGally)
	require(reshape2)
	require(lme4)
	require(compiler)
	require(parallel)
	require(boot)
	
	# basic model, estimate the model and store results in m
		m_plot <- glmer(Surv ~ z + (1|PlotID), data = dfclean, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
		m_site <- glmer(Surv ~ z + (1|SiteID), data = dfclean, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)	
		m_both <- glmer(Surv ~ z + (1|SiteID/PlotID), family=binomial(link="logit"), na.action=na.omit, data=dfclean)

		# difference across models so far 
			fixef(m_plot)
			fixef(m_site)
			fixef(m_both)
			coef(mod.surv.glm)

		
	# Standard error for models, on the logit scale(?)
		# plot only 
		se <- sqrt(diag(vcov(m_plot))); (tab <- cbind(Est = fixef(m_plot), LL = fixef(m_plot) - 1.96 * se, UL = fixef(m_plot) + 1.96 * se))
		# site only 
		se <- sqrt(diag(vcov(m_site))); (tab <- cbind(Est = fixef(m_site), LL = fixef(m_site) - 1.96 * se, UL = fixef(m_site) + 1.96 * se))
		# both nested 
		se <- sqrt(diag(vcov(m_both))); (tab <- cbind(Est = fixef(m_both), LL = fixef(m_both) - 1.96 * se, UL = fixef(m_both) + 1.96 * se))
		# basic glm
		se <- sqrt(diag(vcov(mod.surv.glm))); (tab <- cbind(Est = coef(mod.surv.glm), LL = coef(mod.surv.glm) - 1.96 * se, UL = coef(mod.surv.glm) + 1.96 * se))
		
	# set up function for bootstrapping 
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
		tmp <- sampler(dfclean, "PlotID", reps = 100) # run at 100 for now, but need to bump this up to 1000 or 10,000 once things get all set up
		bigdata <- cbind(tmp, dfclean[tmp$RowID, ])
		
	# Steps: refit model to sample data
		# 1. Store original estimates from previous models. Will use these as 'start values' for bootstrap models
		# 2. Make a local cluster with 4 - nodes (# of processors on my laptop)
		# 3. Export data & lme4 on each cluster.
		# 4. Write a fxn to fit the model and return estimates. 
		# * the call glmer() is wrapped in try to avoid stopping the process if models do not converge. 
		
		f <- fixef(m_both) # store original parameters 
		r <- getME(m_both, "theta") 

		cl <- makeCluster(4)
		clusterExport(cl, c("bigdata", "f", "r"))
		clusterEvalQ(cl, require(lme4))
		
				myboot <- function(i) {
					object <- try(glmer(Surv ~ z + (1|SiteID/PlotID), data = bigdata, subset = Replicate == i, family = binomial,
					 nAGQ = 1, start = list(fixef = f, theta = r)), silent = TRUE)
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
		# results are fairly comparable, but possible to see a bit of bias. 

###########################################################################################
###########################################################################################
###########################################################################################

	# Graph out results to evaluate & check for other possible problems.  
		# conditional partial plots are fine, but run into the problem with MEMs of 
		# group level vs. average. Therefore we really want the average marginal probability.
		# will get the average marginal probability by plotting the average of each group - 
		# so have to caculate the conditional probabilities for each group and then average them.
		
		library(lme4) # should be already loaded 
		summary(dfclean$z); summary(dfclean$z1) # values range from 0 - 10. 
		jvalues <- with(dfclean, seq(from = min(z), to = max(z), length.out = 100)) # sequence of values to plot.
		tmpdat <- dfclean # copy the object over as back up, predictions will be added to this data frame.  
		
		# calculate predicted probabilities and store in a list
			pp <- lapply(jvalues, function(j) {
				tmpdat$z <- j
				predict(m_both, newdata = tmpdat, type = "response")
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
			
			# plot average marginal predicted probabilities
			ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_line() + ylim(c(0, 1))	

			# add on We could also add the lower and upper quartiles. This information shows us the range in which 50 percent of the predicted probabilities fell
			ggplot(plotdat, aes(x = z, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
				ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

###########################################################################################		
		# refit model with year & region & re-plot
			# covnert regions down to just three levels N, C & S and make year a facot variable 
			library(plyr)
			dfclean$Region <- revalue(dfclean$Region, c("C1"="C", "C2"="C", "C3"="C", "N1"="N", "N2"="N", "S1"="S", "S2"="S"))
			dfclean$Year <- as.factor(dfclean$Year)
			m_both <- glmer(Surv ~ z + Year + Region + (1|SiteID/PlotID), family=binomial(link="logit"), na.action=na.omit, data=dfclean)
				# re-run existing code
				tmpdat <- dfclean # copy the object over as back up, predictions will be added to this data frame. 
			
			# calculate predicted probabilities and store in a list
		
				# calculate predicted probabilities and store in a list
					biprobs <- lapply(levels(dfclean$Region), function(stage) {
					  tmpdat$Region[] <- stage
					  lapply(jvalues, function(j) {
						tmpdat$z <- j
						predict(m_both, newdata = tmpdat, type = "response")
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
						
				# plot it out 
				# graph it
					ggplot(plotdat2, aes(x = z, y = PredictedProbability)) +
					  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Region), alpha = .15) +
					  geom_line(aes(colour = Region), size = 2) +
					  ylim(c(0, 1))
				
				# It looks like this distribution is skewed a bit across each of the levels
					# can view distribution of predicted probabilities here:
					# center
					ggplot(data.frame(Probs = biprobs[[1]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75))
					# north
					ggplot(data.frame(Probs = biprobs[[2]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75))
					# south
					ggplot(data.frame(Probs = biprobs[[3]][[100]]), aes(Probs)) + geom_histogram() + scale_x_sqrt(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75))
					
	# even after log transformation our data is still pretty skewed. 				
					
				
###########################################################################################
###########################################################################################
###########################################################################################

# PROPER PROCEDURE FOR:
# Three level mixed effects logistic regression
# http://www.ats.ucla.edu/stat/r/dae/melogit.htm
	
	
	# estimate the model and store results in m_both2
			 m_both2 <- glmer(Surv ~ z + Year + Region + (1|SiteID) + (1|PlotID),
			  na.action=na.omit, data=dfclean, family = binomial, nAGQ=1)
		# check to see if this is the same as m_both with (1|SiteID/PlotID)
		m_both # nestedness defined 
		m_both2 # nestedness undefined - r figures out automatically, just as if random effects were crossed
			# yes it is so just use m_both
		# The standard deviation displayed is simply the square root of the variance
		# not the standard error of the estimate. 
		print(m_both, corr=FALSE)

		# Plot conditional modes PlotID & SiteID level effects 
			par(mfrow=c(1, 2))
			dotplot(ranef(m_both2, which = "PlotID", condVar = TRUE), scales = list(y = list(alternating = 0)))
			dotplot(ranef(m_both2, which = "SiteID", condVar = TRUE))

	###################################
	# RANDOM SLOPE: Try adding a random slope to the model. Can do so for just the SiteID
			 m_both3 <- glmer(Surv ~ z + Year + Region + (1 + z|SiteID) + (1|PlotID),
			  na.action=na.omit, data=dfclean, family = binomial, nAGQ=1)
	
			m_both3 # print results 
		
		dotplot(ranef(m_both3, which = "SiteID", condVar = TRUE), scales = list(y = list(alternating = 0)))
		dotplot(ranef(m_both3, which = "PlotID", condVar = TRUE), scales = list(y = list(alternating = 0)))

	##################################
	
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part 2: Choose the best model & then use bootstrapping for model inference 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
########################################################################################################################################	
	# PROPER RANDOM EFFECT STRUCTURE 
	 # Slightly more complicated because with GLMM, we can't just use REML to optimize the random 
	 # effect structure and component
		# possible alternative random effect structure: 
		# http://glmm.wikidot.com/faq#reml-glmml 
		# Can we drop interactive effect of PlotID nested within SiteID?
			# POSSIBLE ALTERNATIVES:
				#A - (1|SiteID)	random SiteID intercept
				#B - (z|SiteID) = (1+z|SiteID)	random slope of z within SiteID with correlated intercept
				#C - (0+z|SiteID) = (-1+z|SiteID)	random slope of z within SiteID: no variation in intercept
				#D - (1|SiteID) + (0+z|SiteID)	uncorrelated random intercept and random slope within SiteID
				#E - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
				#F - site+(1|site:PlotID)	fized effect of sites plus random variation in intercept among PlotIDs within sites
				#G - (z|site/PlotID) = (z|site)+(z|site:PlotID) = (1 + z|site)+(1+z|site:PlotID)	slope and intercept varying among sites and among PlotIDs within sites

		#A - (1|SiteID)	random SiteID intercept
			testA_full <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
			testA_reduced_just_site <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
			AIC(testA_full)	
			AIC(testA_reduced_just_site)
			anova(testA_full, testA_reduced_just_site)
			# Winner: testA_full
	
		#B - (z|SiteID) = (1+z|SiteID)	random slope of z within SiteID with correlated intercept
			testB_site_slope <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (z|SiteID), na.action=na.omit, data=dfclean, family = binomial)
			testB_reduced_just_site_int <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (1|SiteID), na.action=na.omit, data=dfclean, family = binomial)
			AIC(testB_site_slope)
			AIC(testB_reduced_just_site_int)
			AIC(testA_full)	
			# Winner: testA_full

		#C - (0+z|SiteID) = (-1+z|SiteID)	random slope of z within SiteID: no variation in intercept
		#D - (1|SiteID) + (0+z|SiteID)	uncorrelated random intercept and random slope within SiteID
			testC_site_slope_no_int <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (0 + z|SiteID), na.action=na.omit, data=dfclean, family = binomial)
			testD_site_slope_no_int <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (-1 + z|SiteID), na.action=na.omit, data=dfclean, family = binomial)
			AIC(testB_site_slope)
			AIC(testC_site_slope_no_int)
			AIC(testD_site_slope_no_int)
			# Winner: testA_full

		#E - (1|site/PlotID) = (1|site)+(1|site:PlotID)	intercept varying among sites and among PlotIDs within sites (nested random effects)
			AIC(testA_full)	
			# Winner: testA_full

		#F - site+(1|site:PlotID)	fized effect of sites plus random variation in intercept among PlotIDs within sites
			testF_with_site_fixed <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + SiteID + (1|PlotID), na.action=na.omit, data=dfclean, family = binomial)
			AIC(testA_full)	
			AIC(testF_with_site_fixed)	
			# model won't even converge 

		#G - (z|site/PlotID) = (z|site)+(z|site:PlotID) = (1 + z|site)+(1+z|site:PlotID)	slope and intercept varying among sites and among PlotIDs within sites
			testG_rand_int_slope_site_plot <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (z|SiteID/PlotID), na.action=na.omit, data=dfclean, family = binomial)
			AIC(testA_full)	
			AIC(testG_rand_int_slope_site_plot)	
			# this model with random slopes and intercepts is the winner, but hardly realistic 
			coef(mod.surv.glm)
			fixef(testA_full)
			# will go ahead with  'testA_full'. The model that has random intercept for the plot & site. 

		
########################################################################################################################################	
	# PROPER FIXED EFFECT STRUCTURE 	
	# code from: https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	
	
	
	
	
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Part 3: Model adequacy, check for violation of the basic underlying assumptions 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Chapter 21: GLMM applied on the spatial distribution of Koalas in a fragmented landscape
		# section 24.4.2 (pg. 483)
		# code from: https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
	
	library(glmmML)
	SiteID <- dfclean$SiteID
	glmmML(Surv ~ z + Year + Region, cluster = SiteID, na.action=na.omit, data=dfclean, family = binomial)
	
	# Possible models 
	 #mod1 <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + z*Year*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod2 <- glmer(Surv ~ z + Year + Region + z*Year + z*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod3 <- glmer(Surv ~ z + Year + Region + z*Year + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod4 <- glmer(Surv ~ z + Year + Region + z*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod5 <- glmer(Surv ~ z + Year + Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod6 <- glmer(Surv ~ z + Region + z*Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod7 <- glmer(Surv ~ z + Year + z*Year + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod8 <- glmer(Surv ~ z + Region + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod9 <- glmer(Surv ~ z + Year + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)
	mod10 <- glmer(Surv ~ z + (1|SiteID) + (1|PlotID) + (1|SiteID:PlotID), na.action=na.omit, data=dfclean, family = binomial)

	require(MuMIn)

	# use the mod.sel function to conduct model selection
	# and put output into object out.put
	out_put <- mod.sel(mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10)
	modList <- c("mod2","mod3","mod4","mod5","mod6","mod7","mod8","mod9","mod10")
	# what's it look like, hmm AIC with small sample bias adjustment AICc
	# delta AICc, and the model weights
	out_put

	# we can create a confidence set of models using the 
	# subset command E.G.,
	# select models with delta AICc less than 5
	# IMPORTANT: Weights have been renormalized!!
	subset(out_put, delta <5) # don't use these  values, just a quick check 
	
	# select models using Royall's 1/8 rule for strength of evidence
	# IMPORTANT: Weights have been renormalized!!
	subset(out_put, 1/8 < weight/max(out_put$weight))
	
	# select models 95% cumulative weight criteria 
	# IMPORTANT: Weights have been renormalized!! 
	subset(out_put, cumsum(out_put$weight) <= .95) 
		
		# coerce the object out.put into a data frame 
		# elements 6-10 in out.put have what we want 
		sel.table<-as.data.frame(out_put)[7:11] 
			
		# a little clean-up, lets round things a bit 
		sel.table[,2:3]<- round(sel.table[,2:3],2) 
		sel.table[,4:5]<- round(sel.table[,4:5],3) 
		# that’s better 
		sel.table 

		# how about a little renaming columns to fit proper conventions 
		# number of parameters (df) should be K 
		names(sel.table)[1] = "K" 

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
			
		#############################
		# 	Variable importance (based on AIC scores of all models with variable) 
		# Importance weights for individual predictor variables 
		# calculated using the importance function 
		importance(out_put)
	
		#############################
		# Model averaging 
			# since year is still significant in all top models, should do 
			# model averaging so we don't have to incorporate it into the IPM as separate runs 
			# of the IPM for each year. 
			
			# Model average using all candidate models, always use revised.var = TRUE 
			MA.ests<-model.avg(out_put, revised.var = TRUE) 
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
			write.csv(MA.est.table, "SurvModAvgEst.csv") 
			# Unfortunately estimating response values from these average 
			# parameters won't work for logistic regression.  
			
		#############################
		# Make predictions	
			#extract parameters and weights from confidence model set 
			# using get.models function 
			pred.parms<-get.models(out_put, subset= delta < 5) 
			# predict values using each model, here were just using the 
			# the example dataset, you could use a new dataset 
			model.preds = sapply(pred.parms, predict(type="response"), newdata = dfclean) 			
			
			#One word of warning, do not estimate model averaged parameters
			# for mixed models! You can, however, model average the predictions of GLMM. 	
				
		##########################################################	
		# Need to account for difference between glm & lm
		
			#All of the functions above can be used with objects created by
			#fitting linear models with the glm function. All of the above
			#applies to these glm objects. However, there is one important
			#difference that you should know. GLM’s such as Poisson regression
			#and Logistic regression can often fail to meet statistical
			#assumptions due to extra variance. This is often defined 
			#as over-dispersion (not an issue for normal linear regression!)
			#and requires the use of quasi-AIC for model selection.
		
		# Check for overdispersion 
			library(RVAideMemoire)
			overdisp.glmer(mod4)
			
		##########################################################	
		# Using MuMIn for glmer obejects
		dredge(mod4, rank = BIC) 
		
		#One word of warning, do not estimate model averaged parameters
		# for mixed models! You can, however, model average the predictions of GLMM. 
		
			

			
			