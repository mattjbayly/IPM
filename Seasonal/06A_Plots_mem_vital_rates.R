#---
#title: "6_ Plots of vital rates in mixed effect modelling framework
#author: "Matthew Bayly"
#date: "Sunday, January 02, 2014"
#---

# Correctly plotting mems is tricky 


# INDEX
##
##
##
##
##
##
##
##

#============================================================================================#
# 1. START: SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS
#============================================================================================#

### _Set directories for computer_ ###
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
setwd(path.set)
source("00_SetDirectories.R") # directory script (edit for your own computer). 
setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)

# libraries
library(lme4)
library(nlme)
library(effects)
library(AED)
library(glmmML)
library(sjPlot)

#Open 2014 plant datafile 
setwd(path.dat)
# revised dataframe from previous script
d <- read.csv(file="Data_2014.csv")
dim(d); #colnames(plantdat)
# make site level factor
site <- levels(d$site); site <- as.factor(site) # study sites
# order by grouping
site <- factor(c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# Temporary Growth dataframe with no NA's 
			completeFun <- function(data, desiredCols) {
			completeVec <- complete.cases(data[, desiredCols])
			return(data[completeVec, ])}	
			dtmp <- completeFun(d, c("size2_ln", "start"))

#============================================================================================#
# SURVIVORHSIP - VISUAL PLOTs
#============================================================================================#

	# VISUAL PLOTs
			dtmp <- completeFun(d, c("surv_end", "start"))
			dtmp2 <- dtmp
			#dtmp2$moist_score <- scale(dtmp$moist_score, center = TRUE, scale = TRUE)
			#dtmp2$start <- scale(dtmp$start, center = TRUE, scale = TRUE)
			#dtmp2$ENSEMBLE <- scale(dtmp$ENSEMBLE, center = TRUE, scale = TRUE)

			SurBPSW <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=dtmp2)	

			
					# useful extraction calls	
					# from: http://jaredknowles.com/journal/2014/5/17/mixed-effects-tutorial-2-fun-with-mermod-objects
					class(SurBPSW)
					VarCorr(SurBPSW)
					attr(VarCorr(SurBPSW), "sc")
					slotNames(SurBPSW)
					SurBPSW@call
					methods(class = "merMod")
					isREML(SurBPSW); isGLMM(SurBPSW)
					fixef(SurBPSW) # extracts fixed effects
					confint(SurBPSW, level = 0.95, method="Wald") # rough bootst confint (fast)
					sigma(SurBPSW) # residual standard error
					#standardize fixed effects into “effect sizes” by dividing the fixed effect paramters by the residual standard error,
					fixef(SurBPSW) / sigma(SurBPSW)
					#ranef(SurBPSW)
					re1 <- ranef(SurBPSW, condVar=TRUE) # save the ranef.mer object
							class(re1)
					#attr(re1[[1]], which = "postVar")
					re1 <- ranef(SurBPSW, condVar=TRUE, whichel = "site")
					print(re1)
					dotplot(re1)
					
					# plot for variable values 
					names <- c("(Intercept)", "start", "moist_score", "ENSEMBLE", "site_typeoccupied", "site_typeunoccupied")
					frame <- cbind(names, data.frame(as.matrix(fixef(SurBPSW))), data.frame(as.matrix(confint(SurBPSW, level = 0.95, method="Wald"))))
					colnames(frame) <- c( 'VAR', 'EST', 'LOW', 'UP')
					plot(frame$EST, frame$VAR, xlim=c(-2,2), pch=19)
					
				##################################################################################
					# Plotting tutorial from 	
					#http://www.ats.ucla.edu/stat/r/dae/melogit.htm
					# confidence intervals of fixed effects
					se <- sqrt(diag(vcov(SurBPSW)))	
					(tab <- cbind(Est = fixef(SurBPSW), LL = fixef(SurBPSW) - 1.96 * se, UL = fixef(SurBPSW) + 1.96 * se))
					# for odds ratio (rather than logit scale)
					exp(tab)
			
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(start), to = max(start), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$start <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "start")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = start, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = start, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$start <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "start")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = start, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p1 <- ggplot(plotdat2, aes(x = start, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("Start Size") +
								ylab("Survivorship") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################
		# SOIL MOISTURE PLOT
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(moist_score), to = max(moist_score), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p2 <- ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("Plot-level Moisture") +
								ylab("Survivorship") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################						  
########################################################################################################################

		library(ggplot2)
		library(grid)
		library(gridExtra)
		grid.arrange(p1, p2, ncol = 2, main = "SURVIVORSHIP")			  
	


########################################################################################################################
		# ENSEMBLE score plot		
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(ENSEMBLE), to = max(ENSEMBLE), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$ENSEMBLE <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "ENSEMBLE")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							
							# graph it - type 1
							p3 <- ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = .15) +
							  geom_line(size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("SDM Ensemble Score") +
								ylab("Survivorship") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
							
			grid.arrange(p1, p2, p3, ncol = 3, main = "SURVIVORSHIP")			  

			
			
#============================================================================================#
# pFlower - VISUAL PLOTs
#============================================================================================#
		
		
		# VISUAL PLOTs
			#dtmp2$moist_score <- scale(dtmp$moist_score, center = TRUE, scale = TRUE)
			#dtmp2$start <- scale(dtmp$start, center = TRUE, scale = TRUE)
			#dtmp2$ENSEMBLE <- scale(dtmp$ENSEMBLE, center = TRUE, scale = TRUE)
		
		dtmp <- completeFun(d, c("size2_ln", "pFlower"))
		dtmp2 <- dtmp
		
			SurBPSW <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=dtmp2)	

				##################################################################################
					# Plotting tutorial from 	
					#http://www.ats.ucla.edu/stat/r/dae/melogit.htm
					# confidence intervals of fixed effects
					se <- sqrt(diag(vcov(SurBPSW)))	
					(tab <- cbind(Est = fixef(SurBPSW), LL = fixef(SurBPSW) - 1.96 * se, UL = fixef(SurBPSW) + 1.96 * se))
					# for odds ratio (rather than logit scale)
					exp(tab)
			
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(size2_ln), to = max(size2_ln), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$size2_ln <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "size2_ln")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = size2_ln, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = size2_ln, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$size2_ln <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "size2_ln")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = size2_ln, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p1 <- ggplot(plotdat2, aes(x = size2_ln, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("Fall Size") +
								ylab("pFlower") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################
		# SOIL MOISTURE PLOT
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(moist_score), to = max(moist_score), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p2 <- ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("Plot-level Moisture") +
								ylab("pFlower") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################						  
########################################################################################################################

		library(ggplot2)
		library(grid)
		library(gridExtra)
		grid.arrange(p1, p2, ncol = 2, main = "")			  
	


########################################################################################################################
		# ENSEMBLE score plot		
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(ENSEMBLE), to = max(ENSEMBLE), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$ENSEMBLE <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "ENSEMBLE")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							
							# graph it - type 1
							p3 <- ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = .15) +
							  geom_line(size = 2) +
							  ylim(c(0, 1)) + 
							    xlab("SDM Ensemble Score") +
								ylab("pFlower") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
							
			grid.arrange(p1, p2, p3, ncol = 3, main = "Prob FLOWER")			  

	

			
			
#============================================================================================#
# FECUNDITY - VISUAL PLOTs
#============================================================================================#
		
		
		# VISUAL PLOTs
			#dtmp2$moist_score <- scale(dtmp$moist_score, center = TRUE, scale = TRUE)
			#dtmp2$start <- scale(dtmp$start, center = TRUE, scale = TRUE)
			#dtmp2$ENSEMBLE <- scale(dtmp$ENSEMBLE, center = TRUE, scale = TRUE)
		
		dtmp <- completeFun(d, c("size2_ln", "fec"))
		dtmp2 <- dtmp
	
		
			SurBPSW <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family = poisson, na.action=na.omit, data=dtmp2)	

				##################################################################################
					# Plotting tutorial from 	
					#http://www.ats.ucla.edu/stat/r/dae/melogit.htm
					# confidence intervals of fixed effects
					se <- sqrt(diag(vcov(SurBPSW)))	
					(tab <- cbind(Est = fixef(SurBPSW), LL = fixef(SurBPSW) - 1.96 * se, UL = fixef(SurBPSW) + 1.96 * se))
					# for odds ratio (rather than logit scale)
					exp(tab)
			
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(size2_ln), to = max(size2_ln), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$size2_ln <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "size2_ln")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = size2_ln, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = size2_ln, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$size2_ln <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "size2_ln")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = size2_ln, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p1 <- ggplot(plotdat2, aes(x = size2_ln, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							    xlab("Fall Size") +
								ylab("Fecundity") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################
		# SOIL MOISTURE PLOT
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(moist_score), to = max(moist_score), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							ggplot(plotdat, aes(x = moist_score, y = PredictedProbability)) + geom_linerange(aes(ymin = Lower,
							ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))

							# calculate predicted probabilities and store in a list
							biprobs <- lapply(levels(dtmp$site_type), function(stage) {
							  dtmp$site_type[] <- stage
							  lapply(jvalues, function(j) {
								dtmp$moist_score <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							  })
							})

							# get means and quartiles for all jvalues for each level of site_type
							plotdat2 <- lapply(biprobs, function(X) {
							  temp <- t(sapply(X, function(x) {
								c(M=mean(x), quantile(x, c(.25, .75)))
							  }))
							  temp <- as.data.frame(cbind(temp, jvalues))
							  colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "moist_score")
							  return(temp)
							})

							# collapse to one data frame
							plotdat2 <- do.call(rbind, plotdat2)

							# add cancer stage
							plotdat2$site_type <- factor(rep(levels(dtmp$site_type), each = length(jvalues)))

							# show first few rows
							head(plotdat2)

							# graph it - type 1
							ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .15) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1)) + facet_wrap(~  site_type)

							 # graph it - type 3 ~ just lines
							p2 <- ggplot(plotdat2, aes(x = moist_score, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							    xlab("Plot-level Moisture") +
								ylab("Fecundity") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################################						  
########################################################################################################################

		library(ggplot2)
		library(grid)
		library(gridExtra)
		grid.arrange(p1, p2, ncol = 2, main = "")			  
	


########################################################################################################################
		# ENSEMBLE score plot		
						# sample plot 
						jvalues <- with(dtmp, seq(from = min(ENSEMBLE), to = max(ENSEMBLE), length.out = 100))
							# calculate predicted probabilities and store in a list
							pp <- lapply(jvalues, function(j) {
								dtmp$ENSEMBLE <- j
								predict(SurBPSW, newdata = dtmp, type = "response")
							})

							# get the means with lower and upper quartiles
							plotdat <- t(sapply(pp, function(x) {
										c(M = mean(x), quantile(x, c(0.25, 0.75)))
									}))

							# add in LengthofStay values and convert to data frame
							plotdat <- as.data.frame(cbind(plotdat, jvalues))

							# better names and show the first few rows
							colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "ENSEMBLE")
							head(plotdat)
						
							# plot average marginal predicted probabilities
							library(ggplot2)
							ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) + geom_line() +
									ylim(c(0, 1))
							
							# graph it - type 1
							p3 <- ggplot(plotdat, aes(x = ENSEMBLE, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = .15) +
							  geom_line(size = 2) +
							    xlab("SDM Ensemble Score") +
								ylab("Fecundity") +
							  guides(fill=guide_legend(title=NULL)) +
							  theme(legend.title=element_blank()) +	
							  theme_bw() + 
							   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
									panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
							
			grid.arrange(p1, p2, p3, ncol = 3, main = "FECUNDITY")			  

#============================================================================================#

#============================================================================================#

