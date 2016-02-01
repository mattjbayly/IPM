#---
#title: "3. 03_Growth vital rate exploration of variance heterogeniity within the dataset 
#author: "Matthew Bayly"
#date: "Sunday, November 02, 2014"
#output: html_document
#---

# The follow script explores model structures with the inclusion and exclusion of 
# random effects & fixed effects for each of the four different vital rates. 

# all steps and procedures were followed from Zuur 2009. 
#Mixed Effects Models and Extensions in Ecology with R
# Authors: Zuur, Ieno, Walker, Saveliev and Smith. Publisher: Springer

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


#Open 2014 plant datafile 
setwd(path.dat)
# revised dataframe from previous script
d <- read.csv(file="Data_2014.csv")
dim(d); #colnames(plantdat)
# make site level factor
site <- levels(d$site); site <- as.factor(site) # study sites
# order by grouping
site <- factor(c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# second order polynomial for quadratic regression 
d$start2 <- d$start*d$start
d$size1_ln2 <- d$size1_ln*d$size1_ln


#============================================================================================#
# 2.  GROWTH - Exploring heterogeneity within the dataset 
#============================================================================================#

# The following protocol was adapted from Zurr 2009 Chapter 4 (pg. 90 - 99)
# STEP 1 : Start with full model with as many relevant fixed effects as possible & determine which homogeneity assumptions are violated.
# STEP 2 : Repeat step 1 using gls & REML
# STEP 3 : Plot & choose an appropriate variance structure. 
# STEP 4 : Fit a new gls with selected variance-covariance structure, full fixed parts. 
# STEP 5 : Compare new & old model with liklihood test (AIC), do residual still show signs of heterogenity
# STEP 6 : If residuals show heterogeneity then choose another variance structure; also try a different distribution.
# STEP 7 : Find optimal fixed component structure. 
# STEP 8 : Select top model 
# STEP 9 : Fit final model with REML, graphical validation ect. 
# STEP 10: Present final results in a table for presentation. 

##################################################################################################
# START Start with full model with as many relevant fixed effects as possible & determine which homogeneity assumptions are violated.
# work with gls models &
# choose appropriate variance structure
##################################################################################################
# GROWTH
	M1 <- lm(size2_ln ~ start*ENSEMBLE*moist_score*site_type*site, data=d, na.action=na.omit)
	par(mfrow=c(2,3))
	d2 <- d[complete.cases(d$size2_ln),]; d2 <- d2[complete.cases(d2$start),]
			plot(M1, which=c(1), col=1, add.smooth=FALSE, caption="", na.action=na.omit)
			plot(d2$start, resid(M1), xlab="start", ylab="Residuals")
			plot(d2$moist_score, resid(M1), xlab="moist_score", ylab="Residuals")
			plot(d2$ENSEMBLE, resid(M1), xlab="ENSEMBLE", ylab="Residuals")
			plot(d2$site_type, resid(M1), xlab="site_type", ylab="Residuals")
			plot(d2$site, resid(M1), xlab="site", ylab="Residuals")
				# higher variability for growth of small plants, wetter plots & one site. 
				# sites within the range are more variable

# Fixed variance structure / generalized least squared (gls).
	# "start" - Comparison of variance proportional to start size
		M.lm <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, data=d, na.action=na.omit) # gls=lm when no weight
		vf1Fixed <- varFixed(~start) # variance in model is proportional to starting size
		M.gls1 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, data=d, na.action=na.omit,
					weights=vf1Fixed)
		anova(M.lm, M.gls1) # from logLike, variance with size is not the better model. 
		AIC(M.lm, M.gls1) # can use AIC, because same number of parameters 
			# same general trend if moist_score is substituted
	
	# "site_type" - Comparison of variance proportional to site type (within/beyond)
		vf2 <- varIdent(form = ~ 1|site_type)
		M.gls2 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
						weights=vf2, data=d, na.action=na.omit)
		anova(M.lm, M.gls1, M.gls2) # from logLike ratio model with different variance per site type is better
		AIC(M.lm, M.gls1, M.gls2) # also possible here

				#Variance function: OUTPUT FOR (M.gls2)
				#Structure: Different standard deviations per stratum
				#Formula: ~1 | site_type 
				# Parameter estimates:
					#beyond     occupied   unoccupied  (M.gls2)
					#1.0000000  1.3680926  0.7557849 

	# Try accounting for site identity (Chap 4 pg. 77)
				#vf3 <- varIdent(form = ~ 1|site_type * factor(site))
				#		M.gls3 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
				#		weights=vf3, data=d, na.action=na.omit)
				#anova(M.gls2, M.gls3) # having site-type with site was better, probably because of Coast (hugely variable)

	# Which is better, incorporating diff variance spread per site & siteType, or just different spread per site-type
			# site_type OR site
			plot(M.lm, which=c(1), col=d$site_type, add.smooth=FALSE, caption="site_type", pch=19)
			plot(M.lm, which=c(1), col=as.factor(d$site), add.smooth=FALSE, caption="site", pch=19)
			# residual spread per site
			E <- resid(M.lm)
			coplot(E ~ start|site, data=d2)
			coplot(E ~ start|site_type, data=d2)

# Explore additional variance structures - or "variance covariates"
	# where variance is related to values of a covariate (size or site and or site_type)
	
	# The varPower variance structure (4.1.4 pg. 78 Zuur)
		# "start only" recall previously we used varFixed & varIdent
			vf3 <- varPower(form = ~ start)
			M.gls3 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf3, data=d, na.action=na.omit)
				# estimate: VarPower of start = 0.05335136 
			AIC(M.lm, M.gls1, M.gls2, M.gls3) # site type still works best; but gl3 beats gls1 & lm 
		
		# "start only | site_type" maybe with site type
			vf4 <- varPower(form = ~ start|site_type)
			M.gls4 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf4, data=d, na.action=na.omit)
						#	Parameter estimates for "start|site_type"
						#	beyond    occupied  unoccupied 
						#-0.02420274  0.10286843  0.17835435 
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4) # gls4 beats gls3 
				
	# The varExp Variance Structure (4.1.5 pg. 80 Zuur)
		# "start only" 
			vf5 <- varExp(form = ~ start)
			M.gls5 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf5, data=d, na.action=na.omit)
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5) 
			# varExp not working very well 
		# "start only | site_type" 
			vf5b <- varExp(form = ~ start|site_type)
			M.gls5b <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf5b, data=d, na.action=na.omit)
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5, M.gls5b) 			
	
	# The varConstPower Variance Structure (4.1.6 Zuur)
		# "start only" 
			vf6 <- varConstPower(form = ~ start)
			M.gls6 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf6, data=d, na.action=na.omit)
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5, M.gls5b, M.gls6) # worst yet for start size
		# "start only | site_type" 
			vf7 <- varConstPower(form = ~ start|site_type)
			M.gls7 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf7, data=d, na.action=na.omit)
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5, M.gls5b, M.gls6, M.gls7) 
			# this seems to be the best overall model, despite additional parameters
					#		beyond  occupied   unoccupied
					#const  0.34      0.83       0.00... (relative to 1)
					#power -0.04      0.21       0.20	 (relative to 0) ~ if constant =1 & power =0 then just back to lm

	# The varComb Variance Structure (4.1.7 Zuur)
			vf8 <- varComb(varIdent(form = ~ 1 | site_type))
			M.gls8 <- gls(size2_ln ~ start*ENSEMBLE*moist_score*site_type, 
							weights=vf8, data=d, na.action=na.omit)
			AIC(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5, M.gls5b, M.gls6, M.gls7, M.gls8) 	
				#varComb identical to varIdent in M.gls2
					    #beyond     occupied   unoccupied 
						#1.0000000  1.3680926  0.7557849 
						
	# Test competing models (watch out order matters here)
	anova(M.lm, M.gls1, M.gls2, M.gls3, M.gls4, M.gls5, M.gls5b, M.gls6, M.gls7, M.gls8) 	
				# & just the final models
	anova(M.gls2, M.gls7)
	# the gls7 (varConstPower) variance structure seems to do best, but only by <3 AIC points
	# var(eij) = var^2 * (constant1 + start|site_type^(power))^2
	# might be more complicated to work with
	# since groups site_type are factors & start was weak choose VarIdent here. 
		# allow for increase or decrease in variance across site types
		# interpretation of gls7 is not really relevant ecologically for study
	
		# OVERVIEW
			#VarFixed -  Fixed variance 	gls1
			#VarIdent - Different variances per stratum	gls2 ********
			#VarPower - Power of the variance covariate	gls3 & gls4
			#VarExp - Exponential of the variance covariate	gls5
			#VarConstPower - Constant plus power of the variance covariate	gls6
			#VarComb - A combination of variance functions	gls7
			
# graphical validation of final model output (4.1.9 Zuur)
		# ordinary residuals (observed - fitted values)
		# normalized residuals
		E1 <- resid(M.gls2)
		coplot(E1 ~ start|site_type, data=d2, ylab="Ordinary Residuals")
		# spread of residuals across site_types are allowed & now accounted for in model 
		E2 <- resid(M.gls2, type="normalized")
		coplot(E2 ~ start|site_type, data=d2, ylab="Normalized Residuals")

		
##################################################################################################	
# Find optimal fixed component structure. 
	library(nlme)
	#full fixed effects model formula for global model 
	anova(M.gls2) # crazy non-sense 4-way interactions still in model & can drop
	vf2 <- varIdent(form = ~ 1|site_type)

	# Test all two-way interactions for significant prior to exploring 3-way interactions
	# Forward Selection 
				fFull <- formula(size2_ln ~ start + ENSEMBLE + moist_score + site_type + 
								start:ENSEMBLE + start:moist_score + start:site_type +
								ENSEMBLE:moist_score + ENSEMBLE:site_type + 
								moist_score:site_type)
				M3.Full <- gls(fFull, weights=vf2, method="ML", data=d, na.action=na.omit) # use ML to compare fixed effects
			# Drop moist_score:site_type
				M3.Drop1 <- update(M3.Full, .~. -moist_score:site_type)
				anova(M3.Full, M3.Drop1) # no difference, later had  lower aic (less terms)
			# Drop ENSEMBLE:site_type
				M3.Drop2 <- update(M3.Full, .~. -ENSEMBLE:site_type)
				anova(M3.Full, M3.Drop2) # KEEP, interaction significant *********
			# Drop ENSEMBLE:moist_score
				M3.Drop3 <- update(M3.Full, .~. -ENSEMBLE:moist_score)
				anova(M3.Full, M3.Drop3) # no difference, drop terms
			# Drop start:site_type
				M3.Drop4 <- update(M3.Full, .~. -start:site_type)
				anova(M3.Full, M3.Drop4) # marginal significance (.), but drop as well given multiple series of test
			# Drop start:moist_score
				M3.Drop5 <- update(M3.Full, .~. -start:moist_score)
				anova(M3.Full, M3.Drop5) # no difference 
			# Drop start:ENSEMBLE
				M3.Drop6 <- update(M3.Full, .~. -start:ENSEMBLE)
				anova(M3.Full, M3.Drop6) # no difference 
			# Only keep two-way interaction for ENSEMBLE:site_type & start:site_type			
		
	# Backwards Selection 
			# Update model, 
			fReduced <- formula(size2_ln ~ start + ENSEMBLE + moist_score + site_type + 
								ENSEMBLE:site_type)	
			M4.Full <- gls(fReduced, weights=vf2, method="ML", 
								data=d, na.action=na.omit) 
			# Drop ENSEMBLE:site_type
				M4.Drop1 <- update(M4.Full, .~. -ENSEMBLE:site_type)
				anova(M4.Full, M4.Drop1) # keep interaction in 
			# Drop moist_score
				M4.Drop2 <- update(M4.Full, .~. -moist_score)
				anova(M4.Full, M4.Drop2) # keep moisture in
			# Drop quadratic term
				M4.Drop2 <- update(M4.Full, .~. -moist_score)
				anova(M4.Full, M4.Drop2) # keep moisture in
			# Dose not look like any further terms can be doped. 
			
	# Test quadratic term
			quad <- gls(size2_ln ~ start + I(start^2) + ENSEMBLE + moist_score + site_type + 
								ENSEMBLE:site_type, weights=vf2, method="ML", 
								data=d, na.action=na.omit)
			no_quad <- gls(size2_ln ~ start + ENSEMBLE + moist_score + site_type + 
								ENSEMBLE:site_type, weights=vf2, method="ML", 
								data=d, na.action=na.omit)
			anova(quad, no_quad)
			# Don't include quadratic term
	# Test quadratic term
			MFinal <- M4.Full
			E <- resid(MFinal, type = "normalized")
			Fit <- fitted(MFinal)
			op <- par(mfrow = c(1, 2))
			plot(x = Fit, y = E, xlab = "Fitted values", ylab = "Residuals", main = "Residuals versus fitted values")
			#identify(Fit, E) # to choose & identify particular point in plot
			hist(E, nclass = 15, xlab="Residuals")
			par(op)
						# temp fxn for na's
						completeFun <- function(data, desiredCols) {
						completeVec <- complete.cases(data[, desiredCols])
						return(data[completeVec, ])}	
				dtmp <- completeFun(d, c("size2_ln", "start"))
			# box plot of site types
			boxplot(predict(MFinal) ~ site_type, data = dtmp)
			
###########################################################################################		
				
			