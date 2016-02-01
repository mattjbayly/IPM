#---
#title: "4_ Exploration of Mixed Effect Models and Nested structure (Zuur Chapter 5)
#author: "Matthew Bayly"
#date: "Sunday, January 02, 2014"
#---

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

# Temporary Growth dataframe with no NA's 
			completeFun <- function(data, desiredCols) {
			completeVec <- complete.cases(data[, desiredCols])
			return(data[completeVec, ])}	
			dtmp <- completeFun(d, c("size2_ln", "start"))

	

#============================================================================================#
# 1. Explore basic growth model 
#============================================================================================#

library(nlme)
# For growth
		Mlme1 <- lme(size2_ln ~ start, random= ~1|site/Uplot, 
						data=d, na.action=na.omit)
	#summary(Mgr)
		# random intercept: site level growth plot - random intercept 
		F0 <- fitted(Mlme1, level = 0)
		F1 <- fitted(Mlme1, level = 1)
		I <- order(dtmp$start); starts <- sort(dtmp$start)
		plot(starts, F0[I], lwd = 4, type = "l",
		ylim = c(0, 6.5), ylab = "Final", xlab = "Initial")
		for (i in 1:length(site)){
				current <- site[i]
				x1 <- dtmp$start[dtmp$site == current]
				y1 <- F1[dtmp$site == current]
				K <- order(x1)
				lines(sort(x1), y1[K])
			}
		points(d$start, d$size2_ln, col=as.factor(d$site), cex=0.3, pch=19)
		
		# random slope & intercept (ignor plot nesting now for convergence)
		Mlme1 <- lme(size2_ln ~ start, random= ~ 1 + start|site, 
						data=d, na.action=na.omit)		
		F0 <- fitted(Mlme1, level = 0)
		F1 <- fitted(Mlme1, level = 1)
		I <- order(dtmp$start); starts <- sort(dtmp$start)
		plot(starts, F0[I], lwd = 4, type = "l",
		ylim = c(0, 6.5), ylab = "Final", xlab = "Initial")
		for (i in 1:length(site)){
				current <- site[i]
				x1 <- dtmp$start[dtmp$site == current]
				y1 <- F1[dtmp$site == current]
				K <- order(x1)
				lines(sort(x1), y1[K])
			}
		points(d$start, d$size2_ln, col=as.factor(d$site), cex=0.3, pch=19)
				

#============================================================================================#
# GROWTH - The good approach to model selection (from Zuur Chapter 5, Section 5.8.2)
#============================================================================================#
	
# For the growth model Plots 1 - 5 at Coast Fork are really screwing things up
# they always act as huge residual leverage points. 
# what happened was plants survived but no growth occurred for some strange reason. 
# Waterfowl herbivory? Not clear, but only plants in whole project to have done this. 
# Therefore should exclude for growth function. We saw in the last script 03_Growth_var_heterogeneity_explor
# that these plots alone drove site-type differences in variance structure. 
# we will exclude them temporary and reference a new dataframe called dg for d(growth estimate) 

dg<-d[!(d$Uplot=='COASTP01' | d$Uplot=='COASTP02' | d$Uplot=='COASTP03' | d$Uplot=='COASTP04' | d$Uplot=='COASTP05'),]
# & complete cases for predict function
dg <- completeFun(dg, c("size2_ln", "start"))

# Start with all key relevant fixed effects and their interactions 
# work through to get appropriate random structure and then 
# work through appropriate fixed effects for final model 
		
	#5.10.1 Step 1 of the Protocol: Linear Regression		
	#5.10.2 Step 2 of the Protocol: Fit the Model with GLS	
	#5.10.3 Step 3 of the Protocol: Choose a Variance Structure
	#5.10.4 Step 4: Fit the Model
    #5.10.5 Step 5 of the Protocol: Compare New Model with Old Model
	#5.10.6 Step 6 of the Protocol: Everything Ok?
	#5.10.7 Steps 7 and 8 of the Protocol: The Optimal Fixed Structure
	#5.10.8 Step 9 of the Protocol: Refit with REML and Validate the Model
	#5.10.9 Step 10 of the Protocol
	
	#5.10.1 Step 1 of the Protocol: Linear Regression		
			library(AED)
			M.lm <- lm(size2_ln ~ start + moist_score + site_type + ENSEMBLE, 
						data=dg, na.action=na.omit)
					plot(M.lm, select=c(1))
			# weak evidence of variance heterogeneity
			E <- rstandard(M.lm)
			boxplot(E ~ site_type, data = dg,
					ylim = c(-3, 3))
			abline(0,0); axis(2)						
			
	#5.10.2 Step 2 of the Protocol: Fit the Model with GLS	
			library(nlme)
			Form <- formula(size2_ln ~ start + moist_score + site_type + ENSEMBLE + site_type:ENSEMBLE + site_type:moist_score)
			M.gls <- gls(Form, data=dg)
	
	#5.10.3 Step 3 of the Protocol: Choose a Variance Structure
		# need to explore this again with COAST plots 1 - 5 dropped. 
		# quick check showed gls2 (VarIdent - site_type) to still be the best
	
	#5.10.4 Step 4: Fit the Model
   		# plots nested within sites (random intercept)
			# plot nested within site
			M1.lme <- lme(Form, random= ~ 1|site/Uplot, method="REML", data=dg)
			# site level nesting
			M1b.lme <- lme(Form, random= ~ 1|site, method="REML", data=dg)
			# assume all plots are totally independent
			M1c.lme <- lme(Form, random= ~ 1|Uplot, method="REML", data=dg)

    #5.10.5 Step 5 of the Protocol: Compare New Model with Old Model
			# can use anova since both gls & lme used REML
			anova(M.gls, M1.lme, M1b.lme, M1c.lme)
			anova(M1c.lme, M1.lme) # ignoring site-level produces a slighly better model 
						# but not significant and AIC only 2 points different
						# will run into problems down the road if we treat all plots as independant
						# so keep multilevel nesting for now. 

			anova(M.gls, M1.lme)
			# Random intercept way better! AIC different by allot. 
			# L = 144.2 (df = 9, p < 0.0001)
			
	#5.10.6 Step 6 of the Protocol: Everything Ok?
			E2 <- resid(M1.lme, type = "normalized")
			F2 <- fitted(M1.lme)
			op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
			MyYlab <- "Residuals"
			plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
			boxplot(E2 ~ site_type, data = dg,
				main = "site_type", ylab = MyYlab)
			plot(x = dg$start, y = E, ylab = MyYlab,
			main = "start", xlab = "start")
			plot(x = dg$ENSEMBLE, y = E, ylab = MyYlab,
			main = "m", xlab = "ENSEMBLE")
			par(op)
		# Some heteroscedasticity still noticed for start size - smaller more variable
		# and not accounted for in current nested MEM model structure. 
	
	#5.10.7 Steps 7 and 8 of the Protocol: The Optimal Fixed Structure
			summary(M1.lme) # don't use anova here b/c relies on sequential testing
			# site_type & ENSEMBLE relationships weak to non-existent in data
			
		# Drop terms & test 
			# be sure to change method to ML to compare models with different fixed structures
			M1.Full <- lme(Form, random= ~ 1|site/Uplot, method="ML", data=dg)
			M1.A <- update(M1.Full, .~. -site_type:ENSEMBLE)
			anova(M1.Full, M1.A) # drop interaction
			M1.B <- update(M1.Full, .~. -site_type:moist_score)
			anova(M1.Full, M1.B) # drop interaction
			
			# update form 
			Form <- formula(size2_ln ~ start + moist_score + site_type + ENSEMBLE)
			M2.Full <- lme(Form, random= ~ 1|site/Uplot, method="ML", data=dg)
			M2.A <- update(M2.Full, .~. -ENSEMBLE)
			anova(M2.Full, M2.A) # drop ENSEMBLE score
			M2.B <- update(M2.Full, .~. -site_type)
			anova(M2.Full, M2.B) # keep site_type
			M2.C <- update(M2.Full, .~. -moist_score)
			anova(M2.Full, M2.C) # keep moist_score
			M2.D <- update(M2.Full, .~. -start)
			anova(M2.Full, M2.D) # keep start
			
			# just one more check for the possibility of fitting a quadratic function for growth 
			no_quad <- lme(size2_ln ~ start + moist_score + site_type, random= ~ 1|site/Uplot, method="ML", data=dg)
			quad <- lme(size2_ln ~ start + I(start^2) + moist_score + site_type, random= ~ 1|site/Uplot, method="ML", data=dg)
			anova(no_quad, quad) # no.. its really not improving things at all 
			# growth is linear 
				
	#5.10.8 Step 9 of the Protocol: Refit with REML and Validate the Model
		# FINAL MODEL
		M5 <- lme(size2_ln ~ start + moist_score + site_type, random= ~ 1|site/Uplot, method="REML", data=dg) 
		# need to change method back to REML 
		summary(M5)
		# correlation with intercept is weak (good sign)
		# Visual plot to check assumption & validity of model 
				E2 <- resid(M5, type = "normalized")
				F2 <- fitted(M5)
				op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
				MyYlab <- "Residuals"
				plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
				boxplot(E2 ~ site_type, data = dg,
					main = "site_type", ylab = MyYlab)
				plot(x = dg$start, y = E, ylab = MyYlab,
				main = "start", xlab = "start")
				plot(x = dg$moist_score, y = E, ylab = MyYlab,
				main = "m", xlab = "moist_score")
				par(op)
	
		# Variance explained, correlation at level 
		# ICC plot = variance plot / (variance plot + variance site + variance error)
			h <- VarCorr(M5); den <- sum(as.numeric(c(h[2], h[4], h[5]))) # denominator
			as.numeric(h[4])/den
		# ICC site = variance site / (variance plot + variance site + variance error)
			as.numeric(h[2])/den
		
	#5.10.9 Step 10 of the Protocol - biological interpretation 
		# growth primarily determined by start potential and moisture
		# ENSEMBLE score had no relationship 
		# unoccupid sites did not grow as much as occupied site or 
		# those beyond the range, distinguishing between occupied and 
		# and sites beyond the range was not possible. 
##############################################################################################

# Should multiple variance structure based on site_type be included into growth model?
# Following code from Chapter 20 in Zuur 2009 
		MFinal_b <- lme(size2_ln ~ start + moist_score + site_type, 
				random= ~ 1|site/Uplot, method="REML", data=dg,
				weights = varIdent(form = ~ 1 |site_type)) 

		anova(M5, MFinal_b) # dosn't make a differences, even if we just substitute for site
		# try vaPower & a few others prospective candidates just to be sure
		MFinal_varPower <- lme(size2_ln ~ start + moist_score + site_type, random= ~ 1|site/Uplot, method="REML", data=dg, weights = varPower(form = ~ start|site_type))
		anova(M5, MFinal_varPower) #
		# varpower ~ weights = varPower(form = ~ start|site_type) seems to have a better fit for growth 
		# concerning that varIdent didn't 
		
# PLOT RESPONSE
	library(sjPlot) # refit in lmer to plot
	M5_re <- lmer(size2_ln ~ start + moist_score + site_type + (1|site/Uplot), REML=TRUE, data=dg) 
	sjp.lmer(M5_re, type = "fe.ri", vars = "start")
	

#============================================================================================#
# PLOTING - Partial dependance plots
#============================================================================================#

		library(effects)
	    par(mfrow=c(1,3))
		plot(Effect("start", M5), main="", xlab="Start Size", ylab="Fall Size")
		plot(Effect("moist_score", M5), main="", xlab="Plot Moisture", ylab="Fall Size")
		plot(Effect("site_type", M5), main="", xlab="", ylab="Fall Size", cex=3, lwd=3, cex.axis = 2, cex.lab=2)
	
	M5_re <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|site/Uplot), REML=TRUE, data=dg) 
		plot(Effect("ENSEMBLE", M5_re), main="", xlab="SDM", ylab="Fall Size")

		
		
