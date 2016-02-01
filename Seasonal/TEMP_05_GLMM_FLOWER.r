#---
#title: "5_ Exploration of Mixed Effect Models through GLMM for flowering, fruit & fecundity vital rates
#author: "Matthew Bayly"
#date: "Sunday, January 02, 2014"
#---

#"We use the glmmML function here because it estimates the model
#parameters by maximum likelihood and allows AICs to be calculated" 
# "An alternative would be to use the lmer function in the package 
# lme4 with the Lapacianor adaptive Gauss-Hermite methods"- Zuur


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
# SURVIVORHSIP - Following procedure from Zuur 2009 Chapter 21 - GLMM case study pg. 469 - 497
#============================================================================================#
	
	# glmmML cannot work with nested clusters site/Uplot so will use unique plot ID & revist with glmer
		Glmm_1 <- glmmML(surv_end ~  start + moist_score + site_type + ENSEMBLE + 
				site_type:ENSEMBLE + site_type:moist_score + ENSEMBLE:moist_score,
				cluster = Uplot, data = d, family = binomial, na.action=na.omit)

	# run competing models & rank based on Akaike weights
	#B-basic; P - Plot level; S- SDM; W - within/beyond
	library(qpcR)

		# SURVIVAL 
			SurB <- glmmML(surv_end ~ start, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBP <- glmmML(surv_end ~ start + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBS <- glmmML(surv_end ~ start + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBW <- glmmML(surv_end ~ start + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBWP <- glmmML(surv_end ~ start + moist_score + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBWS <- glmmML(surv_end ~ start + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBSP <- glmmML(surv_end ~ start + moist_score + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
			SurBPSW <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	
			SurBPSW_P.S. <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type + moist_score:ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	
			SurBPSW_P.W. <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type + moist_score:site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	
			SurBPSW_S.W. <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type + ENSEMBLE:site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	

			names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW", "SurBPSW_P.S.", "SurBPSW_P.W.", "SurBPSW_S.W.")
			my_aic <- c(SurB$aic, SurBP$aic, SurBS$aic, SurBW$aic, SurBWP$aic, SurBWS$aic, SurBSP$aic, SurBPSW$aic, SurBPSW_P.S.$aic, SurBPSW_P.W.$aic, SurBPSW_S.W.$aic)
			library(qpcR)
			moo <- akaike.weights(my_aic)
			tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
			tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')
		##########
	#============================================================================================#
			# ALTERNATIVE WITH LMER

			# SURVIVAL 
			SurB <- glmer(surv_end ~ start + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBP <- glmer(surv_end ~ start + moist_score + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBS <- glmer(surv_end ~ start + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBW <- glmer(surv_end ~ start + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBWP <- glmer(surv_end ~ start + moist_score + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBWS <- glmer(surv_end ~ start + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBSP <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
			SurBPSW <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	
			SurBPSW_P.S. <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + moist_score:ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)		
			SurBPSW_P.W. <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + moist_score:site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)		
			SurBPSW_S.W. <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + ENSEMBLE:site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)		

			names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW", "SurBPSW_P.S.", "SurBPSW_P.W.", "SurBPSW_S.W.")
			my_aic <- c(AIC(SurB), AIC(SurBP), AIC(SurBS), AIC(SurBW), AIC(SurBWP), AIC(SurBWS), AIC(SurBSP), AIC(SurBPSW), AIC(SurBPSW_P.S.), AIC(SurBPSW_P.W.), AIC(SurBPSW_S.W.))
			library(qpcR)
			moo <- akaike.weights(my_aic)
			tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
			tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

	#============================================================================================#		
	#============================================================================================#
		# EXPLORE MODEL ADEQUECY (ok?)
			#"Traditionally, the fit of logistic regression models have been assessed using
			#global goodness-of-fit tests based on the deviance or Pearson χ2 statistics" - but sketchy here
	
		# Zuur 2009 - 21.4.2 Model Adequacy pg. 487
		# Approach: quantile-quantile plots & partial residual plots
		# load special Zuur functions
		setwd(path.funct); source("zuur_chap21_qq.R")
						# make temporary df with no NA's
									completeFun <- function(data, desiredCols) {
										completeVec <- complete.cases(data[, desiredCols])
										return(data[completeVec, ])}	
					# Zuur fxn's cant handle NA's
					dtmp <- completeFun(d, c("surv_end", "start"))
						
					# fit the best full model 
					Model_Best <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=dtmp)	
					Fitted<-fitted.glmmML(model=Model_Best,data=dtmp)
					#calculate the ordered residuals for the model
					Resids<-sort(res.glmmML(model=Model_Best,data=dtmp))
					Resids_0<-Resids[Resids<0]
					Resids_1<-Resids[Resids>=0]
					#specify the number of replicates
					Reps<-1000
					#simulate Reps data sets
					Sims<-matrix(rbinom(n=length(Fitted)*Reps,size=1,p=Fitted),nrow=length(Fitted),ncol=Reps)
					#fit the model to each simulated data set
					# should take ~ 2 mins
					Models<-apply(Sims,MARGIN=2,FUN=function(X){return(list(X,glmmML(formula=surv_end ~ start + moist_score + ENSEMBLE + site_type,cluster=Uplot,data=dtmp,family=binomial)))})
					
					#calculate the ordered (interpolated) simulated residuals and the point-wise 95% confidence intervals
					Resids_Sim<-matrix(unlist(lapply(Models,FUN=function(X){TempData<-dtmp;TempData[,"surv_end"]<-X[[1]];return(sort(res.glmmML(X[[2]],TempData)))})),ncol=Reps,nrow=nrow(dtmp))
					Resids_Sim_0<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X<0],ppoints(Resids_0,a=1))}),ncol=Reps,nrow=length(Resids_0))
					Resids_Sim_1<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X>=0],ppoints(Resids_1,a=1))}),ncol=Reps,nrow=length(Resids_1))
					Resids_Sim_0_Median<-apply(Resids_Sim_0,MARGIN=1,FUN=median)
					Resids_Sim_1_Median<-apply(Resids_Sim_1,MARGIN=1,FUN=median)
					Resids_Sim_0_Lower<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.025)})
					Resids_Sim_0_Upper<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.975)})
					Resids_Sim_1_Lower<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.025)})
					Resids_Sim_1_Upper<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.975)})
					#plot the qauntile-quantile plot with 95% confidence intervals and 1:1 line
					plot(Resids_Sim_0_Median, Resids_0, xlim=c(-1,1), ylim=c(-1,1), xlab="simulated quantiles", ylab="fitted quantiles")
					points(Resids_Sim_1_Median,Resids_1)
					lines(Resids_Sim_0_Median,Resids_Sim_0_Lower)
					lines(Resids_Sim_0_Median,Resids_Sim_0_Upper)
					lines(Resids_Sim_1_Median,Resids_Sim_1_Lower)
					lines(Resids_Sim_1_Median,Resids_Sim_1_Upper)
					abline(0,1,lty=3)

		# Points are a little off better explor variable by variable to see whats going on.
					#calculate the partial residuals for each covariate
					Part_Res1<-res.glmmML(Model_Best,dtmp)/(fitted.glmmML(Model_Best,dtmp)*(1-fitted.glmmML(Model_Best,dtmp)))+Model_Best$coefficients[2]*dtmp[,"start"]
					Part_Res2<-res.glmmML(Model_Best,dtmp)/(fitted.glmmML(Model_Best,dtmp)*(1-fitted.glmmML(Model_Best,dtmp)))+Model_Best$coefficients[3]*dtmp[,"moist_score"]
					Part_Res3<-res.glmmML(Model_Best,dtmp)/(fitted.glmmML(Model_Best,dtmp)*(1-fitted.glmmML(Model_Best,dtmp)))+Model_Best$coefficients[4]*dtmp[,"ENSEMBLE"]
					Part_Res4<-res.glmmML(Model_Best,dtmp)/(fitted.glmmML(Model_Best,dtmp)*(1-fitted.glmmML(Model_Best,dtmp)))+Model_Best$coefficients[5]*dtmp[,"site_type"]

					surv_end ~ start + moist_score + ENSEMBLE + site_type
			#plot the partial residuals and the smoothed plot
					split.screen(c(2,2))
					screen(1)
					plot(dtmp[,"start"],Part_Res1,ylab="partial residuals")
					lines(seq(min(dtmp[,"start"]),max(dtmp[,"start"]),length.out=100),predict(loess(formula=Part_Res1~start,data=dtmp),newdata=data.frame(start=seq(min(dtmp[,"start"]),max(dtmp[,"start"]),length.out=100))))
					screen(2)
					plot(dtmp[,"moist_score"],Part_Res2,ylab="partial residuals")
					lines(seq(min(dtmp[,"moist_score"]),max(dtmp[,"moist_score"]),length.out=100),predict(loess(formula=Part_Res2~moist_score,data=dtmp),newdata=data.frame(moist_score=seq(min(dtmp[,"moist_score"]),max(dtmp[,"moist_score"]),length.out=100))))
					screen(3)
					plot(dtmp[,"ENSEMBLE"],Part_Res3,ylab="partial residuals")
					lines(seq(min(dtmp[,"ENSEMBLE"]),max(dtmp[,"ENSEMBLE"]),length.out=100),predict(loess(formula=Part_Res3~ENSEMBLE,data=dtmp),newdata=data.frame(ENSEMBLE=seq(min(dtmp[,"ENSEMBLE"]),max(dtmp[,"ENSEMBLE"]),length.out=100))))
					close.screen(all = TRUE)
			
			#from sjplot lmer
			sjp.glmer(SurBSP, type = "re.qq")

	
	#============================================================================================#	
	#============================================================================================#
	# VISUAL PLOTs
			SurBPSW <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	

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

							  # graph it - type 2
							ggplot(plotdat2, aes(x = start, y = PredictedProbability)) +
							  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = site_type), alpha = .20) +
							  ylim(c(0, 1)) 
							 
							 # graph it - type 3 ~ just lines
							ggplot(plotdat2, aes(x = start, y = PredictedProbability)) +
							  geom_line(aes(colour = site_type), size = 2) +
							  ylim(c(0, 1))


							  
							  
							  
							  
							  
	# VISUAL PLOTs
		sjp.glmer(SurBSP, type = "re.qq")
		
		# for plotting responses only 
		SurBPSW <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	
		sjp.glmer(SurBPSW, type = "ri.pc", facet.grid = FALSE, vars = "start")
		sjp.glmer(SurBPSW, type = "ri.pc", facet.grid = FALSE, vars = "moist_score")
		sjp.glmer(SurBPSW, type = "ri.pc", facet.grid = FALSE, vars = "ENSEMBLE")
		#sjp.glmer(SurBPSW, type = "ri.pc", facet.grid = FALSE, vars = "site_type")

# usefull extraction calls	
# from: http://jaredknowles.com/journal/2014/5/17/mixed-effects-tutorial-2-fun-with-mermod-objects
VarCorr(SurBPSW)
attr(VarCorr(SurBPSW), "sc")
slotNames(SurBPSW)
SurBPSW@call
methods(class = "merMod")
isREML(SurBPSW)

SurBP1 <- glmmML(surv_end ~ start + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBP2 <- glmmML(surv_end ~ start + moist_score, cluster=site, family=binomial, na.action=na.omit, data=d)


	





































# CHAPTER 21 FROM ZURR


##############################################################################################
# 

#============================================================================================#
# 0. load libraries & set directories to computer
#============================================================================================#
 
# set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
setwd(path.set)
source("00_SetDirectories.R") # directory script (edit for your own computer). 
# LIBRARIES
library(IPMpack)
library(lme4)
require(ggplot2)
require(lattice)

#============================================================================================#
# 1. Check data structure & do log transformation (where necessary).
#============================================================================================#

#Open 2014 plant datafile 
setwd(path.dat.raw)
plantdat <- read.csv(file="plantdata.csv") # raw data
dim(plantdat)
setwd(path.dat)
d <- read.csv(file="Data_2014.csv") # revised truncated dataframe from 01DataStruc_check"
d$start.height_ln <- log(d$start.height + 1)
colnames(d) 
site <- levels(d$site); site <- as.factor(site) # study sites
d$fec <- d$fr.fl.total
d$fec[d$fec==0]<- NA
# binary response variable for flowering probability, given survival to fall
d$pFlower[ d$surv2 == 0 ] <- NA
d$pFlower[ d$surv2 == 1 & d$fr.fl.total == 0 ] <- 0
d$pFlower[ d$surv2 == 1 & d$fr.fl.total >= 1 ] <- 1
#summary & str
summary(d)
str(d)

#============================================================================================#
# 2. Starting state variable (what should be used as start size from greenhouse?). 
#============================================================================================#

d$start <- 0.01*(d$start.height_ln) + 0.05*(d$pot) +  0.1*(d$start.height_ln*d$pot)# multiplicative

#============================================================================================#
# 3. VITAL RATE REGRESSION MODESL
#============================================================================================#

###########################################
# GROWTH
# refit model with REML and validate model 
library(nlme) 
Mgr <- lme(size2_ln ~ start, random= ~1|site/Uplot, 
				method="ML", data=d, na.action=na.omit)
summary(Mgr)

###########################################
# FECUNDITY
# Fit in site & plot as random effects.
library(nlme)
Mfr <- glmer(fec ~ size2_ln + (1|site/Uplot), 
	family = poisson, data=d, na.action=na.omit)
summary(Mfr) #

###########################################
# SURVIVAL
library(lme4)
Msu <- glmer(surv2 ~ start + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)
summary(Msu)
				
############################################
# pFLOWER
library(lme4)
Mfl <- glmer(pFlower ~ size2_ln + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)
summary(Mfl)



#============================================================================================#
# X. Attach plot level variables 
#============================================================================================#
setwd(path.dat.raw) # where environmental data is stored
enviro <- read.csv(file="environmental_variables.csv") # plot level enviornmental variables
colnames(enviro) # seems ok
levels(enviro$SITE)[levels(enviro$SITE)=="HUNT"] <- "HUNTER"
levels(d$site); levels(enviro$SITE) 
levels(d$plot); levels(enviro$PLOT) 
enviro$Un.plot <- paste(enviro$SITE, enviro$PLOT, sep="") # unique plot variable to mere with 
###### will need to merge these data frames
head(enviro$Un.plot); head(d$Uplot)

##### MERGE
d <- merge(d, enviro, by.x = "Uplot", by.y = "Un.plot", all.x = TRUE, all.y = FALSE)

#============================================================================================#
# X. Attach SITE level variables 
#============================================================================================#
setwd(path.dat.raw) # where site data is stored
SDMsite <- read.csv(file="occ_site_preds_sept2014.csv") # site level enviornmental variables
colnames(SDMsite) # seems ok
levels(d$site); levels(SDMsite$MERGE) # HUNTER needs to be renamed to match
levels(d$plot); levels(enviro$PLOT) # HUNTER needs to be renamed to match

##### MERGE
d <- merge(d, SDMsite, by.x = "site", by.y = "MERGE", all.x = TRUE, all.y = FALSE)
d$site_type <- d$site
	library(plyr)
	revalue(d$site_type, c("CALAPOOIA" = "beyond")) -> d$site_type
	revalue(d$site_type, c("WILEY" = "beyond")) -> d$site_type
	revalue(d$site_type, c("HUNTER" = "beyond")) -> d$site_type
	revalue(d$site_type, c("THOMAS" = "beyond")) -> d$site_type
	revalue(d$site_type, c("ROCK" = "occupied")) -> d$site_type
	revalue(d$site_type, c("COAST" = "occupied")) -> d$site_type
	revalue(d$site_type, c("LOOK" = "unoccupied")) -> d$site_type
	revalue(d$site_type, c("MOSBY" = "unoccupied")) -> d$site_type
levels(d$site_type)
	
#============================================================================================#
# X. Partial Dependence Plots
#============================================================================================#


###########################################
# Standardize Environmental Variables
d$start_sd <- ((d$start - mean(d$start, na.rm=TRUE))/sd(d$start, na.rm=TRUE))
d$ENSEMBLE_sd <- ((d$ENSEMBLE - mean(d$ENSEMBLE, na.rm=TRUE))/sd(d$ENSEMBLE, na.rm=TRUE))
d$moist_score_sd <- ((d$moist_score - mean(d$moist_score, na.rm=TRUE))/sd(d$moist_score, na.rm=TRUE))
###########################################

###########################################
# GROWTH
Mgr <- glm(size2_ln ~ start_sd + ENSEMBLE + moist_score + ENSEMBLE*moist_score, data=d, na.action=na.omit)

# SURVIVAL
Msu <- glm(surv2 ~ start + ENSEMBLE + moist_score,
				family=binomial(link="logit"), na.action=na.omit, data=d)

# pFlower
Mfl <- glm(pFlower ~ size2_ln + ENSEMBLE + moist_score,
				family=binomial(link="logit"), na.action=na.omit, data=d)

# Fecundity 
Mfr <- glm(fec ~ size2_ln + ENSEMBLE + moist_score,
				family=poisson, data=d, na.action=na.omit)
		
#########################################


#============================================================================================#
# 17. CHAPTER 17
#============================================================================================#

# 17.4 Estimating Common Patterns Using Additive MIXED MODEL

library(mgcv); library(nlme)
M1 <- gam(size2_ln ~ site + Uplot + s(start), data=d) # explore relationship
E <- resid(M1)
F <- fitted(M1)

op <- par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(M1)
plot(F, E, xlab = "Fitted values", ylab = "Residuals")
par(op)

########################################################################
# LIST OF POSSIBLE MODELS
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000) # iteration settings

f1 <- formula(size2_ln ~ s(start) + site)

# A - Homogeneity model added as a reference point
M17.4A <- gamm(f1, random=list(site= ~1), method="REML", control = lmc, data=d)
# B - Heterogeneity per site, but homogeneity within a site given starting size (random intercept) 
M17.4B <- gamm(f1, random=list(site= ~1), method="REML", control = lmc, data=d,
			weights = varIdent(form =~ 1 | site))
# C - Homogeneity across site, but heterogeneity across sizes for each site (random slope?)
M17.4C <- gamm(f1, random=list(site= ~1), method="REML", control = lmc, data=d,
			weights = varPower(form =~ start))
# D - Heterogeneity across sites and heterogeneity across start sizes
M17.4D <- gamm(f1, random=list(site= ~1), method="REML", control = lmc, na.action=na.omit, data=d,
			weights = varComb(varIdent(form =~ 1|site), varPower(form =~start)))
# E - Heterogeneity across sites and site*start
M17.4E <- gamm(f1, random=list(site= ~1), method="REML", control = lmc, data=d,
			weights = varComb(varIdent(form =~ 1|site), 
			varPower(form =~start|site)))

# random slope for site
	AIC(M17.4A$lme, M17.4B$lme, M17.4C$lme)			

################################################
# error throughout rest Zurr script move to chapter 19 & 20

#============================================================================================#
# 20. CHAPTER 20
#============================================================================================#

# is random effect of plot significant
library(nlme)
f1 <- formula(size2_ln ~ start*site)
M1 <- gls(f1, method = "REML", data = d, na.action=na.omit)

M2 <- lme(f1, random =~1 | site / Uplot,
data = d, method = "REML", na.action=na.omit)
anova(M1, M2)
# Plot random effect is significant

# add additonal variance based on site_type and see if improves model 
M3 <- lme(f1, random =~1 | site / Uplot,
data = d, method = "REML", na.action=na.omit,
weights = varIdent(form =~ 1 | site_type))
anova(M2, M3)

################################################


#============================================================================================#
# 21. CHAPTER 21
#============================================================================================#

# Final Growth, Survivorship, Flower and Fecundity Model
Gr <- lme(size2_ln ~ start + moist_score + ENSEMBLE, random= ~1|site/Uplot, 
				method="ML", data=d, na.action=na.omit)
				
Fec <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot),
				family = poisson, data=d, na.action=na.omit)

Flr <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
				
Sur <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
##############
AIC(Sur)


# COMPARE MODELS BASED ON AIC (ZURR CHAP 21) 
# 	Method outlined in Zurr works, but does not allow for nested random effects.
#	Therefor should compare AICs with lmer too to see if we get same results (below)

# standardize variables to improve convergence of fitting algorithm 
	# also same scale, effect size ect. 
# subset
d2 <- d[ ,c("site", "Uplot", "pFlower", "surv_end", "fec", "size2_ln", "start", "moist_score", "site_type", "ENSEMBLE")]


############
# run competing models
#B-basic; P - Plot level; S- SDM; W - within/beyond
library(qpcR)

# SURVIVAL 
SurB <- glmmML(surv_end ~ start, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBP <- glmmML(surv_end ~ start + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBS <- glmmML(surv_end ~ start + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBW <- glmmML(surv_end ~ start + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBWP <- glmmML(surv_end ~ start + moist_score + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBWS <- glmmML(surv_end ~ start + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBSP <- glmmML(surv_end ~ start + moist_score + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBPSW <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	

names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW")
my_aic <- c(SurB$aic, SurBP$aic, SurBS$aic, SurBW$aic, SurBWP$aic, SurBWS$aic, SurBSP$aic, SurBPSW$aic)
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

###########################################################################################################
# pFlower
FlrB <- glmmML(pFlower ~ size2_ln, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBP <- glmmML(pFlower ~ size2_ln + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBS <- glmmML(pFlower ~ size2_ln + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBW <- glmmML(pFlower ~ size2_ln + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBWP <- glmmML(pFlower ~ size2_ln + moist_score + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBWS <- glmmML(pFlower ~ size2_ln + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBSP <- glmmML(pFlower ~ size2_ln + moist_score + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBPSW <- glmmML(pFlower ~ size2_ln + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	

names <- c("FlrB", "FlrBP", "FlrBS", "FlrBW", "FlrBWP", "FlrBWS", "FlrBSP", "FlrBPSW")
my_aic <- c(FlrB$aic, FlrBP$aic, FlrBS$aic, FlrBW$aic, FlrBWP$aic, FlrBWS$aic, FlrBSP$aic, FlrBPSW$aic)
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

###########################################################################################################
# FECUNDITy
FecB <- glmmML(fec ~ size2_ln, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBP <- glmmML(fec ~ size2_ln + moist_score, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBS <- glmmML(fec ~ size2_ln + ENSEMBLE, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBW <- glmmML(fec ~ size2_ln + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBWP <- glmmML(fec ~ size2_ln + moist_score + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBWS <- glmmML(fec ~ size2_ln + ENSEMBLE + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBSP <- glmmML(fec ~ size2_ln + moist_score + ENSEMBLE, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBPSW <- glmmML(fec ~ size2_ln + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)	

names <- c("FecB", "FecBP", "FecBS", "FecBW", "FecBWP", "FecBWS", "FecBSP", "FecBPSW")
my_aic <- c(FecB$aic, FecBP$aic, FecBS$aic, FecBW$aic, FecBWP$aic, FecBWS$aic, FecBSP$aic, FecBPSW$aic)
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

###########################################################################################################
# growth
GrB <- lmer(size2_ln ~ start + (1|Uplot), na.action=na.omit, REML=FALSE, data=d)
GrBP <- lmer(size2_ln ~ start + moist_score + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBS <- lmer(size2_ln ~ start + ENSEMBLE + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBW <- lmer(size2_ln ~ start + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWP <- lmer(size2_ln ~ start + moist_score + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWS <- lmer(size2_ln ~ start + ENSEMBLE + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBSP <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBPSW <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)	

names <- c("GrB", "GrBP", "GrBS", "GrBW", "GrBWP", "GrBWS", "GrBSP", "GrBPSW")
my_aic <- c(AIC(GrB), AIC(GrBP), AIC(GrBS), AIC(GrBW), AIC(GrBWP), AIC(GrBWS), AIC(GrBSP), AIC(GrBPSW))
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#============================================================================================#
# ALTERNATIVE WITH LMER

# SURVIVAL 
SurB <- glmer(surv_end ~ start + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBP <- glmer(surv_end ~ start + moist_score + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBS <- glmer(surv_end ~ start + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBW <- glmer(surv_end ~ start + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBWP <- glmer(surv_end ~ start + moist_score + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBWS <- glmer(surv_end ~ start + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBSP <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBPSW <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	

names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW")
my_aic <- c(AIC(SurB), AIC(SurBP), AIC(SurBS), AIC(SurBW), AIC(SurBWP), AIC(SurBWS), AIC(SurBSP), AIC(SurBPSW))
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

###########################################################################################################
# pFLOWER 
FlrB <- glmer(pFlower ~ size2_ln + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBP <- glmer(pFlower ~ size2_ln + moist_score + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBS <- glmer(pFlower ~ size2_ln + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBW <- glmer(pFlower ~ size2_ln + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBWP <- glmer(pFlower ~ size2_ln + moist_score + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBWS <- glmer(pFlower ~ size2_ln + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBSP <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBPSW <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	

names <- c("FlrB", "FlrBP", "FlrBS", "FlrBW", "FlrBWP", "FlrBWS", "FlrBSP", "FlrBPSW")
my_aic <- c(AIC(FlrB), AIC(FlrBP), AIC(FlrBS), AIC(FlrBW), AIC(FlrBWP), AIC(FlrBWS), AIC(FlrBSP), AIC(FlrBPSW))
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')
						
###########################################################################################################
# Fecundity 
FecB <- glmer(fec ~ size2_ln + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBP <- glmer(fec ~ size2_ln + moist_score + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBS <- glmer(fec ~ size2_ln + ENSEMBLE + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBW <- glmer(fec ~ size2_ln + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBWP <- glmer(fec ~ size2_ln + moist_score + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBWS <- glmer(fec ~ size2_ln + ENSEMBLE + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBSP <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBPSW <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)	

names <- c("FecB", "FecBP", "FecBS", "FecBW", "FecBWP", "FecBWS", "FecBSP", "FecBPSW")
my_aic <- c(AIC(FecB), AIC(FecBP), AIC(FecBS), AIC(FecBW), AIC(FecBWP), AIC(FecBWS), AIC(FecBSP), AIC(FecBPSW))
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')
						
				
###########################################################################################################
# growth
GrB <- lmer(size2_ln ~ start + (1|site/Uplot), na.action=na.omit, REML=FALSE, data=d)
GrBP <- lmer(size2_ln ~ start + moist_score + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBS <- lmer(size2_ln ~ start + ENSEMBLE + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBW <- lmer(size2_ln ~ start + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWP <- lmer(size2_ln ~ start + moist_score + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWS <- lmer(size2_ln ~ start + ENSEMBLE + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBSP <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBPSW <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)	

names <- c("GrB", "GrBP", "GrBS", "GrBW", "GrBWP", "GrBWS", "GrBSP", "GrBPSW")
my_aic <- c(AIC(GrB), AIC(GrBP), AIC(GrBS), AIC(GrBW), AIC(GrBWP), AIC(GrBWS), AIC(GrBSP), AIC(GrBPSW))
library(qpcR)
moo <- akaike.weights(my_aic)
tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

#########			
				










#============================================================================================#
# PLOTS
#============================================================================================#


# Final Growth, Survivorship, Flower and Fecundity Model
Gr <- lme(size2_ln ~ start + moist_score + ENSEMBLE, random= ~1|site/Uplot, 
				method="ML", data=d, na.action=na.omit)
				
Fec <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot),
				family = poisson, data=d, na.action=na.omit)

Flr <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
				
Sur <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
###################################################################

fixef(Sur)
#confint(Sur, level = 0.99)
moo <- ranef(Sur); moo <- moo$Uplot



# contrast over conditional
###############################
library(visreg)
par(mfrow=c(1,2))
visreg(Gr, "ENSEMBLE", by="moist_score", breaks=4, overlay=TRUE,
				partial=TRUE, jitter=TRUE, legend=FALSE, 
				points=list(pch = c(0,6,12,1,7,13,2,8,14,3,9,4,10,5,11,0)[as.numeric(as.factor(d$plot))], col="darkgrey", cex=0.5),
				ylab="Relative Growth Rate", xlab="Predicted Site Suitability From ENM",
				ylim=c(3, 4.2))
				
	
###############		
# tricky to plot glmer in visreg so follow from
# http://www.strengejacke.de/sjPlot/sjp.glmer/
library(sjPlot)
sjp.glmer(Sur)
sjp.glmer(Sur, type = "fe", sort.coef = TRUE)
sjp.glmer(Flr, type = "fe", sort.coef = TRUE)
sjp.glmer(Fec, type = "fe", sort.coef = TRUE)




methods(class="merMod")
Sur <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
###########################################
sjp.lmer(Sur)
moo <- ranef(Sur); moo <- moo$Uplot


###########################################
# JARED KNOWLES - tutorials 
#Fun with merMod Objects
#http://jaredknowles.com/journal/2014/5/17/mixed-effects-tutorial-2-fun-with-mermod-objects




































				
Sur <- glmer(surv_end ~ start + moist_score + (1|site), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)

sjp.glmer(Sur, type="ri.pc",  show.se= TRUE)




