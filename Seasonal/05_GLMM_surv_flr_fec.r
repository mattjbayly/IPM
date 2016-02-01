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

			names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW")
			my_aic <- c(SurB$aic, SurBP$aic, SurBS$aic, SurBW$aic, SurBWP$aic, SurBWS$aic, SurBSP$aic, SurBPSW$aic)
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

			names <- c("SurB", "SurBP", "SurBS", "SurBW", "SurBWP", "SurBWS", "SurBSP", "SurBPSW")
			my_aic <- c(AIC(SurB), AIC(SurBP), AIC(SurBS), AIC(SurBW), AIC(SurBWP), AIC(SurBWS), AIC(SurBSP), AIC(SurBPSW))
			library(qpcR)
			moo <- akaike.weights(my_aic)
			tabby <- data.frame(matrix(my_aic,ncol=1,byrow=TRUE))
			tabby <- cbind(names, tabby, moo[3]); colnames(tabby) <- c( 'MODEL', 'AIC', 'WEIGHT')

#============================================================================================#
# SURVIVORHSIP - EXPLORE MODEL ADEQUECY (ok?)
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

