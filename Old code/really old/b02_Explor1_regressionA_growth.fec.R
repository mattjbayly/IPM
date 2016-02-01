#title: "VITAL RATE REGRESSION FUNCTIONS"
#author: "Matthew Bayly"
#modified: "Thursday, January 16, 2015"

# PURPOSE OF SCRIPT: 
# Develop the functional form of each vital rate regression function,
# for each of the four vital rate regression equations (growth, survival
# pflower & fecundity. Run global models here to justify the functional 
# form of models for each vital rate with size. Ie. Should equations be 
# a simple linear form, quadratic or just an intercept (flat line). 
# Continuous variables: growth & fecundity are in script - 02_Vital_rate_regressionA
# Binary variables: surv & pflower are in script - 03_Vital_rate_regressionB

# To justify the use of a more complex model form (ie. vital rate has a quadratic 
# relationship with size), compare models with AIC & visually by subsampling 
# the data. Separate models will be run later on in the analysis for each site,
# but this script is just to justify the size based relationships. 

# All code used is this analysis was copied and modified from:
# 1. CJE Metcalf, SM McMahon, R Salguero-Gomez, E Jongejans and Cory Merow
		#(2013). IPMpack: Builds and analyses Integral Projection Models
		#(IPMs).. R package version 2.0.
		# http://CRAN.R-project.org/package=IPMpack
# 2. Ellner, S. P. & Rees, M. Integral projection models for species
		# with complex demography. Am. Nat. 167, 410–428 (2006)
		# SUPPLIMENTARY MATERIAL FILES
		
# INDEX: 
	# 0. load libraries & set directories to computer
	# 1. Check data structure & do log transformation (where necessary).
	# 2. Starting state variable (what should be used as start size). 
			# a. from pot size & total height. 
	# 3. Explore global models based on AIC. 
		# 3.1 Growth model (form & AIC)
		# 3.2 Fecundity model (form & AIC)
	# 4. Sensitivity of vital rates to model structure and sample size:
		# 4.1 Growth model (Merow et al 2014: 4.5.2).
		# 4.2 Growth 'increment' model (Merow et al 2014: 4.5.3). #(temporary exploration)
		# 4.3 Fecundity model (Merow et al 2014: 4.5.2).

# RANDOM EFFECTS:
# Individual plants are nested within plots, which themselves are 
# nested within sites. Need to account for nested structure in the
# data when testing relationships using random effects models.  

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
# PLOT GROWTH TO FALL ACROSS SITES
par(mfrow=c(2,4),mar=c(4,4,2,1))
### May to July 
for (i in 1:length(site)){
   current <- d[d$site==site[i], ]
  plot(current$start.height_ln, current$size1_ln, xlab="Green house_ln", ylab="July size_ln", 
          col=as.factor(current$source), main=site[i], pch=19, xlim=c(0,6), ylim=c(0,6), 
          cex=1.2); abline(0, 1)
  }
### July - Sept
for (i in 1:length(site)){
   current <- d[d$site==site[i], ]
  plot(current$size1_ln, current$size2_ln, xlab="Spring size", ylab="Fall size", 
          col=as.factor(current$source), main=site[i], pch=19, xlim=c(0,6), ylim=c(0,6), 
          cex=1.2); abline(0, 1)
  }
  
#summary & str
summary(d)
str(d)
# data already log transformed from previous r-script
# run a quick correlation plot 
library(PerformanceAnalytics)
chart.Correlation(d[, c("start.height", "start.height_ln", "size1",
	"size2", "pot", "size1_ln", "size2_ln", "fr.fl.total")])

	
	
	
#============================================================================================#
# 2. Starting state variable (what should be used as start size from greenhouse?). 
#============================================================================================#

# Weight start height in greenhouse (prior to pruning & pot size) to combine 
# into a single variable. 
explorG = lm(size1_ln ~ start.height_ln + pot + pot*start.height_ln, data=d)
# fitted values of explorS on exponential scale (
coef(explorG)

# Play around with variable combinations 
# Normalize start height

d$startFULLfit <- 1.75927313 + 0.09923*d$start.height_ln + -0.072*d$pot +  0.093*d$start.height_ln*d$pot# multiplicative
d$startFULLfit2 <- log((0.01*d$start.height_ln + 0.05*d$pot +  0.1*d$start.height_ln*d$pot)+1)# multiplicative

d$start1 <- log((0.2*d$start.height_ln + 0.3*d$pot) + 1)
chart.Correlation(d[, c("startFULLfit2",
	"startFULLfit", "pot", "start.height_ln", "size1_ln", "size2_ln", "fr.fl.total")])

# Which combinations work best for survivroship
explorS1 = glm(surv1 ~ d$startFULLfit,  family = binomial(link="logit"), data=d)
explorS2 = glm(surv1 ~ d$startFULLfit2,  family = binomial(link="logit"), data=d) # winning model
anova(explorS1, explorS2)

# Which combination works best for growth
explorG1 = lm(size2_ln ~ startFULLfit, data=d)
explorG2 = lm(size2_ln ~ startFULLfit2, data=d)
explorG3 = lm(size2_ln ~ pot, data=d)
explorG4 = lm(size2_ln ~ start.height_ln, data=d)

summary(explorG1)$adj.r.squared
summary(explorG2)$adj.r.squared
summary(explorG3)$adj.r.squared
summary(explorG4)$adj.r.squared

############################################################################
# start potential seems to work best using an approximation of start height &
# pot size. START = 0.01*startheight_ln + 0.05*pot + 0.1height*pot

d$start <- 0.01*(d$start.height_ln) + 0.05*(d$pot) +  0.1*(d$start.height_ln*d$pot)# multiplicative


#============================================================================================#
# 3. Explore global models based on AIC. 
#============================================================================================#
library(lme4) # main package to use

# add possible quadratic terms to test
d$start2 <- (d$start)^2 # quadratic term
d$start3 <- (d$start)^3 # polynomial term

# Top-down strategy for model selection in mixed effect structure. Start with as many
# variables as possible in fixed component & compare random effect structre. 

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

# 3.1 Growth model, given survivorship (form & AIC)
library(nlme)
# no random effect
		growth_gls <- gls(size1_ln ~ 1 + start + start2 + start3, 
				method="REML", data=d, na.action = na.omit)
# random intercept model, with nested plots within sites
		growth_lme1 <- lme(size1_ln ~ 1 + start + start2 + start3,
				random = ~1|site/Uplot, data=d, method="REML", na.action = na.exclude)
# random intercept & slope model, with nested plots within sites
		growth_lme2 <- lme(size1_ln ~ 1 + start + start2 + start3,
				random = ~1 + start|site/Uplot, data=d, method="REML", na.action = na.exclude)
AIC(growth_gls, growth_lme1, growth_lme2)
# from AIC scores, random int&slope seems to out preform others, but wont have enough 
# data to work with it. Also variance explained by random slope is marginal relative 
# to random intercept model.
summary(growth_lme1)
summary(growth_lme2)
# plot to visualize residual spread
par(mfrow=c(2,4))
plot((lm(size1_ln ~ 1 + start, data=d)), select= c(1))
hist(d$start); hist(d$size1_ln); hist(d$start.height_ln)
# data a little bit hetero - variance much larger with smaller start sizes
# probably because start data is unbalance across size and many more observation 
# across small size classes.
#######################################################
 
# COMPARE RANDOM INTERCEPT MODEL FOR GROWTH
# Fit the model with GLS
library(nlme)
Form <- formula(size1_ln ~ + start + start2 + start3)
M.gls <- gls(Form, data = d, na.action=na.omit)
# need to add random intercept prior to choosing the variance structure
M1.lme <- lme(Form, random= ~1|site/Uplot, 
		method="REML", data=d, na.action=na.omit) 
anova(M1.lme, M.gls)
# random intercept is significant for site and/or slope
# AIC of rand intercept is also smaller.  
# L = 365 (df = 7, p <0.0001)
# plot out model, is everything OK?
E2 <- resid(M1.lme, type="normalized")
F2 <- fitted(M1.lme)
op <- par(mfrow=c(2,2), mar=c(4,4,3,2))
MyYlab <- "Residuals"
plot(x=F2, y=E2, xlab="Fitted Values", ylab=MyYlab)
# some violation of homogeneity in residuals. 
# ASSES OPTIMAL FIXED STRUCTURE
summary(M1.lme)
# neither the quadratic or polynomial terms are significant so drop 
# least sig term & reapply the model, changing method to 'ML' rather than 'REML' 
anova(M1.lme)
# eliminate non-sig fixed terms with the likelihood ratio test 
M1.full <- lme(Form, random = ~ 1|site/Uplot, 
	method="ML", data=d, na.action=na.omit)
M1.A <- update(M1.full, .~. -start2)
M1.B <- update(M1.full, .~. -start3)
anova(M1.full, M1.A) # drop this after start3
anova(M1.full, M1.B) # drop polynomial term for sure!
###########################################
Form2 <- formula(size1_ln ~ start + start2)
M2.full <- lme(Form2, random=~1|site/Uplot, 
	method="ML", data=d, na.action=na.omit)
M2.A <- update(M2.full, .~. -start2)
M2.B <- update(M2.full, .~. -start)
anova(M2.full, M2.A) # drop quadratic term p=0.09, so maybe marginally significant!?
anova(M2.full, M2.B) # keep start (as expected)
###########################################
Form3 <- formula(size1_ln ~ start)
M3.full <- lme(Form3, random=~1|site/Uplot, 
		method="ML", data=d, na.action=na.omit)
M3.A <- update(M3.full, .~. -start)
anova(M3.A) # keep 'start'
# results are the same when size1_ln is substituted with size2_ln
# end of basic form of model selection processes for
# fixed effects. Quadratic term (size2) is moderately sig, so
# worried about completely dropping it at this point. With revisit 
# below in the Merow variable form IPM tutorial to be certain relationship
# is indeed linear. 
###########################################
# refit model with REML and validate model 
M5 <- lme(size1_ln ~ start, random= ~1|site/Uplot, 
				method="ML", data=d, na.action=na.omit)
summary(M5)
# random intercept for site is normally distributed with 
	# mean of 0 & SD of 0.375 (variance 0.375^2)
# random intercept for plot (within site) is norm dist with 
	# mean of 0 & SD of 0.568
# the residual term is also normally distributed with mean of 
	# zero & SD of 0.643, these 
# correlation of observation across the same site
0.3750124^2/(0.5685559^2 + 0.6435375^2 + 0.3750124^2) # ==0.16 (low, but significant)
# correlation of observation across the same plot within a site
0.5685559^2/(0.5685559^2 + 0.6435375^2 + 0.3750124^2) # ==0.36
###########################################
# FINAL CHECKS: LOESS, plots & visualizations
library(lattice)
E2 <- resid(M5, type="normalized")
xyplot(E2 ~ start|site, data=d, ylab="Residuals",
		xlab="Start", 
		panel = function(x,y){
		panel.grid(h=-1, v=2)
		panel.points(x,y, col=1, pch=19, cex=0.1)
		panel.loess(x,y, span=0.5, col="red", lwd=1)})
# with(d, table(d$plot, d$site)) # reminder of plants within sites
library(mgcv)
# quick addative mixed model 
M6 <- gamm(size1_ln ~ s(start), random=list(Uplot=~1),  
			method="ML", data=d, na.action=na.omit)~1+ bio14 	+s(bio14,2) +s(bio14,3)	+s(bio14,4),
# quick addative mixed model 
M6 <- gamm(size1_ln ~ s(start), random=list(Uplot=~1),  
			method="ML", data=d, na.action=na.omit)
#summary(M6)
#summary(M6$gam) # output of smoothers and parameters in model
#anova(M6$gam) # output (warning not sequential testing!)
par(mfrow=c(2,2))
plot(M6$gam)
plot(M6$lme)
dev.off()
# conclusion - we can be pretty happy with the linear term for size
	
# END: Growth model, given survivorship

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


# 3.2 Fecundity model (form & AIC)
# model fruit & flower produced given size in fall
par(mfrow=c(1,1)); dev.off()
plot(d$size2_ln, d$fr.fl.total, pch=19, cex=0.6, 
		col=as.factor(d$site), xlab="Fall Size", ylab="Total FR & FL")
# very zero inflated, most plants don't produce regardless of size. 
# pflower (another vital rate, will model probability of flowering so for now
# we can assume plants have already flowered.
d$fec <- d$fr.fl.total
d$fec[d$fec==0]<- NA
# table per site
aggregate(fec ~ site, data=d, FUN=sum)
# table fruits in plots within sites
aggregate(fr.fl.total ~ plot + site, data=d, FUN=sum)

#######################################
# sketch out basic Poisson model 
d2 <- d[ ,c("fec", "size2_ln", "site")]
dtemp <- d2[complete.cases(d2),]
M1 <- glm(fec ~ size2_ln, family = poisson, data = dtemp)

require(graphics)
d2 <- d[ ,c("fec", "size2_ln", "site")]
dtemp <- d2[complete.cases(d2),]
# Sketch out fitted values, optional 
plot(dtemp$size2_ln, dtemp$fec, pch=19, cex=0.6, col=as.factor(dtemp$site), xlab="Fall Size", ylab="Total FR & FL")
MyData <- data.frame(size2_ln = seq(from = 0, to = 6,
		 by = 0.01))
G <- predict(M1, newdata = MyData, se=TRUE, type = "link", na.action=na.omit)
F <- exp(G$fit)
FSEUP <- exp(G$fit + 1.96 * G$se.fit)
FSELOW <- exp(G$fit - 1.96 * G$se.fit)
lines(MyData$size2_ln, F, lty = 1)
lines(MyData$size2_ln, FSEUP, lty = 2)
lines(MyData$size2_ln, FSELOW, lty = 2)

# plot quick line for each site
for(i in 1:length(site)){
		current <- dtemp[which(dtemp$site==site[i]), ]
		if(length(current$fec) > 4) {
				M1 <- glm(fec ~ size2_ln, family = poisson, data = current)
				G <- predict(M1, newdata = MyData, se=TRUE, type = "link", na.action=na.omit)
				F <- exp(G$fit)
				lines(MyData$size2_ln, F, lty = 1, col="red")
				} else {(print("too few"))
				}
	}
# relationship will not be linear, probably exponential, unless the fruit
# count data is transformed with to a negative binomial distribution 
#######################################################################

# basic poisson model 
d2 <- d[ ,c("fec", "size2_ln", "site")]
dtemp <- d2[complete.cases(d2),]
M1 <- glm(fec ~ size2_ln, family = poisson, data = dtemp)
summary(M1) # pseduo R^2 = 41%
# pseudo R^2 in GLM:  ([null deviance - residual deviance] / null deviance), with null deviance 
# corresponding to the worst possible model. 
((M1$"null.deviance" - deviance(M1))/M1$"null.deviance")*100 # Pseudo R^2

# full model with everything (just exploratory - not accounting for independence)
M2<-glm(fec ~ size2_ln + site, family = poisson, data=d, na.action=na.omit) # + source dropped 
summary(M2)
#Dispersion parameter for poisson family taken to be 1 - quick check 
# [residual deviance / (n - p)] ~ should be around 1, will not work in this case
drop1(M2,test = "Chi") # source is not significant, site is sig. 


# Used Poisson distribution (above), detected overdispersion 
# & corrected the standard errors with a quasi-GLM model (below). 
M4<-glm(fec ~ size2_ln + site, family = quasipoisson, data=d, na.action=na.omit) # + source dropped 
summary(M4)
# Dispersion parameter for quasi-poisson family taken to be 5.648178 (rather than 1 in Poisson dist).
drop1(M4,test = "F") # site differences still holds,
#######################################################################

# Fit in site & plot as random effects.
library(MASS); #library(AED)
library(nlme)
M5 <- glmer(fec ~ size2_ln + site + (1|site/Uplot), 
	family = poisson, data=d, na.action=na.omit)
summary(M5) # corr between int & slope for 
# try with interaction between site and slope
M6 <- glmer(fec ~ size2_ln + site + site*size2_ln + (1|site/Uplot),
	family = poisson, data=d, na.action=na.omit)
anova(M5, M6) # just for fun

#######################################################################
# Fit gam smoother to data to visualize
library(mgcv)
O4.gamm <- gamm(fec ~ s(size2_ln) + source, random= list(Uplot =~1), data=d,
		na.action=na.omit, family=poisson)
summary(O4.gamm, cor=FALSE)
anova(O4.gamm$gam)
#######################################################################	

# END: Fecundity model, given size

##############################################################################################
##############################################################################################
##############################################################################################

#============================================================================================#
# 4. Sensitivity of vital rates to model structure and sample size:
#============================================================================================#

##############################################################################################
##############################################################################################
##############################################################################################

# 4.1 Growth model (Merow et al 2014: 4.5.2).
library(IPMpack)
# set up data according to IPMpack - from GreenHouse to JULY CENSUS
d_back <- d # to get retrieve old variables

startTime <- "Greenhouse"
finishTime <- "Sept"
##################################
d$size <- d$start
d$sizeNext <- d$size2_ln 
d$stage <- "continuous"
d$stageNext <- "continuous"
d$surv <- d$surv_end
d$fec1Bolt[ d$surv2 == 0 ] <- NA
d$fec1Bolt[ d$surv2 == 1 & d$fr.fl.total == 0 ] <- 0
d$fec1Bolt[ d$surv2 == 1 & d$fr.fl.total >= 1 ] <- 1
d$fec3Head <- d$fec
#################################


# make temporary mesh points
x <- seq(from=0, to=10, length=1001) # shorten in plot space
x0 <- data.frame(size=x, size2=x*x)
minSize <- min(d$size, na.rm=T)
maxSize <- max(d$sizeNext, na.rm=T)
#################################
# temporary survival model comp.
survModelComp(dataf = d, makePlot=TRUE, legend="bottomright", mainTitle="Survival")
# as expected, global model with just size is best fit for survival 
# does it also hold true for each site? 
dev.off()
par(mfrow=c(2,4))
for(i in 1:length(site)){
	current <- d[d$site==site[i], ]
	survModelComp(dataf = current, makePlot=TRUE, mainTitle=paste(site[i], " Surv", sep=""))
	s_temp <- makeSurvObj(current, surv~size)
	assign(paste("so_",site[i], sep=""), s_temp)
	} 
# surv~size, wins for every site! ~ so use, assume survival probability is linear with size & not 
# quadratic, will investigate this more rigorously in next script. 
so <- makeSurvObj(d, surv~size)
#################################
# temporary growth model comp.
dev.off(); par(mfrow=c(1,1))
growthModelComp(dataf = d, makePlot=TRUE, legend="bottomright", mainTitle="Growth")
growthModelComp(dataf = d, makePlot=TRUE, 
		legend="bottomright", mainTitle="Growth", regressionType = "changingVar")
dev.off()
# test individually for each site
par(mfrow=c(2,4))
for(i in 1:length(site)){
	current <- d[d$site==site[i], ]
	growthModelComp(dataf = current, 
		makePlot=TRUE, 
		mainTitle=paste(site[i], "Growth", sep=""))
		g_temp <- makeGrowthObj(current, sizeNext~size)
		assign(paste("go_",site[i], sep=""), g_temp)
} 
# regardless of constant or changing Var a form of 
# sizeNext ~ size, seems to be the top model over size + size2 + size3 ect. 
# in only one case at Thomas Creek the the sizeNext = size + size2 win as the 
# top model and only by a very marginal AIC difference.
# perhaps size + size2 would have been possible with more data and much lager
# sizes that were not possible with the transplant study?
go <- makeGrowthObj(d, sizeNext ~ size)
# make P matrix
Pmatrix <- makeIPMPmatrix(survObj=so, 
			growObj=go, 
			minSize=minSize,
			maxSize=maxSize)
# plot with image plot
require(fields)
dev.off(); par(mfrow=c(1,1))
image.plot(Pmatrix@meshpoints, 
			Pmatrix@meshpoints,
			t(Pmatrix), 
			main="Pmatrix: survival and growth",
			xlab= paste("Size: ", startTime, sep=""),
			ylab = paste("Size: ", finishTime, sep=""))
			abline(0,1,lty=2, lwd=2)
#####################################
# plot individual matrix for each site
dev.off(); par(mfrow=c(2,4))
for(i in 1:length(site)){
	current <- d[d$site==site[i], ]
	so_temp <- get(paste("so_", site[i], sep=""))
	go_temp <- get(paste("go_", site[i], sep=""))
	PmatrixT <- makeIPMPmatrix(survObj=so_temp, 
			growObj=go_temp, 
			minSize=minSize,
			maxSize=maxSize)
	image.plot(PmatrixT@meshpoints, 
			PmatrixT@meshpoints,
			t(PmatrixT), 
			main=site[i],
			xlab= paste("Size: ", startTime, sep=""),
			ylab = paste("Size: ", finishTime, sep=""))
			abline(0,1,lty=2, lwd=2)
} 	

# Visualization of the P-matrix: Individuals on the line neither grow or shrink (stasis)
# Individuals above the line grow & individuals below the line shrink (over the time period) 
# color shading denotes probability mass (not number of individuals)

# diagnose
#diagnosticsPmatrix(Pmatrix, growObj=go, survObj=so, correction="constant")

# END: Growth model sensitivity 
##############################################################################################
##############################################################################################
##############################################################################################

# 4.2 Growth 'increment' model (Merow et al 2014: 4.5.3). #(temporary exploration)

# END: Growth 'increment' model sensitivity to sample size
##############################################################################################
##############################################################################################
##############################################################################################

# 4.3 Fecundity model (Merow et al 2014: 4.5.2).
# set up time period
startTime <- "Greenhouse"
finishTime <- "July"
##################################
d$size <- d$start
d$sizeNext <- d$size2_ln 
d$stage <- "continuous"
d$stageNext <- "continuous"
d$surv <- d$surv_end
d$fec1Bolt[ d$surv2 == 0 ] <- NA
d$fec1Bolt[ d$surv2 == 1 & d$fr.fl.total == 0 ] <- 0
d$fec1Bolt[ d$surv2 == 1 & d$fr.fl.total >= 1 ] <- 1
d$fec3Head <- d$fec
#################################
# make fecundity kernel 
# compare pflower model
#fo1 <- makeFecObj(d, Formula=fec3Head~1, Family="binomial", Transform="none") # int only
#fo2 <- makeFecObj(d, Formula=fec1Bolt~size, Family="binomial") # linear

# END: Fecundity model sensitivity to sample size
