#title: "VITAL RATE REGRESSION FUNCTIONS"
#author: "Matthew Bayly"
#modified: "Thursday, January 16, 2015"

# PART B - BINARY RESPONSE VARIABLES 
	# PROBABILITY OF Survival
	# PROBABILITY OF FLOWERING
	

# INDEX: 

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

# SURVIVAL MODELS 
# fit inital model assuming all independant observations
Form <- formula(surv1 ~ start + start2 + start3 + pot)
M1 <- glm(Form, data=d,  family = binomial(link="logit"), na.action = na.exclude)
drop1(M1, test="Chi")			
# intercept fits best?				

M2 <- glm(surv1 ~ start, 
			data=d,  family = binomial(link="logit"), na.action = na.exclude)
# visualize
plot(jitter(d$surv1, amount=0.05) ~ d$start, 
	pch=19, cex=0.5, col=as.factor(d$site), xlab="Start potential", ylab="Survival")

###########################################################
###########################################################

# 3.2 RANDOM INTERCEPT MODEL 
# Fit model used random effects
library(MASS)
M3 <- glmmPQL(surv1 ~ start + start2 + start3, 
				random = ~ 1|site/Uplot, family=binomial(link="logit"),
				na.action=na.omit, data=d)
summary(M3)
# the standard error of the site & plot are high in comparison 
# to the standard error of the residual, this is because these are not Pearson residuals. 
# a little bit sceptical about significance of higher order terms. 		
s <- summary(M3)
# str(s) # extract values from model 	
			
# alternative method using LMER	
library(lmer)
M4 <- glmer(surv1 ~ start + start2 + start3 + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)
###############################################################
M4.a <- glmer(surv1 ~ start + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
				
				
# Try with glmmML, can't seem to get a nested structure to work
library(glmmML)
M5 <- glmmML(surv1 ~ start + start2 + start3,
			cluster = site, family=binomial(link="logit"),
			data=d)
summary(M5)


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

# pFLOWR MODELS 

# binary response variable for flowering probability, given survival to fall
d$pFlower[ d$surv2 == 0 ] <- NA
d$pFlower[ d$surv2 == 1 & d$fr.fl.total == 0 ] <- 0
d$pFlower[ d$surv2 == 1 & d$fr.fl.total >= 1 ] <- 1

# fit initial exploratory model assuming all independent observations
Form <- formula(pFlower ~ start + start2 + start3 + pot)
M1 <- glm(Form, data=d,  family = binomial(link="logit"), na.action = na.exclude)
drop1(M1, test="Chi")			
# intercept fits best?				

M2 <- glm(pFlower ~ start, 
			data=d,  family = binomial(link="logit"), na.action = na.exclude)
# visualize
plot(jitter(d$pFlower, amount=0.05) ~ d$start, 
	pch=19, cex=0.5, col=as.factor(d$site), xlab="Start potential", ylab="Survival")

###########################################################
#  RANDOM INTERCEPT MODEL FOR PFLOWER
# Fit model used random effects
library(MASS)
M3 <- glmmPQL(pFlower ~ start + start2 + start3, 
				random = ~ 1|site/Uplot, family=binomial(link="logit"),
				na.action=na.omit, data=d)
summary(M3)
##########################################################		
s <- summary(M3)
# str(s) # extract values from model 	
			
# alternative method using LMER	
library(lme4)
M4 <- glmer(pFlower ~ start + start2 + start3 + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)
###############################################################
M4.a <- glmer(pFlower ~ start + (1|site/Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
				
				
# Try with glmmML, can't seem to get a nested structure to work
library(glmmML)
M5 <- glmmML(pFlower ~ start + start2 + start3,
			cluster = site, family=binomial(link="logit"),
			data=d)
summary(M5)


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


#============================================================================================#
# X. Fitting in environmental covariates  
#============================================================================================#


# Compare models with environmental variables 
library(glmmML)
M6 <- glmmML(surv1 ~ start + pot, cluster = site, data =d,
			na.action=na.omit, family = binomial) 
# cluster argument is the high level groups, but current not for plots. 
# direction from Section 21.4 pg. 481 (Zurr et al 
# "Mixed Effects Models and Extensions in Ecology with R"
