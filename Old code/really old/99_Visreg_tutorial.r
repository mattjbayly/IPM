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
		
###############################
library(visreg)
par(mfrow=c(1,2))
visreg(Mgr, "ENSEMBLE", by="moist_score", breaks=4, overlay=TRUE,
				scale="response", partial=TRUE, jitter=TRUE, 
				points=list(pch = c(0,6,12,1,7,13,2,8,14,3,9,4,10,5,11,0)[as.numeric(as.factor(d$plot))], col="darkgrey", cex=0.5),
				ylab="Relative Growth Rate", xlab="Predicted Site Suitability From ENM",
				ylim=c(3, 4.2))
				
visreg(Msu, "ENSEMBLE", by="moist_score", breaks=4, overlay=TRUE,
				scale="response", partial=TRUE,
				points=list(pch = c(0,6,12,1,7,13,2,8,14,3,9,4,10,5,11,0)[as.numeric(as.factor(d$plot))], col="darkgrey", cex=0.5),
				ylab="Probability of Survival", xlab="Predicted Site Suitability From ENM", 
				ylim=c(0,1))

###############################
par(mfrow=c(1,2))

visreg(Mfl, "ENSEMBLE", by="moist_score", breaks=4, overlay=TRUE,
				scale="response", partial=TRUE, jitter=TRUE,
				points=list(pch = c(0,6,12,1,7,13,2,8,14,3,9,4,10,5,11,0)[as.numeric(as.factor(d$plot))], col="darkgrey", cex=0.5),
				ylab="Probability of Flowering", xlab="Predicted Site Suitability From ENM", ylim=c(0,1))

visreg(Mfr, "ENSEMBLE", by="moist_score", breaks=4, overlay=TRUE,
				scale="response", partial=TRUE, jitter=TRUE, 
				points=list(pch = c(0,6,12,1,7,13,2,8,14,3,9,4,10,5,11,0)[as.numeric(as.factor(d$plot))],
				cex=0.5, col="darkgrey"),
				ylab="Relative Fecundity", xlab="Predicted Site Suitability From ENM",
				ylim=c(2,10))

###############################
###############################
###############################
###############################
###############################
###############################
###############################
###############################
###############################
###############################




# SITE PLOT
Mgr <- glm(size2_ln ~ start_sd + site + moist_score, data=d, na.action=na.omit)
visreg(Mgr, "start_sd", by = "site")


# 3D PLOT
fit1 <- lm(size2_ln ~ start_sd + I(start_sd^2) + Densi_cover + moist_score
			+ moist_score*Densi_cover + moist_score*start_sd + I(Densi_cover^2) + moist_score*start_sd*I(Densi_cover^2), data=d, na.action=na.omit)
			
		
visreg2d(fit1, "Densi_cover", "moist_score", plot.type = "image")
visreg2d(fit1, "Densi_cover", "moist_score", plot.type = "persp")
visreg2d(fit1, "Densi_cover", "moist_score", plot.type = "rgl")


visreg(fit1, "start_sd")
visreg(fit2, "start_sd", trans = exp)
visreg(fit2, "start_sd")

#
visreg(Mgr, "SUBSTRATE", type = "contrast")

visreg(Mgr, "SUBSTRATE", type = "conditional")
