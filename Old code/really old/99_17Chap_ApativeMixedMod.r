# CHAPTER 17 FROM ZURR
#	MAKE MULTIPANNEL GRAPHS OF STATIONS
# 17.4 Estimating Common Patterns Using Additive
# 18 - PLANKTON TIME SERIES 
# 19 - TWO WAY NESTED DATA (similar to chapter 5?)
# 20 - three way nested data. 


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
# 19. CHAPTER 19
#============================================================================================#




























































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




















































				
Sur <- glmer(surv_end ~ start + moist_score + (1|site), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)

sjp.glmer(Sur, type="ri.pc",  show.se= TRUE)






fit2,
          type = "ri.pc",
          show.se = TRUE)
