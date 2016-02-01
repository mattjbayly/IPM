#---
#title: "8_ AIC tables of MEMs with alternative variables
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
library(WriteXLS)


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
# COMPARE MODELS BASED ON AIC (ZURR CHAP 21) 
# make tables
#============================================================================================#

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
tabby$Structure <- "glmmML Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Surv_Uplot <- tabby[order(-tabby$WEIGHT),]

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
tabby$Structure <- "glmmML Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Flr_Uplot <- tabby[order(-tabby$WEIGHT),]


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
tabby$Structure <- "glmmML Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Fec_Uplot <- tabby[order(-tabby$WEIGHT),]


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
tabby$Structure <- "lmer Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Gr_Uplot <- tabby[order(-tabby$WEIGHT),]


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
tabby$Structure <- "glmer Site/Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Surv_Uplot_Site <- tabby[order(-tabby$WEIGHT),]


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
tabby$Structure <- "glmer Site/Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Flr_Uplot_Site <- tabby[order(-tabby$WEIGHT),]
						
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
tabby$Structure <- "glmer Site/Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Fec_Uplot_Site <- tabby[order(-tabby$WEIGHT),]
						
				
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
tabby$Structure <- "lmer Site/Uplot"; tabby <- tabby[,c("Structure", "MODEL", "AIC", "WEIGHT")]
Gr_Uplot_Site <- tabby[order(-tabby$WEIGHT),]
###########################################################################################################				

# reminder index sheet 
header <- c("CODE", "variable", "scale", "description")
a1 <- c("B", "start", "plant", "Start Size (start for growth & surv, but flr & fec use size2_ln")
a2 <- c("P", "moist_score", "plot", "Plot level moisture score")
a3 <- c("S", "ENSEMBLE", "site", "Ensemble SDM from occupancy manuscript")
a4 <- c("W", "site_type", "region", "Categorical for site (beyond range, within range, within range unoccupied)")
var_codes <- rbind(a1, a2, a3, a4); var_codes <- data.frame(var_codes); colnames(var_codes) <- header



#============================================================================================#
# Save tables to excel file
#============================================================================================#
library(WriteXLS)
testPerl(perl = "perl", verbose = TRUE) # make sure pearl is installed properly.
# for windows install proper version of perl from 
#http://www.activestate.com/activeperl/

setwd(path.obj)
# merge glmer & GlmmML tables together for final tables 
Surv <- rbind(Surv_Uplot_Site, Surv_Uplot)
Gr <- rbind(Gr_Uplot_Site, Gr_Uplot)
Flr <- rbind(Flr_Uplot_Site, Flr_Uplot)
Fec <- rbind(Fec_Uplot_Site, Fec_Uplot)

# all the AIC summary tables that will enter the excel spreadsheet. 
# save excel spreadhseet
WriteXLS(c("var_codes", "Surv", "Gr", "Flr", "Fec"), "Vital_AICw_tables.xlsx", c("README", "Survival", "Growth", "Flowering", "Fecundity"), 
BoldHeaderRow = TRUE, AdjWidth=TRUE)


