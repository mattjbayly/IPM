#title: "MEM TUTORIAL"
#author: "Matthew Bayly"
#date: "Tuesday, November 25, 2014"

### _Set directories for computer_ ###
########################### 
## set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"
	
setwd(path.set)
source("01_SetDirectories.R")
# LIBRARIES
library(IPMpack)
library(psych)
library(lme4)
require(ggplot2)
require(GGally)
require(reshape2)
require(compiler)
require(parallel)
require(boot)

########################### 
#Open 2014 plant datafile 
setwd(path.dat.raw)
plantdat <- read.csv(file="plantdata.csv")
dim(plantdat);
#colnames(plantdat)


#############################
#############################
#############################
### START DATA INITIALIZATION

# remove seeds rows  
plantdat <- plantdat[ !grepl("SEED", plantdat$pot) , ] 
#plantdat <- plantdat[ !grepl("HERB", plantdat$mid.condition) , ]
summary(plantdat$pot) 
################################################
# Set up Merow plant parameters according to Merow TUTORIALs
d <- plantdat # dataframe
d$size <- d$start.height # start height is hight in greenhouse for now
d$sizeNext <- d$FALL.total.height # mid summer height is 
d$sizeNext_ln <- log(d$sizeNext + 1)
##########
d$surv2 <- d$FALL.total.height # assuming only AUGUST census
d$surv2[d$surv2>0] <- 1 # Make 0/1
##########
d$surv1 <- d$MID.total.height # assuming only JULY census
d$surv1[d$surv1>0] <- 1 # Make 0/1
##########
d$surv0 <- d$post.trans.condition # assuming only JULY census
d$surv0 <- sub("A", "1", d$surv0)
d$surv0 <- sub("D", "0", d$surv0)
d$surv0 <- as.numeric(d$surv0)
##########
d$fec.seed <- d$FALL_fruit # sum of fruits in FALL
#fec.seed <- rowSums(d[, c("FALL.FR.devo.", "FALL.FR.full.", "FALL.FR.dehis.")]) # adds fall fruit columns
fec.flower <- rowSums(d[, c("FALL.plus.devo", "FALL.flower.fail")]) # adds fall flower columns
##########
# Include stem numbers
d$MID.stem.number[d$MID.stem.number<1] <- NA # Make NA if die
d$FALL.stem.number[d$FALL.stem.number<1] <- NA # Make NA if die
d$MID.total.height[d$MID.total.height<1] <- NA # Make NA if die in spring too 
d$sizeNext[d$sizeNext<1] <- NA # Make NA if die
##########
d <- cbind(d, fec.flower) # bring them together

# Include other study variable (growth, survive ect).
site <- levels(d$site); site <- as.factor(site) # study sites
# filter down dataframe
d <- d[ ,c('size', 'sizeNext', 'sizeNext_ln', 'FALL.stem.number', 'surv2', 'surv1',
			'surv0', 'fec.seed', 'fec.flower', 'site', 
			'plot', 'source','MID.total.height', 
			'post.trans.condition', 'pot', 'start.bsd', 'MID.stem.number', 'MID.max.stem')]
#########
# order sites for plot pannels 
d$site <- factor(d$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
site <- levels(d$site); site <- as.factor(site) #sites
###########
# Convert pot values to sequential numeric numbers
# replace values
d$pot <- as.character(d$pot)
d$pot[d$pot=="XL"]<-5
d$pot[d$pot=="XS"]<-1
d$pot[d$pot=="S"]<-2
d$pot[d$pot=="L"]<-4 # despite a few outliers model fit better with L as 4
d$pot[d$pot=="M"]<-3
d$pot <- as.numeric(d$pot) # BACK TO NUMERIC
summary(d$pot) 
#######
d$start.bsd <- log(d$start.bsd)
d$MID.total.height_ln <- log(d$MID.total.height)
d$Uplot <- paste(d$site, d$plot, sep="")

### END DATA INITIALIZATION

###########################
###########################
###########################
###########################
###########################
###########################

### START PREDICTOR VARIABLE EXPLORATION plots (STARTING POTENTIAL)

# best predictors from greenhouse to spring growth and survivorship 
library(psych)
pairs.panels(d[c('size', 'pot', 'start.bsd', 'MID.stem.number', 'MID.max.stem','MID.total.height', 'surv1')], 
	bg=c("blue","red","yellow","pink","purple","green","black","orange")[d$site], 
	pch=21, lm=TRUE)
# THOMAS    HUNTER    WILEY     CALAPOOIA MOSBY     COAST     ROCK      LOOK 

### END PREDICTOR VARIABLE EXPLORATION (STARTING POTENTIAL)

###########################
###########################
###########################
###########################
###########################
###########################

#### START EXPLORATORY LINEAR MIXED EFFECTS MODELS

# GROWTH_spring ~ START_potential + SITE + (1|PLOT) + ε
# SURVIRORSHIP_spring ~ START_potential + SITE + (1|PLOT) + ε
### additional quick exploration plots
boxplot(size ~ source*pot,d) # should potentially switch 3 & 4 pot size Medium & Large
boxplot(MID.stem.number ~ source*pot,d) # SHOULD SWITCH
boxplot(MID.total.height_ln ~ source*pot,d) # SHOULD SWITCH
boxplot(MID.total.height ~ source*pot,d) # SHOULD SWITCH, some outliers in large!
boxplot(sizeNext ~ source*pot,d) # Large tended to do a bit better in the fall.
### end additional explor plots

##########################################################################
##### 	QUICK MODEL RUNS (EXPLORATION)
JulStNum1 = lmer(MID.stem.number ~ size + pot + site + (1|Uplot), data=d)
JulStNum2 = lmer(MID.stem.number ~ site + (1|Uplot), data=d)
anova(JulStNum1, JulStNum2) # MUST INCLUDE POT & SIZE FOR STM #
# JulStNum2 - thomas; calapooia; rock; hunt; mosby; wiley; look; coast  
# JulStNum1 - thomas; calapooia; rock; mosby; wiley; hunt; look; coast
JulTHNum_full = lmer(MID.total.height_ln ~ site + pot*size + (1|Uplot), data=d, REML=FALSE) 
JulTHNum1 = lmer(MID.total.height_ln ~ size + pot + site + (1|Uplot), data=d, REML=FALSE) # 1ST BEST, df=12 (8 sites - 1)
JulTHNum2A = lmer(MID.total.height_ln ~ site + pot + (1|Uplot), data=d, REML=FALSE) # 2ND PLACE, 
JulTHNum2B = lmer(MID.total.height_ln ~ site + size + (1|Uplot), data=d, REML=FALSE) # 3RD PLACE (ok), df=7
JulTHNum3 = lmer(MID.total.height_ln ~ site + (1|Uplot), data=d, REML=FALSE) # 4TH PLACE
JulTHNum4 = lmer(MID.total.height_ln ~ (1|Uplot), data=d, REML=FALSE) # WOREST MODEL, df=4
JulTHNum5 = lmer(MID.total.height_ln ~ source + pot + size + site + (1|Uplot), data=d, REML=FALSE) # Plants from Rock creek appear to do slighlty better
JulTHNum6 = lmer(MID.total.height_ln ~ pot + size + site + (1|Uplot), data=d, REML=FALSE) # Plants from Rock creek appear to do slighlty better
anova(JulTHNum5, JulTHNum6)

# ESTIMATE REPRESENTING CHANGE FROM THOMAS TO SPECIFIC GROUP
# JulTHNum1 - thomas; rock; calapooia; hunt ----  mosby; wiley; look; coast
# JulStNum1 - thomas; calapooia; rock; mosby ---- wiley; hunt; look; coast
###############################
# look at likelihood ratio test
anova(JulTHNum_full, JulTHNum1, JulTHNum2A, JulTHNum2B, JulTHNum3, JulTHNum4)
# coef(JulTHNum1)
# RANDOM SLOPE & RANDOM INTERCEPT MODEL
JulTHNum1A = lmer(MID.total.height_ln ~ size + pot + site + (1|Uplot), data=d, REML=FALSE)
JulTHNum1B = lmer(MID.total.height_ln ~ size + pot + site + (1+pot|Uplot), data=d, REML=FALSE)
anova(JulTHNum1A, JulTHNum1B)
# coef(JulTHNum1B) # effects of random slope of pot w/ plot all still positive

#### END EXPLORATORY LINEAR MIXED EFFECTS MODELS

###########################
###########################
###########################
###########################
###########################
###########################

#### START EXPLORATORY MIXED EFFECTS 'LOGISTIC REGRESSION' MODELS
# SOURCE: http://www.ats.ucla.edu/stat/r/dae/melogit.htm

# for survivorship 0-1
 
require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)

# model for spring survivorship
m <- glmer(surv1 ~ size + pot + site +(1 | Uplot), data = d, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
print(m, corr = FALSE)

se <- sqrt(diag(vcov(m)))
# SITE estimates with 95% CI (Coefficients on **logit** scale)
(tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 * se))
# SITE estimates with 95% CI (Coefficients on **Log-odds ratio** scale)
exp(tab)

# Recall variable relationships & check details (look for errors & wierd values)
pairs.panels(d[c("sizeNext_ln", "FALL.stem.number", "MID.total.height_ln", "MID.stem.number", "pot", "size")], 
	bg=c("blue","red","yellow","pink","purple","green","black","orange")[d$site], 
	pch=21, lm=TRUE)
	
#### END EXPLORATORY MIXED EFFECTS 'LOGISTIC REGRESSION' MODELS

###########################
###########################
###########################
###########################
###########################
###########################

#### START D.SCHLUTER FIT MIXED EFFECTS MOD R-TIPS PG.
# source: https://www.zoology.ubc.ca/~schluter/R/fit-model/#lme
##### 	QUICK MODEL RUNS (EXPLORATION)
JulStNum1 = lmer(MID.stem.number ~ size + pot + site + (1|Uplot), data=d, REML=FALSE)
JulStNum2 = lmer(MID.stem.number ~ site + (1|Uplot), data=d, REML=FALSE)
anova(JulStNum1, JulStNum2) # MUST INCLUDE POT & SIZE FOR STM
summary(JulStNum1)                # parameter estimates, fit, tests of fixed effects
plot(JulStNum1)                   # plot of residuals against predicted values

plot(resid(JulStNum1))            # residuals
plot(fitted(JulStNum1))           # best linear unbiased predictors (BLUPs)
VarCorr(JulStNum1)                # variance components for random effects  

anova(JulStNum1)                  # Tests of fixed effects (fitted sequentially)
anova(JulStNum1, type="marginal") # Test of fixed effects using "drop-one" testing
anova(JulStNum1, JulStNum2)       # [don't use this with lme to compare model fits]

library(lme)

StNum1 <- lme(MID.stem.number ~ 1, random= ~ 1 | Uplot, data = d, na.action=na.exclude)
summary(StNum1)
intervals(StNum1)
VarCorr(StNum1)
plot(StNum1)
stripchart(MID.stem.number ~ Uplot, vertical = TRUE, pch = 1, data = d)
stripchart(fitted(MID.stem.number) ~ Uplot, vertical = TRUE, add = TRUE, pch = "---", data = d)
StNum1 <- lme(MID.stem.number ~ 1, random= ~ 1 | site/Uplot, data = d, na.action=na.exclude)
StNum1 <- lme(MID.stem.number ~ site + pot + source + size, random= ~ 1 | Uplot, data = d, na.action=na.exclude)
intervals(StNum1) # 95% confidence intervals
anova(StNum1)                    # for sequential testing
anova(StNum1, type = "marginal") # for drop-one testing

#### END D.SCHLUTER FIT MIXED EFFECTS MOD R-TIPS PG.

###########################
###########################
###########################
###########################
###########################
###########################

#### START DECIDE STARTING POTENTIAL (size plot multiplier)

# define starting potential 
JulStNum = lmer(MID.stem.number ~ size + pot + pot*size + (1|Uplot), data=d, REML=FALSE)
JulHi = lmer(MID.total.height ~ size + pot + pot*size + (1|Uplot), data=d, REML=FALSE)
JulMax = lmer(MID.max.stem ~ size + pot + pot*size + (1|Uplot), data=d, REML=FALSE)

XS <- d[ which(d$pot=='1'), ]
S <- d[ which(d$pot=='2'), ]
M <- d[ which(d$pot=='3'), ]
L <- d[ which(d$pot=='4'), ]
XL <- d[ which(d$pot=='5'), ]
sizes <- c("XS", "S", "M", "L", "XL")

par(mfrow=c(2, 3))
for(i in 1:length(sizes)){
	temp <- get(sizes[i])
		plot(MID.total.height ~ size, pch = 19, data = temp, main=sizes[i], col=as.factor(temp$site))
	rm(temp)
}

JulStNum1 = lmer(MID.stem.number ~ size + pot + (1|Uplot), data=d, REML=FALSE)
JulStNum2 = lmer(MID.stem.number ~ pot + (1|Uplot), data=d, REML=FALSE)
JulStNum3 = lmer(MID.stem.number ~ size + (1|Uplot), data=d, REML=FALSE)
anova(JulStNum1, JulStNum2, JulStNum3)
JulHi = lmer(MID.total.height ~ size + pot + (1|Uplot), data=d, REML=FALSE)
JulMax = lmer(MID.max.stem ~ size + pot + (1|Uplot), data=d, REML=FALSE)

# DIFFICULT TO SEE STRONG POSITIVE TREND POT VS SIZE
#### END DECIDE STARTING POTENTIAL (size plot multiplier)

###########################
###########################
###########################
###########################
###########################
###########################



