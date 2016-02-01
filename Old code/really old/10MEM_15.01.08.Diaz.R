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
###########################

#### START DECIDE STARTING POTENTIAL (size plot multiplier)

# define starting potential 

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
anova(JulStNum1, JulStNum2); anova(JulStNum1, JulStNum3)
JulHi1 = lmer(MID.total.height ~ size + pot + (1|Uplot), data=d, REML=FALSE)
JulHi2 = lmer(MID.total.height ~ pot + (1|Uplot), data=d, REML=FALSE)
JulHi3 = lmer(MID.total.height ~ size + (1|Uplot), data=d, REML=FALSE)
anova(JulHi1, JulHi2); anova(JulHi1, JulHi3)

# DIFFICULT TO SEE STRONG POSITIVE TREND POT VS SIZE
#### END DECIDE STARTING POTENTIAL (size plot multiplier)

###########################
###########################
###########################
###########################
###########################
###########################

# MIXED EFFECT MODELS FOR GROWTH
# SPRING 
GSpring = lmer(MID.total.height ~ size + pot + site + (1|Uplot), data=d, REML=FALSE)
summary(GSpring)
# coef(GSpring)
library(effects)
ef <- effect("site", GSpring)
summary(ef)
#ef$lower; ef$upper; ef$estimate
x <- as.data.frame(ef)
# order sites for plot pannels 
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
######
# SPRING FACTOR LEVEL SITE ESTIMATES FOR SPRING GROWTH (given survivorship)
ggplot(x, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)


# FALL growth, given surv
GFall = lmer(sizeNext ~ MID.total.height + site + (1|Uplot), data=d, REML=FALSE)
summary(GFall)
ef <- effect("site", GFall)
summary(ef)
x <- as.data.frame(ef)
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
######
# SPRING FACTOR LEVEL SITE ESTIMATES FOR SPRING GROWTH (given survivorship)
ggplot(x, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)

######################################
# WERE SITE LEVEL DIFFERENCES SIGNIFICANT?
GSpring1 = lmer(MID.total.height ~ size + pot + site + (1|Uplot), data=d, REML=FALSE)
GSpring2 = lmer(MID.total.height ~ size + pot + (1|Uplot), data=d, REML=FALSE) # drop site
# SPRING DIFFERENCES ~ YES
anova(GSpring1, GSpring2)

GFall1 = lmer(sizeNext ~ MID.total.height + site + (1|Uplot), data=d, REML=FALSE)
GFall2 = lmer(sizeNext ~ MID.total.height + (1|Uplot), data=d, REML=FALSE)
# FALL DIFFERENCES ~ YES
anova(GFall1, GFall2)

###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################

# START SURVIVORSHIP MEMs
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

d0 <- d

# model for spring survivorship
m <- glmer(surv0 ~ size + pot + site +(1 | Uplot), data = d0, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
print(m, corr = FALSE)
se <- sqrt(diag(vcov(m)))
# SITE estimates with 95% CI (Coefficients on **logit** scale)
(tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 * se))
# SITE estimates with 95% CI (Coefficients on **Log-odds ratio** scale)
exp(tab)
ef <- effect("site", m)
x <- as.data.frame(ef)
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# SURVIVORSHIP POST TRANSPLANT 
ggplot(x, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)

###################################################################

# spring survival (remove post-trans mortalities from dataframe)
d1 <- d
d1$surv0[d1$surv0<1] <- NA # Make NA if die
d1 <- d1[complete.cases(d1$surv0),]
###################
# model for SPRING survivorship
m <- glmer(surv1 ~ size + pot + site +(1 | Uplot), data = d1, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m)
se <- sqrt(diag(vcov(m)))
# SITE estimates with 95% CI (Coefficients on **logit** scale)
(tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 * se))
# SITE estimates with 95% CI (Coefficients on **Log-odds ratio** scale)
exp(tab)
ef <- effect("site", m)
x <- as.data.frame(ef)
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# SURVIVORSHIP TO SPRING FROM POST TRANS 
ggplot(x, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)


###################################################################

# fall survival (remove spring mortalities from dataframe)
d2 <- d
d2$surv1[d2$surv1<1] <- NA # Make NA if die
d2 <- d2[complete.cases(d2$surv1),]
###################
# model for FALL survivorship
m <- glmer(surv2 ~ size + pot + site + size*site + (1 | Uplot), data = d2, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
print(m, corr = FALSE)
summary(m)
se <- sqrt(diag(vcov(m)))
# SITE estimates with 95% CI (Coefficients on **logit** scale)
(tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 * se))
# SITE estimates with 95% CI (Coefficients on **Log-odds ratio** scale)
exp(tab)
ef <- effect("site", m)
x <- as.data.frame(ef)
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# SURVIVORSHIP TO FALL FROM SPRING 
ggplot(x, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)


# END SURVIVORSHIP MEMs

###########################
###########################
###########################
###########################
###########################
######


# START FECUNDITY MEMs
(Ffall1 <- glmer(fec.seed ~ sizeNext + site + (1|Uplot), family="poisson",data=d))
(Ffall2 <- glmer(fec.seed ~ sizeNext + (1|Uplot), family="poisson",data=d))
anova(Ffall1, Ffall2)
#####################
ef <- effect("site", Ffall1)
x <- as.data.frame(ef)
x$site <- factor(x$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# FECUNDITY IN FALL (SITE-level effect terms)
ggplot(x, aes(site, fit)) + geom_point() + theme_bw(base_size=12)
x2 <- x[c(-8),]
ggplot(x2, aes(site, fit)) + geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.4) + theme_bw(base_size=12)

# END FECUNDITY MEMs


###########################
###########################
###########################
###########################
###########################
######