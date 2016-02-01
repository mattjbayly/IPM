#title: "MEM TUTORIAL"
#author: "Matthew Bayly"
#date: "Thursday, January 8, 2015"

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
##########
d$surv2 <- d$FALL.total.height # assuming only AUGUST census
d$surv2[d$surv2>0] <- 1 # Make 0/1
##########
d$surv1 <- d$MID.total.height # assuming only JULY census
d$surv1[d$surv1>0] <- 1 # Make 0/1
##########
d$surv0 <- d$post.trans.condition # assuming only JUNE census
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
d$size_ln <- log(d$start.height)
d$sizeNext_ln <- log(d$FALL.total.height)
##########

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
pairs.panels(d[c('size_ln', 'pot', 'start.bsd', 'MID.stem.number', 'MID.max.stem','MID.total.height_ln', 'surv1')], 
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

setwd(path.dat.raw)


# DIFFICULT TO SEE STRONG POSITIVE TREND POT VS SIZE
#### END DECIDE STARTING POTENTIAL (size plot multiplier)

###########################
###########################
###########################
###########################
###########################
###########################



setwd(path.dat)
write.csv(d, file="Data_2014.csv")
### END DATA INITIALIZATION

###########################
###########################
###########################
###########################
###########################
###########################
