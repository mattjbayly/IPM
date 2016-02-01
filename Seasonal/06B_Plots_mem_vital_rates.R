#---
#title: "6B_ Plots of vital rates in mixed effect modelling framework
#author: "Matthew Bayly"
#date: "Sunday, January 02, 2014"
#---

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

# END: SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS
##############################################################################################

#============================================================================================#
# 2. START: MEM simple plots with 'effects' package
#============================================================================================#

library(effects)
setwd(path.funct)
# run all models from the main list
	source("run_all_models.R"); setwd(path.obj)

# examine coefficents from model, but would not want to plot this values since they don't take the random effects into account(?)
	coef(summary(SurBPSW2)) # from the glmmML pretty much same results 'coef(summary(SurBPSW))'

# using the effects packge 
# to convert these parameter estimates into condition mean and SE estimates
	ef <- effect("site_type", SurBPSW2)
	summary(ef) # for the categorical levels 

	# for plotting convert the effect list object into a data frame
	x <- as.data.frame(ef); x
# plot & add on error bars 
	plot(x$site_type, x$fit, ylim=c(0,1), col=c("red", "blue", "blue"), pch=19)
	# add error bars
	arrows(c(1:3), x$lower, c(1:3), x$upper, angle=90, code=3)
	
###########################################################################################
# Jared Knowland Explore MerMod Objects
	# main slots
		slotNames(SurBPSW2) 
	# all options in the mermod object
	methods(class = "merMod")
		#e.g's
		SurBPSW2@call # model call
		head(SurBPSW2@frame)
	# fixed effects of model 
		fixef(SurBPSW2)
		