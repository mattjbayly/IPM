#title: "Hierarchical Spatial Sampling: Plots, Sites, Region"
#author: "Matthew Bayly"
#date: "Thursday, January 8, 2015"

# Individual plants are nested within plots, which themselves are 
# nested within sites. Need to account for nested structure in the
# data and explore where the most variation originated (ie. between plots
# or across sites. 


### _Set directories for computer_ ###
########################### 
## set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
setwd(path.set)
source("01_SetDirectories.R")
# LIBRARIES
library(IPMpack)
library(lme4)
require(ggplot2)
require(lattice)

########################### 
#Open 2014 plant datafile 
setwd(path.dat.raw)
plantdat <- read.csv(file="plantdata.csv") # raw data
dim(plantdat)
setwd(path.dat)
d <- read.csv(file="Data_2014.csv") # revised manipulated dataframe from "02_State_Variable_15.01.08.R"
colnames(d) 
d$size_ln <- log(d$size + 1)

#============================================================================================#
#  SPECIFY STARTING POTENTIAL (TEMPORARY - COMPRAMISE BETWEEN GROWTH & SURV)
#============================================================================================#

d$start <- 1 + ((0.3425*d$size_ln + 0.28032*d$size_ln)/2) + ((0.7916*d$pot + 0.38542*d$pot)/2)

#============================================================================================#
#  MODEL SIMPLIFICATION IN HEIRARICAL SAMPLING
#============================================================================================#

# CAN WE LEAVE OUT THE RANDOM EFFECT OF SITE (is there a significant difference)
# USING LMER
library(lme4)
GSpring1 = lmer(MID.total.height_ln ~ size_ln + pot + (1|site/Uplot), data=d, REML=FALSE)
GSpring2 = lmer(MID.total.height_ln ~ size_ln + pot + (1|Uplot), data=d, REML=FALSE)

anova(GSpring1, GSpring2)
# no site level effect is significant (but not considering enviro effects yet). 

### compare different fixed effect structures
	 # full facotrial interaction
m1 = lmer(MID.total.height_ln ~ start + Uplot + (1|site/Uplot), data=d, REML=FALSE)
m2 = lmer(MID.total.height_ln ~ start + (1|site/Uplot), data=d, REML=FALSE)
m3 = lmer(MID.total.height_ln ~ 1 + (1|site/Uplot), data=d, REML=FALSE)
anova(m1, m2, m3)

# Model with unique slope and intercept for each plot
int_and_slope = lmer(MID.total.height_ln ~ start + (start|site/Uplot), data=d, REML=FALSE)
#abline(allvars)
# Model with only uniue intercept for each plot
int_only = lmer(MID.total.height_ln ~ start + (1|site/Uplot), data=d, REML=FALSE)
anova(int_and_slope, int_only)

################################################################################
#  RANDOM EFFECT STRUCTURE
# Can we drop site?
	# 1. RANDOM INTERCEPT: regression shifted up or down based on random effect of plot nested within site. 
		int_and_slope1 = lmer(MID.total.height_ln ~ start + (1|site/Uplot), data=d, REML=FALSE)
	# 2. RANDOM SLOPE & INTERCEPT: In addtion to the random intercept, there is also now a random effect for slope. 
				# this mean the rate of growth differs between plots, assuming fixed effect is positive
				# plants in plots grow more or less quickly than the global mean depending on whether 
				# the random effect is positive or negative.
		int_and_slope2 = lmer(MID.total.height_ln ~ start + (1 + start|site/Uplot), data=d, REML=FALSE)
	# 3. RANDOM SLOPE & INTERCEPT: I
		int_and_slope3 = lmer(MID.total.height_ln ~ start*site + (start|site/Uplot:site), data=d, REML=FALSE)

		lmer(ERPindex ~ practice*context + (practice|participants) + 
                (practice|participants:context), data=base) 
		
		
		int_and_slope_site_main = lmer(MID.total.height_ln ~ start + site + (start|site/Uplot), data=d, REML=FALSE)
anova(int_and_slope1, int_and_slope2, int_and_slope_site_main)



#============================================================================================#
#  ERROR CHECKING PLOT RESIDUALS
#============================================================================================#

d2 <- d[which(d$start > 0 & d$MID.total.height_ln > 0), ]
int_and_slope = lmer(MID.total.height_ln ~ start + (start|site/Uplot), data=d2, REML=FALSE)
#plot(int_and_slope) #useless
hs <- groupedData(MID.total.height_ln ~ start|site/Uplot, outer=~start, data=d2)
par(mfrow=c(1,3))
plot(int_and_slope, start~resid(.))
plot(int_and_slope, site~resid(.))
plot(int_and_slope, Uplot~resid(.))

# check for heteroscedasticity
plot(int_and_slope, resid(.,type="pearson")~fitted(.)|site)

# normality plots
qqnorm(resid(int_and_slope), main="Q-Q plot for residuals")
qqnorm(ranef(int_and_slope)$Uplot$"(Intercept)", main="Q-Q plot for the random effect" )


#============================================================================================#
#  VISUALIZE PLOTS
#============================================================================================#

int_and_slope = lmer(MID.total.height_ln ~ start + (start|site/Uplot), data=d2, REML=FALSE)

par(mfrow=c(1,1))
mycoefs <- coef(int_and_slope)
sites <- c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK")
plot(d$start, d$MID.total.height_ln, xlab="Start", ylab="JULY size", 
          col=as.factor(d$site), pch=19, 
          cex=0.2); abline(0, 1)

# LM FIT LINES FOR SITES (red)
for(i in 1:length(sites)){		  
	temp <- sites[i]
	d3 <- d[which(d$site==temp), ] 
	z <- lmer(MID.total.height_ln ~ start + (start|Uplot), data=d3) 
	frammy <- do.call(rbind.data.frame, coef(z)[1])
	print(c(temp, mean(frammy[ ,1]), mean(frammy[ ,2])))
		for(j in 1:dim(frammy)[1]){
			# plot all intercepts & slopes for plot level 
			#abline(frammy[j, 1], frammy[j, 2], col=as.factor(d3$site), lty=3, lwd=0.2)
			}
	abline(mean(frammy[ ,1]), mean(frammy[ ,2]), col=as.factor(d3$site), lwd=4)
}



#============================================================================================#
#  VARIANCE AMOUNG AND VARIANCE WITHIN...
#============================================================================================#

z1 <- lmer(MID.total.height_ln ~ start + (start|site/Uplot) + (1|site/Uplot), data=d) 
z2 <- lmer(MID.total.height_ln ~ start + (start|site/Uplot), data=d) 
z3 <- lmer(MID.total.height_ln ~ start + site + (1|site/Uplot), data=d) 
confint(z3, method="Wald")
VarCorr(z3)


# Extract Variance and Correlation Components
VarCorr(z2)
z2.confint <- confint(z2, method="Wald")
