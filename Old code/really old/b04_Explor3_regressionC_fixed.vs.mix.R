#title: "VITAL RATE REGRESSION FUNCTIONS"
#author: "Matthew Bayly"
#modified: "Thursday, January 16, 2015"

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
# X. Starting state variable (what should be used as start size from greenhouse?). 
#============================================================================================#
d$start <- 0.01*(d$start.height_ln) + 0.05*(d$pot) +  0.1*(d$start.height_ln*d$pot)# multiplicative

#============================================================================================#
# X. SITE SPECIFC FIXED MODEL VS RANDOM EFFECTS MODEL FITTED (for visualization)
#============================================================================================#

##############################################################################################
##############################################################################################
##############################################################################################

	# SURVIVORSHIP PLOT
plot(jitter(d$surv_end, amount=0.05) ~ d$start, 
		pch=19, cex=0.5, data=d, xlab="Start potential", ylab="Survival")



dev.off(); par(mfrow=c(2,4))
for(i in 1:length(site)){
   current <- d[d$site==site[i], ]
   current_post <- d2_post[d2_post$site==site[i], ]
		surv.reg0 = glm(current_post$surv0 ~ current_post$potential, family=binomial(logit), data=current_post)
		surv.reg1 = glm(current$surv1 ~ current$potential, family=binomial(logit), data=current)
		surv.reg2 = glm(current$surv2 ~ current$potential, family=binomial(logit), data=current)
			lines(current$potential, surv.reg2$fitted, type="l", col="red")
			#lines(current$potential, surv.reg1$fitted, type="l", col="blue")
			#lines(current_post$potential, surv.reg0$fitted, type="l", col="green")
		text(x=-5, y=logit(coefficients(surv.reg2)[1]), site[i])
		assign(paste("surv.reg0",site[i], sep=""), surv.reg0)
		assign(paste("surv.reg1",site[i], sep=""), surv.reg1)
		assign(paste("surv.reg2",site[i], sep=""), surv.reg2)			
  }

##############################################################################################  
##############################################################################################  
##############################################################################################
  
# FALL GROWTH LINE FOR EACH SITE
d2_sN <- d2[complete.cases(d2[,c("sizeNext")]),]
d2_sN <- d2_sN[order(d2_sN$potential),]
plot(d2_sN$potential,d2_sN$sizeNext,xlab="Start potential", ylab="FALL Size")
lines(d2_sN$potential, growth.reg$fitted, type="l", col="black")
for (i in 1:length(site)){
   current <- d2_sN[d2_sN$site==site[i], ]
		growth.reg=lm(sizeNext~potential,data=current)
			lines(current$potential, growth.reg$fitted, type="l", col="red")
		assign(paste("growth.reg",site[i], sep=""), growth.reg)
}

##############################################################################################
##############################################################################################
##############################################################################################

# FECUNDITY LINE FOR EACH SITE
plot(d2$potential,d2$fec.seed,xlab="Start potential",ylab="Fall Fruits")
lines(d2$potential, seed.reg$fitted, type="l", col="black")
for (i in 1:length(site)){
   current <- d2[d2$site==site[i], ]
		seed.reg=glm(fec.seed~potential,data=current,family=poisson())
			lines(current$potential, seed.reg$fitted, type="l", col="red")
		assign(paste("seed.reg",site[i], sep=""), seed.reg)
}

