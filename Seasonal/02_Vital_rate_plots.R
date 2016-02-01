#---
#title: "2. Vital Rate Visual Exploration Plots"
#author: "Matthew Bayly"
#date: "Sunday, November 02, 2014"
#output: html_document
#---

# INDEX
# 1. SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS
# 2. GROWTH PLOTS: Given survivorship across sites
	# 2A Growth (given survivorship) over entire season from May planting (start) to September census (size2_ln)
	# 2B Spring growth from May planting (start) to July census (size1_ln)
	# 2C Fall growth only from July census to Sept census (size1_ln - size2_ln)
# 3. Survivorship plots across sites
	# 3A Survivorship over the entire season (given post transplant survivorship)  
	# 3B Survivorship to the spring May - JulY
	# 3C Survivorship post-transplant only
# 4. Probability of FLOWERING in the fall, based on initial size & given survivorship 
# 5. Fecundity based on size given plant flowered in the fall (also given survivorship in the fall)




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


#Open 2014 plant datafile 
setwd(path.dat)
# revised dataframe from previous script
d <- read.csv(file="Data_2014.csv")
dim(d); #colnames(plantdat)
# make site level factor
site <- levels(d$site); site <- as.factor(site) # study sites
# order by grouping
site <- factor(c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))

# second order polynomial for quadratic regression 
d$start2 <- d$start*d$start
d$size1_ln2 <- d$size1_ln*d$size1_ln


#============================================================================================#
# 2. GROWTH PLOTS: Given survivorship across sites
#============================================================================================#

##############################################################################################
# 2A Growth (given survivorship) over entire season from May planting (start) to September census (size2_ln)
xaxis <- "Initial starting potential"
yaxis <- "Final ln(size) in Sept"

setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_1Growth_may_to_sept.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$start, d$size2_ln, type="n", xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(0, 6.5), xlim=c(0, 3.5)) # basice fixed empty plotting frame
		points(current$start, current$size2_ln, 
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
		
		# make light grey lines for each site within a pannel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					mod1 <- lm(size2_ln ~ start, data=tempy, na.action=na.omit)
					#mod2 <- lm(size2_ln ~ start + start2, data=tempy, na.action=na.omit)
					#print(anova(mod1, mod2))
					#print(AIC(mod1, mod2))
					# no improvement from quadratic & more support here for linear model (not quadratic)
					abline(lm(size2_ln ~ start, data=tempy, na.action=na.omit), col="lightgrey", lwd=0.8) # fixed lm for each site	
			}
	
	   # make single thick blue line specific site within the pannel
			abline(lm(size2_ln ~ start, data=current, na.action=na.omit), col="blue", lwd=2) # fixed lm for each site	

}
dev.off(); dev.off()
			
##############################################################################################
# 2B Spring growth from May planting (start) to July census (size1_ln)
xaxis <- "Initial starting potential"
yaxis <- "Spring July Size"

setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_2Growth_may_to_july.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$start, d$size1_ln, type="n", xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(0, 6.1), xlim=c(0, 3.5)) # basice fixed empty plotting frame
		points(current$start, current$size1_ln, 
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
		
		# make light grey lines for each site within a pannel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					mod1 <- lm(size1_ln ~ start, data=tempy, na.action=na.omit)
					#mod2 <- lm(size1_ln ~ start + start2, data=tempy, na.action=na.omit)
					#print(anova(mod1, mod2))
					#print(AIC(mod1, mod2))
					# no improvement again from quadratic & more support here for linear model (not quadratic)
					abline(lm(size1_ln ~ start, data=tempy, na.action=na.omit), col="lightgrey", lwd=0.8) # fixed lm for each site	
			}
	
	   # make single thick blue line specific site within the pannel
			abline(lm(size1_ln ~ start, data=current, na.action=na.omit), col="green", lwd=2) # fixed lm for each site	

}
dev.off(); dev.off()			

##############################################################################################
# 2C Fall growth only from July census to Sept census (size1_ln - size2_ln)
xaxis <- "Spring July Size"
yaxis <- "Fall Sept Size"

setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_3Growth_july_to_sept.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$size2_ln, d$size1_ln, type="n", xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(0, 6.5), xlim=c(0, 6.1)) # basice fixed empty plotting frame
		points(current$size1_ln, current$size2_ln, 
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
		
		# make light grey lines for each site within a pannel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					mod1 <- lm(size2_ln ~ size1_ln, data=tempy, na.action=na.omit)
					#mod2 <- lm(size2_ln ~ size1_ln + size1_ln2, data=tempy, na.action=na.omit)
					#print(anova(mod1, mod2))
					#print(AIC(mod1, mod2))
					# again no improvement from quadratic & more support here for linear model (not quadratic)
					abline(lm(size2_ln ~ size1_ln, data=tempy, na.action=na.omit), col="lightgrey", lwd=0.8) # fixed lm for each site	
			}
	
	   # make single thick blue line specific site within the pannel
			abline(lm(size2_ln ~ size1_ln, data=current, na.action=na.omit), col="purple", lwd=2) # fixed lm for each site	

}
dev.off(); dev.off()			
# END: GROWTH PLOTS: Given survivorship across sites
		
	
#============================================================================================#
# 3. START: Survivorship plots across sites
#============================================================================================#
		
##############################################################################################
# 3A Survivorship over the entire season (given post transplant survivorship)  
xaxis <- "Initial starting potential"
yaxis <- "Survivorship to the fall"


setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_4Surv_may_to_sept.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$start, d$surv2, type="n", yaxt='n', xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(-0.1, 1.1), xlim=c(0, 3.5)) # basice fixed empty plotting frame
		points(current$start, jitter(current$surv2, amount=0.05),
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
			axis(2, at=c(0, 1))#

		# make light grey lines for each site within a panel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					mod1 <- glm(surv2 ~ start, family = binomial(link="logit"), data=tempy, na.action=na.omit)
					#mod2 <- glm(surv2 ~ start + I(start^2), family = binomial(link="logit"), data=tempy, na.action=na.omit)
					#print(anova(mod1, mod2))
					#print(AIC(mod1, mod2))
					# again no improvement from quadratic & more support here for linear model (not quadratic)
					s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
					p <- predict(mod1, data.frame(start = s), type = "response") 
					lines(p ~ s, col ="lightgrey") 
			}
		#plot site line for plannel
			mod_curr <- glm(surv2 ~ start, family = binomial(link="logit"), data=current, na.action=na.omit)
			s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
			p <- predict(mod_curr, data.frame(start = s), type = "response") 
					lines(p ~ s, col ="blue", lwd=2) 
}
dev.off(); dev.off()
# END: Whole season survivorship plots across sites

	
##############################################################################################
# 3B Survivorship to the spring May - July
xaxis <- "Initial starting potential"
yaxis <- "Survivorship to the July"

setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_5Surv_may_to_july.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$start, d$surv1, type="n", yaxt='n', xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(-0.1, 1.1), xlim=c(0, 3.5)) # basice fixed empty plotting frame
		points(current$start, jitter(current$surv1, amount=0.05),
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
			axis(2, at=c(0, 1))#

		# make light grey lines for each site within a panel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					if (min(na.omit(tempy$surv1))==0 & max(na.omit(tempy$surv1))==1) { 
							mod1 <- glm(surv1 ~ start, family = binomial(link="logit"), data=tempy, na.action=na.omit)
							s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
							p <- predict(mod1, data.frame(start = s), type = "response") 
							lines(p ~ s, col ="lightgrey") 
					}else{
							print(paste("model fail: ", site[j], sep=""))
							}

			}
			
		#plot site line for panel
					if (min(na.omit(current$surv1))==0 & max(na.omit(current$surv1))==1) { 
						mod_curr <- glm(surv1 ~ start, family = binomial(link="logit"), data=current, na.action=na.omit)
						s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
						p <- predict(mod_curr, data.frame(start = s), type = "response") 
						lines(p ~ s, col ="green", lwd=2) 
					}else{
						print(paste("model fail: ", site[i], sep=""))
							}

			}					
	dev.off(); dev.off()	
# END: Spring only survivorship plots across sites

	
##############################################################################################
# 3C Survivorship post-transplant only
xaxis <- "Initial starting potential"
yaxis <- "Post transplant survivorship"
setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_6Surv_may_to_june.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$start, d$surv0, type="n", yaxt='n', xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(-0.1, 1.1), xlim=c(0, 3.5)) # basice fixed empty plotting frame
		points(current$start, jitter(current$surv0, amount=0.05),
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
			axis(2, at=c(0, 1))#

		# make light grey lines for each site within a panel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					if (min(na.omit(tempy$surv0))==0 & max(na.omit(tempy$surv0))==1) { 
							mod1 <- glm(surv0 ~ start, family = binomial(link="logit"), data=tempy, na.action=na.omit)
							s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
							p <- predict(mod1, data.frame(start = s), type = "response") 
							lines(p ~ s, col ="lightgrey") 
					}else{
							print(paste("model fail: ", site[j], sep=""))
							}

			}
			
		#plot site line for panel
					if (min(na.omit(current$surv0))==0 & max(na.omit(current$surv0))==1) { 
						mod_curr <- glm(surv0 ~ start, family = binomial(link="logit"), data=current, na.action=na.omit)
						s <- seq(min(na.omit(d$start)), max(na.omit(d$start)), by = 0.1) 
						p <- predict(mod_curr, data.frame(start = s), type = "response") 
						lines(p ~ s, col ="red", lwd=2) 
					}else{
						print(paste("model fail: ", site[i], sep=""))
							}

			}					
dev.off(); dev.off()		
# END: Survivorship post-transplant only



#============================================================================================#
# START: Probability of FLOWERING in the fall, based on final size & given survivorship 
#============================================================================================#

##############################################################################################
# Probability of flowering
xaxis <- "Final size in Fall"
yaxis <- "Probability of Flowering"
setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_7Flr_size2.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$size2_ln, d$pFlower, type="n", yaxt='n', xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(-0.1, 1.1), xlim=c(0, 6.5)) # basice fixed empty plotting frame
		points(current$size2_ln, jitter(current$pFlower, amount=0.05),
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 
			axis(2, at=c(0, 1))#

		# make light grey lines for each site within a panel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
					if (min(na.omit(tempy$pFlower))==0 & max(na.omit(tempy$pFlower))==1) { 
							mod1 <- glm(pFlower ~ size2_ln, family = binomial(link="logit"), data=tempy, na.action=na.omit)
							s <- seq(min(na.omit(d$size2_ln)), max(na.omit(d$size2_ln)), by = 0.1) 
							p <- predict(mod1, data.frame(size2_ln = s), type = "response") 
							lines(p ~ s, col ="lightgrey") 
					}else{
							print(paste("model fail: ", site[j], sep=""))
							}

			}
			
		#plot site line for panel
					if (min(na.omit(current$pFlower))==0 & max(na.omit(current$pFlower))==1) { 
						mod_curr <- glm(pFlower ~ size2_ln, family = binomial(link="logit"), data=current, na.action=na.omit)
						s <- seq(min(na.omit(d$size2_ln)), max(na.omit(d$size2_ln)), by = 0.1) 
						p <- predict(mod_curr, data.frame(size2_ln = s), type = "response") 
						lines(p ~ s, col ="red", lwd=2) 
					}else{
						print(paste("model fail: ", site[i], sep=""))
							}

			}					
	dev.off(); dev.off()	
# END: Probability of flowering



#============================================================================================#
# START: Fecundity based on size given plant flowered in the fall (also given survivorship in the fall)
#============================================================================================#

##############################################################################################
# Probability of flowering
xaxis <- "Final size in Fall"
yaxis <- "Fecundity"
setwd(paste0(path.fig, "/02_Vitals"))
pdf(file="02_8Fec_size2.pdf", width=11, height=8.5)
par(mfrow=c(2, 4)) # for 8 sites
for (i in 1:length(site)){
		current <- d[d$site==site[i], ] # for every site 
		plot(d$size2_ln, d$fec, type="n", xlab=xaxis, ylab=yaxis, main=site[i], ylim=c(0, 60), xlim=c(0, 6.5)) # basice fixed empty plotting frame
		points(current$size2_ln, current$fec,
			# annoying notation to get unique pch for each plot
			pch = c(0,6,1,2,3,4,5,7,8,9,10,12,13,14,24,25)[as.numeric(as.factor(as.character(current$Uplot)))], cex=0.75, col="darkgrey") # points of plants 

		# make light grey lines for each site within a panel
			for (j in 1:length(site)){
				tempy <- d[d$site==site[j], ] # for every site 
							mod1 <- glm(fec ~ size2_ln, family = "poisson", data=tempy, na.action=na.omit)
							s <- seq(min(na.omit(d$size2_ln)), max(na.omit(d$size2_ln)), by = 0.1) 
							p <- predict(mod1, data.frame(size2_ln = s), type = "response") 
							lines(p ~ s, col ="lightgrey") 
			}
			
		#plot site line for panel
						mod_curr <- glm(fec ~ size2_ln, family = "poisson", data=current, na.action=na.omit)
						s <- seq(min(na.omit(d$size2_ln)), max(na.omit(d$size2_ln)), by = 0.1) 
						p <- predict(mod_curr, data.frame(size2_ln = s), type = "response") 
						lines(p ~ s, col ="purple", lwd=2)
			}					
	dev.off(); dev.off()
# END: Fecundity based on size given plant flowered in the fall (also given survivorship in the fall)
