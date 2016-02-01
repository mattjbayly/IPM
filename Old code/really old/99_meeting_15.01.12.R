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
source("00_SetDirectories.R")
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
#  VITAL RATES FROM IPM-PACK
#============================================================================================#

z1 <- lmer(MID.total.height_ln ~ start + (start|site/Uplot) + (1|site/Uplot), data=d) 
z2 <- lmer(MID.total.height_ln ~ start + (start|site/Uplot), data=d) 
z3 <- lmer(MID.total.height_ln ~ start + site + (1|site/Uplot), data=d) 
confint(z3, method="Wald")
VarCorr(z3)

################################################################################################
# from July (take out any post transplant failures
d$sizeNext_ln <- log(d$sizeNext + 1)
fromJul <- d[which(d$surv0 > 0), ] 
par(mfrow=c(2,2))
hist(fromJul$start); hist(fromJul$MID.total.height_ln); 
hist(fromJul$sizeNext_ln)
fromJul$surv <- fromJul$surv1
fromJul$size <- fromJul$start
fromJul$sizeNext <- fromJul$MID.total.height_ln


# VISUAL PLOTS SURVIVAL
par(mfrow=c(2,4))
for(i in 1:length(sites)){		  
	temp <- sites[i]
	d3 <- fromJul[which(fromJul$site==temp), ] 
	plot(d3$start,jitter(d3$surv1), xlab="Start", ylab="Surv to July", pch=19, cex=0.5, main=temp) # jittered
		z <- glmer(surv1 ~ start + (1|Uplot), data = d3, family = binomial, 
			control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
	frammy <- do.call(rbind.data.frame, coef(z)[1])
	}
dev.off()

########################
# mesh points
x<-seq(from=0,to=10,length=1001)
x0<-data.frame(size=x,size2=x*x)
minSize<-min(d$start,na.rm=T)
maxSize<-max(d$start,na.rm=T)

# Global - SURVIVAL
	survModelComp(dataf = fromJul, makePlot = TRUE, legendPos = "bottomright",
	mainTitle = "Survival")
	dev.off()
# SITE SPECIFIC SURVIVAL 
par(mfrow=c(2,4))
		for(i in 1:length(sites)){		  
		temp <- sites[i]
		d3 <- fromJul[which(fromJul$site==temp), ] 
		survModelComp(dataf = d3, makePlot = TRUE, legendPos = "bottomright",
		mainTitle = paste("Survival:", temp, sep=" "))
		# fit MEM, to account for plots
		z1 <- glmer(surv1 ~ 1 + (1|Uplot), data = d3, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
		z2 <- glmer(surv1 ~ start + (1|Uplot), data = d3, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
		z3 <- glmer(surv1 ~ start + start*start + (1|Uplot), data = d3, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
		}
	# models tested: 
		#surv ~ 1       # poor model 
		#surv ~ size	# best
		#surv ~ size	# does not converge at all sites... 
	

########################################################################
# Global - GROWTH
	dev.off()
	growthModelComp(dataf = fromJul, makePlot = TRUE, legendPos = "bottomright",
	mainTitle = "Growth")
	dev.off()
# SITE SPECIFIC GROWTH 
par(mfrow=c(2,4))
		for(i in 1:length(sites)){		  
		temp <- sites[i]
		d3 <- fromJul[which(fromJul$site==temp), ] 
		growthModelComp(dataf = d3, makePlot = TRUE, legendPos = "bottomright",
		mainTitle = paste("Growth:", temp, sep=" "))
		# fit MEM, to account for plots
		z1 <- lmer(sizeNext ~ 1 + (1|Uplot), data = d3, REML=FALSE)
		z2 <- lmer(sizeNext ~ start + (1|Uplot), data = d3, REML=FALSE)
		z3 <- lmer(sizeNext ~ start + start*start + (1|Uplot), data = d3, REML=FALSE)
		}
	# models tested:
		#sizeNext ~ 1       # poor model 
		#sizeNext ~ size	# best
		#sizeNext ~ size^2	# does not converge at all sites...

		
		
		
# NOT ACCOUNTING FOR PLOTS or SITES (GLOBAL P-MATRIX) 
dev.off()
so<-makeSurvObj(fromJul,surv~size)
go<-makeGrowthObj(fromJul,sizeNext~size+size2)

Pmatrix<-makeIPMPmatrix(survObj=so,growObj=go,
minSize=minSize,
maxSize=maxSize)

require(fields); dev.off()
image.plot(Pmatrix@meshpoints,
Pmatrix@meshpoints,
t(Pmatrix),
main = "Pmatrix: survival and growth",
xlab = "Size at t",
ylab = "Size at t+1")
abline(0,1,lty=2,lwd=3)
diagnosticsPmatrix(Pmatrix, growObj=go, survObj=so, correction="constant")
dev.off()

# SITE SPECIFIC P-MATRIX
par(mfrow=c(2,4))
		for(i in 1:length(sites)){		  
		temp <- sites[i]
		d3 <- fromJul[which(fromJul$site==temp), ] 
			so1<-makeSurvObj(d3,surv~size)
			go1<-makeGrowthObj(d3,sizeNext~size+size2)
			Pmatrix<-makeIPMPmatrix(survObj=so1,growObj=go1,
				minSize=minSize,
				maxSize=maxSize)
		image.plot(Pmatrix@meshpoints,
			Pmatrix@meshpoints,
			t(Pmatrix),
			main = paste("Pmatrix:", temp, sep=" "),
			xlab = "Size at t",
			ylab = "Size at t+1")
			abline(0,1,lty=2,lwd=3)
		}
		
########################################################################
#######################################################################
# GLOBAL - FECUNDITY
fromJul$flower <- fromJul$fec.seed
fromJul$flower <- ifelse(fromJul$flower>0.1,1,0)

fo1<-makeFecObj(fromJul, Formula=flower~1, Family = "binomial") # Intercept only model
fo2<-makeFecObj(fromJul, Formula=flower~size, Family = "binomial")
fo3<-makeFecObj(fromJul, Formula=flower~size+size2, Family = "binomial")

fs <- order(fromJul$size)
fs.fec <- (fromJul$flower)[fs]
fs.size <- (fromJul$size)[fs]
pfz <- tapply(fs.size, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE)
ps <- tapply(fs.fec, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE)
plot(as.numeric(pfz), as.numeric(ps), pch = 19, cex=2, col="blue",ylim=c(0,1),
xlab="size", ylab="proportion flowering", main="")
y0<-predict(fo1@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="red")
y0<-predict(fo2@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="green")
y0<-predict(fo3@fitFec[[1]],newdata=x0,type="response")
lines(x,y0,col="blue")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"),
c("AIC"),c(AIC(fo1@fitFec[[1]]),AIC(fo2@fitFec[[1]]),AIC(fo3@fitFec[[1]]))),
col = c(2:4),lty = 1, xjust = 1, bg = "white")

	#LOOKS LIKE sizeNext ~ size^2
	# models tested:
		#sizeNext ~ 1       # poor model 
		#sizeNext ~ size	# second best
		#sizeNext ~ size^2	# BEST, but maybe not enough data per site
		
		
# NUMBER OF FRUITS, CONDITIONAL ON FLOWERING
fo1<-makeFecObj(fromJul, Formula=c(flower~size+size2,fec.seed~1),
Family = c("binomial","poisson"), Transform = c("none","none"))
fo2<-makeFecObj(fromJul, Formula=c(flower~size+size2,fec.seed~size),
Family = c("binomial","poisson"), Transform = c("none","-1"))
fo3<-makeFecObj(fromJul, Formula=c(flower~size+size2,fec.seed~size+size2),
Family = c("binomial","poisson"), Transform = c("none","-1"))
