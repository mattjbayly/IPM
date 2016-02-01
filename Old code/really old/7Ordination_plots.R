#title: "ORDINATION PLOTS"
#author: "Matthew Bayly"
#date: "Tuesday, November 25, 2014"
#**HOW DIFF & SIMILAR ARE PLOTS**

# SOURCE CODE:
# http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf

### _Set directories for computer_ ###
########################### 
## set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"
	
setwd(path.set)
source("00_SetDirectories.R")

# LIBRARIES
library(vegan)
library(MASS)
library(psych)
library(PerformanceAnalytics)
library(rgl)
########################### 
#Open 2014 plant datafile 
setwd(path.dat.raw)
enviro <- read.csv(file="environmental_variables.csv")
dim(enviro);
#colnames(enviro)

###########################
###########################

# START ADJUST DATAFRAME
# retain only useful variables
enviro <- enviro[ , c("SITE","PLOT","SUBSTRATE","Asp_Expo","Densi_cover","sand","boulder","moist_score","Vert._dist_wat","Horiz_dist_to_wat","tot_area_Seasonal","plot_slope","WatAREA__to_bankAREA")] 
summary(enviro) # do all variable values seem ok?
library(PerformanceAnalytics)
chart.Correlation(enviro[,4:13]) # range only includes numerical data
# down to 13 - 2 = 11 causal enviro variables
# replace main substrate to quantitative median particle size (in mm)
enviro$SUBSTRATE <- sub("silt", "0.016", enviro$SUBSTRATE)
enviro$SUBSTRATE <- sub("sand", "0.354", enviro$SUBSTRATE)
enviro$SUBSTRATE <- sub("gravel", "8", enviro$SUBSTRATE)
enviro$SUBSTRATE <- sub("cobble", "43", enviro$SUBSTRATE)
enviro$SUBSTRATE <- sub("boulder", "150", enviro$SUBSTRATE)
enviro$SUBSTRATE <- sub("slab", "500", enviro$SUBSTRATE)
enviro$SUBSTRATE <- as.numeric(enviro$SUBSTRATE)
#summary(enviro); str(enviro)
# END ADJUST DATAFRAME

###########################
###########################
###########################

### START PRELIMINARY EXPLORATION PLOTS
library(psych)
pairs.panels(enviro[,4:8], bg=c("blue","red","yellow","pink","purple","green","black","orange")[enviro$SITE], pch=21,lm=TRUE)
pairs.panels(enviro[8:ncol(enviro)], bg=c("blue","red","yellow","pink","purple","green","black","orange")[enviro$SITE], pch=21,lm=TRUE)

# Transformations needed on variables
# SUBSTRATE; boulder; sand
enviro$SUBSTRATE <- log(enviro$SUBSTRATE + 1)
enviro$sand <- log(enviro$sand + 1)
enviro$boulder <- log(enviro$boulder + 1)

### END PRELIMINARY EXPLORATION PLOTS

###########################
###########################
###########################
###########################
###########################
###########################

### START VEGAN: AN INTRODUCTION TO ORDINATION
# Jari Oksanen Nov 17, 2014

library(vegan)
library(MASS)
enviro2 <- enviro[, 3:ncol(enviro)]
### DCA - Detrended Correspondance Analysis
ord <- decorana(enviro2)
ord; plot(ord)
### NMDS - Non-metric multidimensional scaling
ord <- metaMDS(enviro[, 3:ncol(enviro)])
#### BUILD PLOT from ordination
plot(ord, type="n")
points(ord, display = "sites", cex=0.8, pch=21, col="red", bg=as.factor(enviro$SITE)) # sites just all plots
text(ord, display="spec", cex=0.7, col="blue")
# very high stress = 0.23!; poor representation of data

#######################################################
#### Create environmental matrix (groups; sites ect.)
enviro.env <- enviro[, c("SITE","moist_score")]
mod <- cca(enviro2 ~ SITE, enviro.env)
plot(mod, type="n", scaling = 3)
points(mod, display = "sites", scaling = 3, cex=0.8, pch=21, col="red", bg=as.factor(enviro$SITE)) # sites just all plots
## Catch the invisible result of ordihull...
ordihull(mod, enviro.env$SITE, scaling = 3, label = TRUE)

###################
attach(enviro.env)
ordihull(mod, SITE, scaling = 3) 
ordiellipse(mod, SITE, scaling = 3, col=3, lwd=2)
plot(mod, type="n", scaling = 3)
ordispider(mod, SITE, scaling = 3, col="red", label=TRUE)
points(mod, display = "sites", scaling = 3, cex=0.8, pch=21, col="red", bg=as.factor(enviro$SITE)) # sites just all plots
ordiellipse(mod, SITE, scaling = 3, col=3, lwd=2)

###################
# Fit enviromental variables to data
ord <- metaMDS(enviro2)
ord.fit <- envfit(ord ~ SITE + moist_score, data=enviro.env, perm=999)
plot(ord, type="n", scaling = 3)
plot(ord.fit, scaling = 3)
ordisurf(ord, moist_score, add=T, scaling=3)

###################
# Constrained Ordination
ord <- cca(enviro2 ~ SITE + moist_score, data=enviro.env)
cca(enviro2 ~ ., data=enviro.env)

###################
# Signifigance Test 
anova(ord)
anova(ord, by="term")
anova(ord, by="mar")
anova(ord, by="axis")

###################
# Constrained Ordination 
ord <- cca(enviro2 ~ SITE + Condition(moist_score), data=enviro.env)
anova(ord, by="term")

### END VEGAN: AN INTRODUCTION TO ORDINATION
# Jari Oksanen Nov 17, 2014


###########################
###########################
###########################
###########################
###########################
###########################

### START "Multivariate Analysis of Ecological Communities in R: vegan tutorial"
# Jari Oksanen Feb 8, 2013

library(vegan)
library(MASS)
library(MASS)
enviro2 <- enviro[, 3:ncol(enviro)]
enviro.env <- enviro[, c("SITE","moist_score")]
attach(enviro.env)
##########################
# 2. Ordination basic method
	# 2.1 Non-metric multidimensional scaling
		# problem with row 78 & 81
		vare.dis <- vegdist(enviro2)
		vare.mds0 <- metaMDS(vare.dis)
		stressplot(vare.mds0)
		
		
###########################
###########################
###########################
###########################
###########################
###########################

# TUTORIAL FROM R TIPS PAGE
#D. Schluter, University of BC
#https://www.zoology.ubc.ca/~schluter/R/multivariate/

fit <- princomp(enviro[, 3:ncol(enviro)], cor=TRUE) # correlation matrix rather than covariance matrix
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
par(mfrow=c(1,2))
plot(fit,type="lines") # scree plot 
#fit$scores # the principal components
biplot(fit)

# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation 
fit <- factanal(enviro[, 3:ncol(enviro)], 3, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(enviro[, 3:ncol(enviro)]),cex=.7) # add variable names

# PCA Variable Factor Map 
library(FactoMineR)
result <- PCA(enviro[, 3:ncol(enviro)]) 
result$eigen[2]
barplot(result$eig[,1],main="Eigenvalues",names.arg=1:nrow(result$eig))
summary(result)
