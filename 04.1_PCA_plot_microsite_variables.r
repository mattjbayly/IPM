#title: "PCA of plots with plot-level environmental variables"
#author: "Matthew Bayly"
#modified: "Saturday, January 17, 2015"

# PURPOSE OF SCRIPT: 
# From the plot-level environmental variables collected:
	#i.	Were any sites significantly different from each other? 
		#(run perMANOVA or MRPP on plots with site as a factor)
	#ii. Were the selected transplant plots similar to environmental 
		# variables collected from wild M. cardinalis? 
		#(use MRPP because design heavily unbalanced; wild 
		# individual measurements (~ 20 – 30/creek) nested 
		# within specific creeks (Coast, Rock, Calapooya, Canton)).
		
# If only a very small number of plots are extreme environmental outliers
# then it may be easier to exclude them rather than trying to 
#incorporate microsite variables into final analysis. 

# INDEX: 
	# 0. load libraries & set directories to computer
	# 1. Load plot-level data (check details)
	# 2. Explore attributes of environmental data 
		# 2.1 Correlation matrix (normality met?)
		# 2.2 Make necessary transformations for normality
		# 2.3 Normalize data (is spread too great).
	# 3. Run a simple PCA for all plots
		# 3.1 visualize sites (colour)
		# 3.2 visualize regions (colour border *or use pch) 
	# 4. Diagnose PCA
		# 4.1 Variance explained per axis (histogram)
		# 4.2 
		# 4.3  
		# ..... 
		
	# 5. Site differences 
		# 5.1 Were any sites significantly different from each other overall (perMANOVA)
		# 5.2 For each site, which variable contributed most 
			# strongly to difference from the mean. 
		# 5.3 Which site was most variable & which site most uniform?
	
#============================================================================================#
# 0. load libraries & set directories to computer
#============================================================================================#
### SET DIRECTORIES & LOAD FILES ###
	path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
	setwd(path.set)
	source("00_SetDirectories.R") # directory script (edit for your own computer). 
	setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)
	# Open 2015/2014 plant datafiles 
	setwd(path.dat); dir()
	dfclean <- read.csv(file="IPMData_transplant.csv")
	dim(dfclean); colnames(dfclean)
	head(dfclean) # looks good ;)
	# log transformation of basic size data 
	dfclean$z <- log(dfclean$z ); summary(dfclean$z)
	dfclean$z1 <- log(dfclean$z1); summary(dfclean$z1)
	# round down fecundity data to actual fruit counts:
	dfclean$Fec <- round(dfclean$Fec, digits = 0) # round 
	dfclean$Fec <- as.integer(dfclean$Fec) # make integer for Poisson regression later.
	plot.range <- c(0, 9) # range of size data to be plotted (~ 0 - 9 log(total stem length))
 
# LIBRARIES
library(IPMpack)
library(lme4)
require(ggplot2)
require(lattice)

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
colnames(enviro)

###########################
###########################

	# START ADJUST DATAFRAME
	# retain only useful variables
	enviro <- enviro[ , c("SITE","PLOT","SUBSTRATE","Asp_Expo","Densi_cover","Particle_mean","Particle_even","moist_score","Vert._dist_wat","Horiz_dist_to_wat","tot_area_Seasonal","plot_slope","WatAREA__to_bankAREA")] 
	summary(enviro) # do all variable values seem ok?
	library(PerformanceAnalytics)
	chart.Correlation(enviro[,4:13]) # range only includes numerical data
	# "Asp_Expo" ~ no trans needed           
	# "Densi_cover" ~ no trans needed        
	# "Particle_mean" ~ log trans***************
	# "Particle_even"~ sqrt trans **************     
	# "Vert._dist_wat" ~ log trans**************    
	# "Horiz_dist_to_wat"~ no trans needed 
	# "tot_area_Seasonal" ~ log trans***********    
	# "plot_slope" ~ no trans needed           
	# "WatAREA__to_bankAREA" ~ log trans********   

	# TRANSFORMATIONS NEEDED FOR NORMALITY
	enviro$Particle_mean <- log(enviro$Particle_mean + 1)
	enviro$Particle_even <- sqrt(enviro$Particle_even + 1)
	enviro$Vert._dist_wat <- log(enviro$Vert._dist_wat + 1)
	enviro$tot_area_Seasonal <- log(enviro$tot_area_Seasonal + 1)
	enviro$WatAREA__to_bankAREA <- log(enviro$WatAREA__to_bankAREA + 1)

	chart.Correlation(enviro[,4:13]) # range only includes numerical data
	# NORMAL?
	shapiro.test(enviro$Asp_Expo)
	shapiro.test(enviro$Densi_cover)
	shapiro.test(enviro$Particle_mean)
	shapiro.test(enviro$Particle_even)
	shapiro.test(enviro$moist_score)
	shapiro.test(enviro$Vert._dist_wat)
	shapiro.test(enviro$Horiz_dist_to_wat)
	shapiro.test(enviro$tot_area_Seasonal)
	shapiro.test(enviro$plot_slope)
	shapiro.test(enviro$WatAREA__to_bankAREA)

	summary(enviro); str(enviro)
	# END ADJUST DATAFRAME

	# NEED TO NORMALIZE - PCA BASED ON CORRELATION MATRIX (NOT COVARIANCE MATRIX)
	# http://www.researchgate.net/post/What_is_the_best_way_to_scale_parameters_before_running_a_Principal_Component_Analysis_PCA
	# Since each of the variables were measured on difference scales, we should standardize or normalize
	# There is no dominant species driving community assemblage here
	# but rather than subtracting the mean and dividing by the standard deviation 
	# we can just look at the correlation matrix insted 


#============================================================================================#
# PCA Plots across sites within & beyond range & in relation to variables 
#============================================================================================#

	# PCA OF CORRELATION MATRIX 
	# drop out wild observations:
	enviro <- subset(enviro, SITE != "WILD")
	
	fit <- princomp(enviro[, 4:13], cor=TRUE) # correlation matrix rather than covariance matrix (variables different units)
	summary(fit) # print variance accounted for 
	loadings(fit) # pc loadings 
	par(mfrow=c(1,2))
	plot(fit,type="lines") # scree plot 
	#fit$scores # the principal components
	biplot(fit, expand=0.5)

	# load custom PCA plotting function
	setwd(path.funct)
	# thanks to Boris!
	source("my_Biplot.R") # from Boris https://stat.ethz.ch/pipermail/r-help//2014-April/374098.html
	library(MASS)
	fit <- prcomp(enviro[, 4:13], scale. = TRUE) # correlation matrix rather than covariance matrix (variables different units)
	# scale. scales variables to have unit variance (mean = 0, stdev =1)
	var_exp <- (fit$sdev)^2 / sum(fit$sdev^2) # Percent variance explained by each axis
	myBiplot(fit, choices=1:2, type = "t", pch=20, col="red", col.arrows = "#FF6600")

	#################################################################################
	# color by moisture score
	r <- range(enviro[,8]); n <- 15; cv <- rev(topo.colors(n)) # make color range
	c <- cv[cut(enviro[,8],n)]  
	myBiplot(fit, choices=1:2, type = "p", pch=20, col=c, col.arrows = "#FF6600")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	#################################################################################
	# Plot sites only 
	setwd(path.fig)	
	pdf(file = "04.1_MAP_herb.reps.pdf", width=(8.5 - (1.25 + 0.75)), height=(3.7), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/

	par(mfrow=c(1,2))
	myBiplot(fit, choices=1:2, type = "n", expand=1, xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	# Plot out site colors
	points(fit$x[  1: 13, 1:2], pch=22, col="black", cex=1, bg="darkgoldenrod1") # TOM
	points(fit$x[ 14:27, 1:2], pch=21, col="black", cex=1, bg="darkgoldenrod4") # HUNT
	points(fit$x[28:41, 1:2], pch=23, col="black", cex=1, bg="aquamarine4") # WILEY
	points(fit$x[42:49, 1:2], pch=24, col="black", cex=1, bg="aquamarine") # Calapooia
	points(fit$x[50:63, 1:2], pch=25, col="black", cex=1, bg="blueviolet") # Mosby
	points(fit$x[64:78, 1:2], pch=21, col="black", cex=1, bg="blue") # Coast
	points(fit$x[79:90, 1:2], pch=22, col="black", cex=1, bg="brown4") # Look
	points(fit$x[91:105, 1:2], pch=23, col="black", cex=1, bg="brown2") # rock

	#################################################################################
	# Within vs beyond range plots 
	myBiplot(fit, choices=1:2, type = "n", expand=1, 
			xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	# Plot out site colors
	points(fit$x[  1:49, 1:2], pch=24, cex=1, col="black", bg="aquamarine") # Beyond
	points(fit$x[50:105, 1:2], pch=21, cex=1, col="black", bg="brown2") # Within
	dev.off() # close site level plots
	
# repeat above with wmf
	win.metafile(filename = "04.1_MAP_herb.reps.wmf", width = (8.5 - (1.25 + 0.75)), height = 3.7, family="Times")
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(1,2))
	myBiplot(fit, choices=1:2, type = "n", expand=1, xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	points(fit$x[  1: 13, 1:2], pch=22, col="black", cex=1, bg="darkgoldenrod1") # TOM
	points(fit$x[ 14:27, 1:2], pch=21, col="black", cex=1, bg="darkgoldenrod4") # HUNT
	points(fit$x[28:41, 1:2], pch=23, col="black", cex=1, bg="aquamarine4") # WILEY
	points(fit$x[42:49, 1:2], pch=24, col="black", cex=1, bg="aquamarine") # Calapooia
	points(fit$x[50:63, 1:2], pch=25, col="black", cex=1, bg="blueviolet") # Mosby
	points(fit$x[64:78, 1:2], pch=21, col="black", cex=1, bg="blue") # Coast
	points(fit$x[79:90, 1:2], pch=22, col="black", cex=1, bg="brown4") # Look
	points(fit$x[91:105, 1:2], pch=23, col="black", cex=1, bg="brown2") # rock
	myBiplot(fit, choices=1:2, type = "n", expand=1, xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	points(fit$x[  1:49, 1:2], pch=24, cex=1, col="black", bg="aquamarine") # Beyond
	points(fit$x[50:105, 1:2], pch=21, cex=1, col="black", bg="brown2") # Within
	dev.off() # close site level plots
	
	# add in a correlation chart
		chart.Correlation(enviro[,4:13]) # range only includes numerical data
	
	

#============================================================================================#
# Plot-level vital rates plotted in environmental space 
#============================================================================================#

# load models
	library(nlme)
	library(lme4)
	setwd(path.obj)
	par(mfrow=c(2,2))
	surv_model <- get(load("SurvMod6.rda")) #
		mm <- ranef(surv_model)		
		mm <- mm[1]
		mm <- mm$PlotID
		colnames(mm) <- "Int"
		hist(mm$Int)
		mm$Join <- rownames(mm)
		SurvInt <- mm
		colnames(SurvInt) <- c("SurvInt", "Join")
	grow_model <- get(load("GrowMod.rda")) #	
		mm <- ranef(grow_model)		
		mm <- mm[1]
		mm <- mm$PlotID
		colnames(mm) <- "Int"
		hist(mm$Int)
		mm$Join <- rownames(mm)
		GrowInt <- mm
		colnames(GrowInt) <- c("GrowInt", "Join")
	Repr_model <- get(load("Reprmod9L.rda")) 
		mm <- ranef(Repr_model)		
		mm <- mm[1]
		mm <- mm$PlotID
		colnames(mm) <- "Int"
		hist(mm$Int)
		mm$Join <- rownames(mm)
		ReprInt <- mm
		colnames(ReprInt) <- c("ReprInt", "Join")
	Fec_model <- get(load("FecMod_NB.rda")) 
		mm <- ranef(Fec_model)		
		mm <- mm[1]
		mm <- mm$PlotID
		hist(mm)
		FecInt <- mm
		FecInt <- data.frame(FecInt)
		FecInt$Join <- rownames(FecInt)
		colnames(FecInt) <- c("FecInt", "Join")
	
	
# Merge intercepts to enviro frame
	levels(enviro$SITE)[levels(enviro$SITE)=="HUNT"] <- "HUNTER"
	enviro$Join <- paste0(enviro$SITE, enviro$PLOT)

	# join int values together
	enviro <- merge(enviro, SurvInt, by.x='Join', by.y='Join', all.y=TRUE, all.x=TRUE)
	enviro <- merge(enviro, GrowInt, by.x='Join', by.y='Join', all.y=TRUE, all.x=TRUE)
	enviro <- merge(enviro, ReprInt, by.x='Join', by.y='Join', all.y=TRUE, all.x=TRUE)
	enviro <- merge(enviro, FecInt, by.x='Join', by.y='Join', all.y=TRUE, all.x=TRUE)
	enviro_new <- enviro
	
	
# SAVE FINAL PLOTS
				library(MASS)
				setwd(path.fig)	
	# get legend fxn 
	setwd(path.funct)
	source("legend.col.R")
		
# start plot
setwd(path.fig)	
	pdf(file = "04.2_PCA_vital_PLOTranef_pannels.pdf", width=(8.5 - (1.25 + 0.75)), height=(11-2-(0.75*2)), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
		
		par(mfrow=c(2,2)) # plot all vital rates on one page
				# Plot and color by renf of plots
					vitals <- c("SurvInt", "GrowInt", "ReprInt", "FecInt")
					my_titles <- c("Survival", "Growth (Given: Survivorship)", "Flowering (Given: Fall size)", "Fecundity (Given: Fall size)")
						for(i in 1:length(vitals)){
						r <- range(enviro_new[,c(vitals[i])], na.rm=TRUE); n <- 20; cv <- rev(topo.colors(n)) # make color range
						c <- cv[cut(enviro_new[,c(vitals[i])],n)]  
						colr <- rev(topo.colors(n))
						c[is.na(c)] <- "lightgrey"
						myBiplot(fit, choices=1:2, type = "p", pch=21, col="black", 
								 xaxt='n', yaxt='n', 
								 xlab=paste("PCA1 ", round(((var_exp[1])*100), digits = 1),"%", sep=""), 
								 ylab=paste("PCA2 ", round(((var_exp[2])*100), digits = 1),"%", sep=""), 
								bg=c, col.arrows = "#FF6600", main=my_titles[i])	
						
						# add on legend
						nonNa <- enviro_new[,c(vitals[i])]
						nonNa <- nonNa[!is.na(nonNa)]
						legend.col(col = colr, lev = nonNa)
						}
				dev.off()# close all plots

				
		
				
###########################
###########################
###########################
###########################
###########################
###########################

#============================================================================================#
# FIND VARIABLE MOST RESPONSIBLE FOR PLOT LEVEL DIFFERENCES
#============================================================================================#

PCs <- fit$x
PC1 <- PCs[,1]
PC2 <- PCs[,2]
PC3 <- PCs[,3]
enviro$PC1 <- PC1 
enviro$PC2 <- PC2 
enviro$PC3 <- PC3 

library(MASS)
# SURVIVAL 
	MSurv <- lm(SurvInt~SUBSTRATE + Asp_Expo + Densi_cover + Particle_mean + 
					Particle_even + moist_score + Vert._dist_wat + 
					Horiz_dist_to_wat + tot_area_Seasonal + plot_slope + 
					WatAREA__to_bankAREA + PC1 + PC2 + PC3,
	data=enviro, na.action=na.omit)
	step1Surv <- stepAIC(MSurv , direction="both")
	MSurv <- lm(SurvInt~moist_score + PC1 + PC2, data=enviro, na.action=na.omit)
	step2Surv <- stepAIC(MSurv, direction="backward")
# GROWTH 
	MGrow <- lm(GrowInt~SUBSTRATE + Asp_Expo + Densi_cover + Particle_mean + 
					Particle_even + moist_score + Vert._dist_wat + 
					Horiz_dist_to_wat + tot_area_Seasonal + plot_slope + 
					WatAREA__to_bankAREA + PC1 + PC2 + PC3,
	data=enviro, na.action=na.omit)
	step1Grow <- stepAIC(MGrow , direction="both")
	MGrow <- lm(GrowInt~moist_score + PC1 + PC2, data=enviro, na.action=na.omit)
	step2Grow <- stepAIC(MGrow, direction="backward")
	AIC((lm(GrowInt~moist_score + PC1 + PC2, data=enviro, na.action=na.omit)), (lm(GrowInt~moist_score, data=enviro, na.action=na.omit)))[2]
# Repr
	MRepr <- lm(ReprInt~SUBSTRATE + Asp_Expo + Densi_cover + Particle_mean + 
					Particle_even + moist_score + Vert._dist_wat + 
					Horiz_dist_to_wat + tot_area_Seasonal + plot_slope + 
					WatAREA__to_bankAREA + PC1 + PC2 + PC3,
	data=enviro, na.action=na.omit)
	step1Repr <- stepAIC(MRepr , direction="both")
	MRepr <- lm(ReprInt~moist_score + PC1 + PC2, data=enviro, na.action=na.omit)
	step2Repr <- stepAIC(MRepr, direction="backward")
	AIC((lm(ReprInt~moist_score + PC1 + PC2, data=enviro, na.action=na.omit)), (lm(ReprInt~moist_score, data=enviro, na.action=na.omit)))[2]
#
###
#####
#####
##### USE MOISTURE SCORE AS A KEY VARIABLE!!!!
#####
####
###
##


#============================================================================================#
# BOXPLOTS - of each variable averaged across sites 
#============================================================================================#


# VARIATION IN PLOT VARIABLES ACROSS SITES
setwd(path.dat.raw)
enviro_sor <- read.csv(file="environmental_variables.csv")
enviro_sor <- enviro_sor[ , c("SITE","PLOT","SUBSTRATE","Asp_Expo","Densi_cover","Particle_mean","Particle_even","moist_score","Vert._dist_wat","Horiz_dist_to_wat","tot_area_Seasonal","plot_slope","WatAREA__to_bankAREA")] 
enviro_sor$SITE <- factor(enviro_sor$SITE, levels = c("LOOK", "ROCK", "COAST", "MOSBY", "CALAPOOIA", "WILEY", "HUNT", "THOMAS"))
enviro_sor <- enviro_sor[order(enviro_sor$SITE),]

# set up plotting frame
par(mfrow=c(3,3))
library(agricolae)
library(multcomp)

# Save barplot of variables
setwd(path.fig)
pdf(file = "04.3_Barplot_microsite_vars.pdf", width=(8.5 - (1.25 + 0.75)), height=(8), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(3,3), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(3,3), par(mar=c(4,2,2,1)))

# moisture
z <- lm(Asp_Expo~SITE, data=enviro_sor)
z <- lm(moist_score~SITE, data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(moist_score~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$moist_score, xaxt='n', ylab="",
	main=paste("Moisture (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))
	
# Asp_Expo
z <- lm(Asp_Expo~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(Asp_Expo~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$Asp_Expo, xaxt='n',   ylab="", 
	main=paste("Aspect - Southness (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))

# Densi_cover
z <- lm(Densi_cover~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(Densi_cover~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$Densi_cover, xaxt='n',   ylab="",
	main=paste("Cover (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep="" ), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))

# plot_slope
z <- lm(plot_slope~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(plot_slope~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$plot_slope, xaxt='n',   ylab="",
	main=paste("Slope (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))

# Particle_mean
z <- lm(Particle_mean~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(Particle_mean~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$Particle_mean, xaxt='n',   ylab="", , col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"),
	main=paste("Mean grain size (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))
			
	
# tot_area_Seasonal
z <- lm(tot_area_Seasonal~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(tot_area_Seasonal~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$tot_area_Seasonal, xaxt='n',   ylab="", col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"),
	main=paste("Flow Seasonality (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))
	
# WatAREA__to_bankAREA
z <- lm(log(WatAREA__to_bankAREA+1)~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(log(WatAREA__to_bankAREA+1)~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$WatAREA__to_bankAREA, xaxt='n',   ylab="",
	main=paste("Bank Position (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))
	
# Vert._dist_wat
z <- lm(Vert._dist_wat~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(Vert._dist_wat~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$Vert._dist_wat, xaxt='n',   ylab="",
	main=paste("Vert to Wat (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))
	
# Horiz_dist_to_wat
z <- lm(Horiz_dist_to_wat~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(Horiz_dist_to_wat~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(enviro_sor$SITE, enviro_sor$Horiz_dist_to_wat, xaxt='n',   ylab="",
	main=paste("Horiz to Wat (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""), col=c("red", "red", "red", "red", "blue",  "blue", "blue", "blue"))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))

dev.off(); dev.off() # close site level plots
	
###########################
###########################
###########################
###########################
###########################
###########################
				