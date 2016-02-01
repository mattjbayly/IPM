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
# set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
setwd(path.set)
source("00_SetDirectories.R") # directory script (edit for your own computer). 
# LIBRARIES
library(IPMpack)
library(lme4)
require(ggplot2)
require(lattice)


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

	#################################################################################
	# Plot sites only 
	setwd(path.fig)
	pdf(file="MAP_herb.reps.pdf", width=11, height=8.5)
	par(mfrow=c(1,2))
	myBiplot(fit, choices=1:2, type = "n", expand=1, xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	# Plot out site colors
	points(fit$x[  1: 13, 1:2], pch=22, col="black", cex=2, bg="darkgoldenrod1") # TOM
	points(fit$x[ 14:27, 1:2], pch=21, col="black", cex=2, bg="darkgoldenrod4") # HUNT
	points(fit$x[28:41, 1:2], pch=23, col="black", cex=2, bg="aquamarine4") # WILEY
	points(fit$x[42:49, 1:2], pch=24, col="black", cex=2, bg="aquamarine") # Calapooia
	points(fit$x[50:63, 1:2], pch=25, col="black", cex=2, bg="blueviolet") # Mosby
	points(fit$x[64:78, 1:2], pch=21, col="black", cex=2, bg="blue") # Coast
	points(fit$x[79:90, 1:2], pch=22, col="black", cex=2, bg="brown4") # Look
	points(fit$x[91:105, 1:2], pch=23, col="black", cex=2, bg="brown2") # rock

	#################################################################################
	# Within vs beyond range plots 
	myBiplot(fit, choices=1:2, type = "n", expand=1, 
			xlab="PC1 28.7%", ylab="PC2 19.0%", yaxt='n', xaxt='n')  # no points
	# Plot out site colors
	points(fit$x[  1:49, 1:2], pch=24, cex=1.5, col="black", bg="aquamarine") # Beyond
	points(fit$x[50:105, 1:2], pch=21, cex=1.5, col="black", bg="brown2") # Within
	dev.off(); dev.off() # close site level plots



#============================================================================================#
# Plot-level vital rates plotted in environmental space 
#============================================================================================#

# Plot-level vital rates in PCA space
		#Open 2014 plant datafile 
		setwd(path.dat)
		# revised dataframe from previous script
		d <- read.csv(file="Data_2014.csv")
		library(nlme); library(lme4)
									completeFun <- function(data, desiredCols) {
										completeVec <- complete.cases(data[, desiredCols])
										return(data[completeVec, ])}	
										dtmp <- completeFun(d, c("size2_ln", "start"))
	# Make vital rate models that include only size as a predictor variable
	# Just want to visualize plot level effects so nothing else! ~ do not nest with site here either	

						Gr <- lmer(size2_ln ~ start + (1|Uplot), REML=FALSE, 
										data=d, na.action=na.omit) # growth given survivorship 
										
						Fec <- glmer(fec ~ size2_ln + (1|Uplot),
										family = poisson, data=dtmp, na.action=na.omit) # fecundity given plant flowered

						Flr <- glmer(pFlower ~ size2_ln + (1|Uplot), 
										family=binomial(link="logit"),
										na.action=na.omit, data=dtmp)				 # probability of flowering (given survive and growth to s2) 
										
						Sur <- glmer(surv_end ~ start + (1|Uplot), 
										family=binomial(link="logit"),
										na.action=na.omit, data=d)					 # survival		
						##############				
					# from: http://jaredknowles.com/journal/2014/5/17/mixed-effects-tutorial-2-fun-with-mermod-objects

						library(arm) 
						methods(class="merMod"); fixef(Sur)
						#confint(Sur, level = 0.95) # sample CI 
						fixef(Sur) / sigma(Sur) # standardize effect size of predictor variables
							# for later: 
					# extract ranef from models
								models <- c("Gr","Sur","Flr","Fec")
								for(i in 1:length(models)){
								temppy <- get(models[i])
								plot_frame <- ranef(temppy)[1] # random effect estimate of levels (intercepts)
								plot_frameb <- data.frame(ID=substr(names(plot_frame),1,1), Obs=plot_frame); plot_frameb$names <- rownames(plot_frameb)
								colnames(plot_frameb) <- c("del", "ranef", "Uplot"); plot_frameb <- plot_frameb[,c("ranef", "Uplot")]
								colnames(plot_frameb)[1] <- paste(models[i], "_ranef", sep="")
								assign(paste(models[i], "_ranef", sep=""), plot_frameb)
							}
		# need to consider if ranefs are really the best deception of plot level effect 
		# over for eg average residuals or other measure?
			# merge dfs together 
					levels(enviro$SITE)[levels(enviro$SITE)=="HUNT"] <- "HUNTER" # fix naming error
					enviro$Un.plot <- paste(enviro$SITE, enviro$PLOT, sep="") # unique plot variable to mere with 
					enviro_a <- merge(Sur_ranef, enviro, by.x = "Uplot", by.y = "Un.plot", all.x = TRUE, all.y = TRUE)
					enviro_b <- merge(Gr_ranef, enviro_a, by.x = "Uplot", by.y = "Uplot", all.x = TRUE, all.y = TRUE)
					enviro_c <- merge(Flr_ranef, enviro_b, by.x = "Uplot", by.y = "Uplot", all.x = TRUE, all.y = TRUE)
					enviro_new <- merge(Fec_ranef, enviro_c, by.x = "Uplot", by.y = "Uplot", all.x = TRUE, all.y = TRUE)
					enviro_new <- droplevels( enviro_new[-which(enviro_new$SITE== "WILD"), ] ) # drop out wild observations
	
			# Redo PCA 
				library(MASS)
				fit <- prcomp(enviro_new[, 9:18], scale. = TRUE) # check other settings	
				var_exp <- (fit$sdev)^2 / sum(fit$sdev^2) # Percent variance explained by each axis

				setwd(path.fig)	
				#pdf(file="PCA_vital_Plotranef_all.pdf", width=11, height=8.5)
				par(mfrow=c(2,2)) # plot all vital rates on one page
				pdf(file="PCA_vital_PLOTranef_pannels.pdf", width=11, height=8.5)
				# Plot and color by renf of plots
					vitals <- c("Sur_ranef", "Gr_ranef", "Flr_ranef", "Fec_ranef")
					my_titles <- c("Survival", "Growth (Given: Survivorship)", "Flowering (Given: Fall size)", "Fecundity (Given: Fall size)")
						for(i in 1:length(vitals)){
						r <- range(enviro_new[,c(vitals[i])]); n <- 20; cv <- rev(topo.colors(n)) # make color range
						c <- cv[cut(enviro_new[,c(vitals[i])],n)]  
						myBiplot(fit, choices=1:2, type = "p", pch=21, col="black", 
								 xaxt='n', yaxt='n', 
								 xlab=paste("PCA1 ", round(((var_exp[1])*100), digits = 1),"%", sep=""), 
								 ylab=paste("PCA2 ", round(((var_exp[2])*100), digits = 1),"%", sep=""), 
								bg=c, col.arrows = "#FF6600", main=my_titles[i])				
						}
				dev.off(); dev.off(); dev.off()  # close all plots

###########################
###########################
###########################
###########################
###########################
###########################

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
pdf(file="Barplot_microsite_vars.pdf", width=11, height=8.5)
par(mfrow=c(3,3))

# moisture
z <- lm(moist_score~SITE,data=enviro_sor)
#str(summary(z)); str(anova(z))
r_sq <- summary(z)$adj.r.squared
ppp <- anova(z)[5]; ppp <- round(ppp[1,1], 4)
a1 <- aov(moist_score~SITE,data=enviro_sor)
tuk <- glht(a1, linfct = mcp(SITE = "Tukey"))
#summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Moisture (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Aspect - Southness (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Cover (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Slope (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Mean grain size (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Flow Seasonality (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Bank Position (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Vert to Wat (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
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
plot(tuk.cld, xaxt='n', col="lightgrey", ylab="",
	xlab=paste("Horiz to Wat (r^2 =", round(r_sq, 2), ", p= ", ppp,")",sep=""))
	axis(1, 1:8, c("L", "R", "C", "M", "C", "W", "H", "T"))

dev.off(); dev.off() # close site level plots
	
###########################
###########################
###########################
###########################
###########################
###########################
				