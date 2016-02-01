# ADJUST STRUCTURE OF RAW TRANSPLANT DATAFRAME
# TO WORK WITH EXISTING SCRIPTS CREATED FOR IPM MODELS
# WITH THE DEMOGRAPHY DATASET 

# Matthew - August 27th 2015 (Thursday) 
# Description: The original data file output from the seasonal 
# transplant data needs to be adjusted to match the following output for the demography 
# dataset. 

# IDEAL STRUCTURE: 
#Surv	z	z1	Repr	Fec	SiteID	PlotID	Region	ID	NewPlot_11	NewPlot_12	Year
#1	0.5	0.5	0	NA	Little North Fork Middle Fork	74	C3	14843	FALSE	FALSE	1011
#1	0.5	0.5	0	NA	Little North Fork Middle Fork	74	C3	14873	FALSE	FALSE	1011
#1	0.5	0.5	0	NA	Little North Fork Middle Fork	74	C3	54524	TRUE	FALSE	1112
#1	0.5	0.5	0	NA	Little North Fork Middle Fork	74	C3	54525	TRUE	FALSE	1112

# EXISTING STRUCTURE 
#ID	site	Uplot	surv0	surv1	surv2	surv3	surv4	surv_end15	start	size1	size2	size3	size4 pFlower2014 fec2014 pFlower2015 fec2015 
#tran_518	CALAPOOIA	CALAPOOIAP01	1	1	1	0	NA	0	NA	11	22	NA	NA
#tran_522	CALAPOOIA	CALAPOOIAP01	1	0	NA	NA	NA	NA	NA	NA	NA	NA	NA
#tran_525	CALAPOOIA	CALAPOOIAP01	0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
#tran_523	CALAPOOIA	CALAPOOIAP01	0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA

	# In order to fix the data format the following changes will need to be made. 
		# size4 = z1 (size at time t + 1)
		# size2 = z  (size at time t) 
		# Surv = surv_end15 (or surv4) (survival to time t + 1)
		# pFlower2014 = Repr (flowering at time t) 
		# fec2014 = Fec (fecundity at time t) 
		# pFlower2015 = Repr (flowering at time t, but need to set size4 to z for these rows and need to set surv and z1 to "NA") 
		# fec2015 = Fec (fecundity at time t, but need to set size4 to z for these rows and need to set surv and z1 to "NA")
#---------------------------------------------------------------------------------------


### SET DIRECTORIES & LOAD FILES ###
	path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
	setwd(path.set)
	source("00_SetDirectories.R") # directory script (edit for your own computer). 
	setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)
	# Open 2015/2014 plant datafiles 
	setwd(path.dat); dir()
	Data2015 <- read.csv(file="Data_2015.csv")
	dim(Data2015); colnames(Data2015)
	head(Data2015) # looks good ;)

## working directory must be set here, so the source()'s below run
	setwd("C:/Users/DW/Desktop/transplant_analysis/Planning_Docs/2.IPM_tutorials/Rees_2014_how to IPM/Reese example")
# run the utility functions
	source("./Standard Graphical Pars.R")
## run the ungulate IBM
	source("./Ungulate Demog Funs.R") # but will not use these. 
	setwd(path.dat)
## Load extra libraries
	require(doBy)
	require(car)
	require(mgcv)
	library(coefplot2)	
	

#---------------------------------------------------------------------------------------


### BUILD BASIC IPM DATAFRAME AND ADD 2014-2015 MEASURMENTS (standard t & t+1) 
	# We will need to split up data for the separate years & reformat. We can use the 
	# 2015 July data can be repurposed for the flowering model (probability of flowering based on size) and 
	# the fecundity model (# of seeds, given flowered, based on size). This extra data will be year = 2015 - 2016
	# although no observations come from 2016 so we set Surv & z1 to NA for this data. 
	##############################################################
	# Create empty dataframe according to ideal structure for IPM. 
	IPMData <- data.frame(matrix(NA, nrow = 1, ncol = 12))
	colnames(IPMData) <- c("Surv","z","z1","Repr","Fec","SiteID","PlotID","Region","ID","NewPlot_13","NewPlot_14","Year")
	# Add 2014 to 2015 standard single time step measurements to IPMData. 
	tempFrame <- Data2015[,c("surv_end15", "size2", "size4", "pFlower2014", "fec2014", "site", "Uplot", "ID")]
	tempFrame$Region <- "NA" # will have at add within RL & beyond RL later
	tempFrame$NewPlot_13 <- 0
	tempFrame$NewPlot_14 <- 1
	tempFrame$Year <- "2014to2015"
	# need to IPMData out plants that did died prior to the official t1 measurement in Sept of 2014. 
	before <- dim(tempFrame)[1]
	tempFrame <- tempFrame[!is.na(tempFrame$surv_end15),]
	after <- dim(tempFrame)[1]
	before - after # number of plants that died in the summer of 2014
	after # number of plants remaining for model. 
	# Adjust column names to match conventional IPM format 
	colnames(tempFrame) <- c("Surv","z","z1","Repr","Fec","SiteID","PlotID","ID","Region","NewPlot_13","NewPlot_14","Year")
	tempFrame <- tempFrame [,c("Surv","z","z1","Repr","Fec","SiteID","PlotID","Region", "ID","NewPlot_13","NewPlot_14","Year")] # adjust ordering of ID & region columns 
	head(tempFrame, 2)
	str(tempFrame) # looks good ;) 
	# Adjust "Region" column names according to whether site is within or beyond the RL. 
	tempFrame$Region <- ifelse(tempFrame$SiteID == "CALAPOOIA" | tempFrame$SiteID == "WILEY" | tempFrame$SiteID == "HUNTER" | tempFrame$SiteID == "THOMAS","Beyond","Within")
	head(tempFrame)
	tempFrame$Region <- as.factor(tempFrame$Region); levels(tempFrame$Region) # looks ok ;)
	IPMFrame20142015 <- tempFrame
	
	#########################################################################################################
	# Add on extra plant sizes from the 2015 measurements for fecundity models (remember to keep z1 & Surv as NA)
	# so we are just changing the year here as a factor variable & keeping Surv & z2 as "NA" 
	tempFrame2 <- Data2015[,c("size4", "pFlower2015", "fec2015", "site", "Uplot", "ID")]
	tempFrame2$Region <- "NA" # will have at add within RL & beyond RL later
	tempFrame2$NewPlot_13 <- 0
	tempFrame2$NewPlot_14 <- 1
	tempFrame2$Year <- "2015to2016"
	# need to IPMData out plants that died prior to our t1 measurement here of July 2015. 
	before <- dim(tempFrame2)[1]
	tempFrame2 <- tempFrame2[!is.na(tempFrame2$size4),] # use size4 here since it will also exclude NA values for fec2015 & pFlower2015
	after <- dim(tempFrame2)[1]
	before - after # number of plants that did not flower and/or died prior to the final census 
	after # number of plants that flowered in the final census. 
	# Adjust column names to match conventional IPM format 
	colnames(tempFrame2) <- c("z","Repr","Fec","SiteID","PlotID","ID","Region","NewPlot_13","NewPlot_14","Year")
	tempFrame2$Surv <- "NA"
	tempFrame2$z1 <- "NA"
	tempFrame2 <- tempFrame2[,c("Surv","z","z1","Repr","Fec","SiteID","PlotID","Region","ID","NewPlot_13","NewPlot_14","Year")] # adjust ordering of ID & region columns 
	head(tempFrame2, 2)
	str(tempFrame2) # looks good ;) 
	dim(tempFrame2)
	# Adjust "Region" column names according to whether site is within or beyond the RL. 
	tempFrame2$Region <- ifelse(tempFrame2$SiteID == "CALAPOOIA" | tempFrame2$SiteID == "WILEY" | tempFrame2$SiteID == "HUNTER" | tempFrame2$SiteID == "THOMAS","Beyond","Within")
	head(tempFrame2)
	tempFrame2$Region <- as.factor(tempFrame2$Region); levels(tempFrame2$Region) # looks ok ;)
	IPMFrame20152016 <- tempFrame2

	#########################################################################################################
	# Now merge the two data frames together with all data from the 2014 - 2015 transitions & t1 measurments for the partial 
	# 2015 to 2016 transitions. 
	IPMData <- rbind(IPMFrame20142015, IPMFrame20152016)
	#str(IPMData)

	# Need to fix the variable type, some got mixed incorrectly factor/numeric/character.
	IPMData$Surv <- strtoi(IPMData$Surv) # convert to numeric from character 
	IPMData$z1 <- as.numeric(IPMData$z1)
	IPMData$Year <- as.factor(IPMData$Year)
	head(IPMData, 2)
	summary(IPMData)
	str(IPMData) # looks ok ;)
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------


### BASIC QA/QC CHECKS LOOKING FOR WEIRD ANOMALIES BY PLOTTING BASIC VITAL RATES 	
	# Look for weird residuals, inspect basic plots visually ect. 
	# will probably need to log transform size
	setwd(paste(path.fig, "/01_DataStruct", sep=""))
	pdf(file="01B_OutlierPlots.pdf", width=11, height=8.5)
	# Basic boxplots 
		par(mfrow=c(2,2))
		hist(IPMData$z, main="z - total stem length (cm)")
		hist(IPMData$z1, main="z1 - total stem length (cm)")
		hist(log(IPMData$z), main="z - log transformed")
		hist(log(IPMData$z1), main="z1 - log transformed")
	# log transformation will be necessary, but first explore plots without the transformation
		library(car)
	# Initial size boxplot
		par(mfrow=c(2,2))
		Boxplot(IPMData$z, IPMData$Region, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Region", ylab="Initial Size")
		Boxplot(IPMData$z, IPMData$Site, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Site", ylab="Initial Size")
		Boxplot(log(IPMData$z), IPMData$Region, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Region", ylab="Initial Size")
		Boxplot(log(IPMData$z), IPMData$Site, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Site", ylab="Initial Size")
	# Final size boxplot
		par(mfrow=c(2,2))
		Boxplot(IPMData$z1, IPMData$Region, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Region", ylab="Final Size")
		Boxplot(IPMData$z1, IPMData$Site, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Site", ylab="Final Size")
		Boxplot(log(IPMData$z1), IPMData$Region, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Region", ylab="Final Size")
		Boxplot(log(IPMData$z1), IPMData$Site, labels=IPMData$ID, data=IPMData, id.method="y",
			xlab="Site", ylab="Final Size")
	# Basic plots of growth	
		par(mfrow=c(1,2))
		plot(IPMData$z, IPMData$z1, xlab="Initial Size", ylab="Final Size", bg=as.factor(IPMData$SiteID), pch=21, main="basic growth")
				abline(a = 0, b = 1, col="red", lty=2)
		plot(log(IPMData$z), log(IPMData$z1), xlab="Initial Size", ylab="Final Size", bg=as.factor(IPMData$SiteID), pch=21, main="log(basic growth(")
			abline(a = 0, b = 1, col="red", lty=2)
	# Basic plots of fecundity 
		par(mfrow=c(1,3))
		plot(IPMData$z, IPMData$Fec, xlab="Initial size", ylab="fecundity", bg=as.factor(IPMData$SiteID), pch=21, main="fecundity of a given size z")
		plot(log(IPMData$z), IPMData$Fec, xlab="Log - size", ylab="fecundity", bg=as.factor(IPMData$Site), pch=21, main="fecundity of a given size log(z)")
		plot(log(IPMData$z), log(IPMData$Fec), xlab="Log - size", ylab="Log - fecundity", bg=as.factor(IPMData$Site), pch=21, main="log(fecundity) of a given size log(z)")
 	
	
##################################################	
# Basic plots of Reproduction (0/1)
		par(mfrow=c(1,2))
		drop <- IPMData[!is.na(IPMData$z),]
		drop <- drop[!is.na(drop$Repr),]
				drop <- drop[order(drop$z),]	# sort smallest -> largest (z) 
				z <- drop$z # stand along object 
				# Plotting frame 
					Repr.plot.dataF <- within(drop, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = drop)
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataF)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot.range <- c(0, 2000)
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="")
			points(z, jitter(drop$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			lines(fitted(mod.Repr.glm) ~ z, data = drop, col = "red")
			#############################################################
			# try with a log transformation 
			drop$z <- log(drop$z)
			z <- drop$z # stand along object 
				# Plotting frame 
					Repr.plot.dataF <- within(drop, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = drop)
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataF)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot.range <- c(0, 8)
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("log(Size) , "*italic(z)), ylab = "Prob Reproductive", main="")
			points(z, jitter(drop$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			lines(fitted(mod.Repr.glm) ~ z, data = drop, col = "red")

##################################################	
# Basic plots of survivorship (0/1)		
		par(mfrow=c(1,2))
		drop <- IPMData[!is.na(IPMData$z),]
		drop <- drop[!is.na(drop$Surv),]
				drop <- drop[order(drop$z),]	# sort smallest -> largest (z) 
				z <- drop$z # stand along object 
				# Plotting frame 
					Surv.plot.dataF <- within(drop, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			mod.Surv.glm <- glm(Surv ~ z  , family = binomial, data = drop)
			Surv.ps <- summaryBy(z + Surv ~ z.classes, data = Surv.plot.dataF)
			Surv.ps # (actual values for plotting vs. predicted) 
			plot.range <- c(0, 2000)
			plot(Surv.mean ~ z.mean, data = Surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Survival", main="")
			points(z, jitter(drop$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			lines(fitted(mod.Surv.glm) ~ z, data = drop, col = "red")
			#############################################################
			# try with a log transformation 
			drop$z <- log(drop$z)
			z <- drop$z # stand along object 
				# Plotting frame 
					Surv.plot.dataF <- within(drop, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			mod.Surv.glm <- glm(Surv ~ z  , family = binomial, data = drop)
			Surv.ps <- summaryBy(z + Surv ~ z.classes, data = Surv.plot.dataF)
			Surv.ps # (actual values for plotting vs. predicted) 
			plot.range <- c(0, 8)
			plot(Surv.mean ~ z.mean, data = Surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("log(Size) , "*italic(z)), ylab = "Prob Survival", main="")
			points(z, jitter(drop$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			lines(fitted(mod.Surv.glm) ~ z, data = drop, col = "red")
#
##
###
#### CLOSE PLOTTING DEVICE 
###
##
#
			
dev.off()

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# SAVE FINAL DATA FRAME:
setwd(path.dat); dir()
write.csv(IPMData, file="IPMData_transplant.csv", row.names = FALSE)
