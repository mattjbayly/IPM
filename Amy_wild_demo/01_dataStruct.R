#============================================================================================#
# ANALYSIS OF WILD CARDINALIS TRANSPLANT DATA 
#============================================================================================#
# April 23, 2015
# SCRIPT 1: Preliminary check of data structure and brief data exploration 
# file name: 01_dataStruct.R


# INDEX: 
#------------
# Check basic data structure & reformat layout for use in basic IPM analysis
# Outliers - boxplots & cleveland dotplots, brief check, but more rigorous check will be preformed during model diagnostics
# Homogeneity - conditional boxplots (again brief & preliminary)
# Normality - histograms and QQ-plots (again brief & preliminary)

#============================================================================================#
# part 1: Set directories, load libraries & import data
# Check basic data structure & reformat layout for use in basic IPM analysis
#============================================================================================#

### _Set directories for computer_ ###
path.root="C:/Users/DW/Desktop/transplant_analysis"	
setwd(path.root)
	path.dat=paste(path.root, "/Data/Amy_wild_demo_data", sep="") # MAIN DATA FILES
	path.code=paste(path.root, "/Rcode/Amy_wild_demo", sep="") # Source Code
	path.obj=paste(path.root, "/Robjects/Amy_wild_demo", sep="") # To store temporary objects
	path.fig=paste(path.root, "/Figures/Amy_wild_demo", sep="") # Polished & exploratory figures
	setwd(path.dat); setwd(path.code); setwd(path.fig); setwd(path.obj)
	path.funct = paste(path.root,"/Rcode/functions", sep="")
	setwd(path.dat); dir()

# load in 2010-2012 data (two time step observations) 
rawdf <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
	summary(rawdf)
	dim(rawdf) # observations = 9123
	
	# Weird things noticed in raw data immediately
		# "TotFr_10" - An individual with 7667 fruits? A few XXXL individuals in Sweet Water & Little Jamison 
		# "PercFr_11" & "PercFr_12" - Percent flower in 2011, there are two individuals with values above 1 (2 & 2.33) this must have been a typo during data entry. 
		# "TotFl_11" - Max = 24,680.60, this is ridiculous for a single individual 
			# will have to deal with the super giants later, wont even use the % flowering for this project 
	
# Reformat data frame for basic IPM setup to match Merow & Reese examples 
		start_keep <- c("A", "J")
		year_end_keep <- c("A", "J", "D")
		yr1011  <- rawdf[rawdf$Class_10 %in% start_keep,]
	# drop individuals excluded in the next year  
		yr1011  <- rawdf[rawdf$Class_11 %in% year_end_keep,]
		yr1011 <- yr1011[which(yr1011$TotStLn_10 > 0), ]
		# survival to next census
		yr1011$Surv <- 1
		yr1011$Surv <- ifelse(is.na(yr1011$TotStLn_11), 0, 1)
	# redo for annual transition from 2011 - 2012
		yr1112  <- rawdf[rawdf$Class_11 %in% start_keep,]
		yr1112  <- rawdf[rawdf$Class_12 %in% year_end_keep,]
		yr1112 <- yr1112[which(yr1112$TotStLn_11 > 0), ]
		# survival to next census
		yr1112$Surv <- 1
		yr1112$Surv <- ifelse(is.na(yr1112$TotStLn_12), 0, 1)
	# rename variables to match conventional IPM notation in Reese et al 2014 & Merow papers
		colnames(yr1011)[which(names(yr1011) == "TotStLn_10")] <- "z" # initial size 
		colnames(yr1011)[which(names(yr1011) == "TotStLn_11")] <- "z1" # final size 
	# rename variables to match conventional IPM notation in Reese et al 2014 & Merow papers
		colnames(yr1112)[which(names(yr1112) == "TotStLn_11")] <- "z" # initial size 
		colnames(yr1112)[which(names(yr1112) == "TotStLn_12")] <- "z1" # final size 

# Will assign reproductive variables to the first year. So flowering or fecundity 
# will be regressed against size of the plant at the time of census. e.g. flowers in 
# 2011 on a plant are associated with the plant size in 2011 (z, time). This means that for year 1112 
# our flower count is from the year 2011 z (or time 1) and for year 1011 our flower count is from year
# 2010. The problem with this approach is that we are missing data for 2012. We will have to 
# add 2012 data onto the dataframe for scripts "04.2.3_Flower.R" and "04.2.4_Fecundity.R"

# binary (flower yes/no) rename to "Repr"
		yr1011$Repr[ is.na(yr1011$z)] <- NA  # it didn't survive to the census, or it was from another time period
		yr1011$Repr[ yr1011$z > 0] <- 0  # it survived, but was a juvenile and no flowers 
		yr1011$Repr[ yr1011$z > 0 & yr1011$TotFlStLn_10 > 0] <- 1 # it survived, and recorded flowers
		tail(yr1011[,c("Repr", "TotFlStLn_10", "TotNFStLn_10", "z", "z1")], 10) # good check 
	
# binary (flower yes/no) rename to "Repr"
		yr1112$Repr[ is.na(yr1112$z)] <- NA  # it didn't survive to the census, or it was from another time period
		yr1112$Repr[ yr1112$z > 0] <- 0  # it survived, but was a juvenile and no flowers 
		yr1112$Repr[ yr1112$z > 0 & yr1112$TotFlStLn_11 > 0] <- 1 # it survived, and recorded flowers
		tail(yr1112[,c("Repr", "TotFlStLn_11", "TotNFStLn_11", "z", "z1")], 10) # good check 
		summary(yr1112$Repr); summary(yr1011$Repr)
		
	# fruits production	
		colnames(yr1011)[which(names(yr1011) == "TotFr_10")] <- "Fec" # seed output 
		colnames(yr1112)[which(names(yr1112) == "TotFr_11")] <- "Fec" # seed output 
		
	# make another data frame for the 2012/2013 data (don't have 2013 data so will just use 2012 data 
	# for flowering & fecundity estimates. 
		yr1213 <- rawdf[rawdf$Class_12 %in% start_keep,]
		yr1213 <- yr1213[which(yr1213$TotStLn_12 > 0), ]
		yr1213$Surv <- NA
		colnames(yr1213)[which(names(yr1213) == "TotStLn_12")] <- "z"
		yr1213$z1 <- NA	
		yr1213$Repr[ is.na(yr1213$z)] <- NA  # it didn't survive to the census, or it was from another time period
		yr1213$Repr[ yr1213$z > 0] <- 0  # it survived, but was a juvenile and no flowers 
		yr1213$Repr[ yr1213$z > 0 & yr1213$TotFlStLn_12 > 0] <- 1 # it survived, and recorded flowers
		tail(yr1213[,c("Repr", "TotFlStLn_12", "TotNFStLn_12", "z", "z1")], 10) # good check 
		colnames(yr1213)[which(names(yr1213) == "TotFr_12")] <- "Fec" # seed output 
	
	# main objects 
		# Surv, z, z1, Repr, Fec, 
	
	# Missing seedling size distribution, will need to calculate this from existing plot. Unfortunately we
	# only have data for local recruitment. The seedling size distribution in a plot will be all new individual in 
	# the following year based on the number of seeds dropped in the previous year. For 2010 - 2011, recruitment is all 
	# new individuals recorded in plots in 2011 (at the site level?) based on the total seeds dropped by all individuals 
	# the previous year (at the site level?). 
	# Will come back to this in another script ("03_RecruitmentModelComp.R"), but first should do some exploratory analysis with the data. 
	
	#temporary merge of dfs
	yr1011sub <- yr1011[,c("Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_11", "NewPlot_12")]
	yr1112sub <- yr1112[,c("Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_11", "NewPlot_12")]
	yr1213sub <- yr1213[,c("Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_11", "NewPlot_12")]

	yr1011sub$Year <- "1011"
	yr1112sub$Year <- "1112"
	yr1213sub$Year <- "1213"
	
	explor <- rbind(yr1011sub, yr1112sub) # but don't bind on the 1213 data here
	
# END: part 1: Set directories, load libraries & import data




#============================================================================================#
# part 2: Check for outliers & weird points
# Outliers - boxplots & cleveland dotplots, brief check, but more rigorous check will be preformed during model diagnostics
#============================================================================================#

library(car)
	## save plots
	setwd(path.fig)
	pdf(file="01_OutlierPlots.pdf", width=11, height=8.5)
	par(mfrow=c(1,3))
		Boxplot(explor$z, explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
			xlab="Region", ylab="Initial Size")
		Boxplot(explor$z1, explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
			xlab="Region", ylab="Final size")
		Boxplot(explor$Fec, explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
			xlab="Region", ylab="Fecundity")
	# super giants difficult to deal with, but excluding too many outliers could take away from the actual 
	# potentially strong ecological influence of these in wild populaitons.  
			par(mfrow=c(1,3))
				Boxplot(log(explor$z+1), explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
					xlab="Region", ylab="LOG - Initial Size")
				Boxplot(log(explor$z1+1), explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
					xlab="Region", ylab="LOG - Final size")
				Boxplot(log(explor$Fec+1), explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
					xlab="Region", ylab="LOG - Fecundity")
	# growth increment
		par(mfrow=c(1,2))
		library(car)
		Boxplot((explor$z1)-(explor$z), explor$Region, labels=paste0(explor$ID, "-", explor$Year), data=explor, id.method="y",
				xlab="Region", ylab="Growth Increment")
				# notice plant ID14655 decreases one year & the next increases greatly 
		Boxplot((log(explor$z1)-log(explor$z)), explor$Region, labels=explor$ID, data=explor, id.method="y",
				xlab="Region", ylab="Log - Growth Increment")
		#**** major single growth increment outlier 14655
	# by site ID
		# but subset off 14655 to see rest of values 
				par(mfrow=c(1,1))
				drop <- explor[-which(explor$ID==14655), ]
				Boxplot((drop$z1)-(drop$z), drop$SiteID, labels=paste0(drop$ID, "-", drop$Year), data=drop, id.method="y",
				ylab="Growth Increment, without 14655",  xlab = "", main="Growth Increment, without 14655")
				Boxplot(drop$Fec, drop$SiteID, labels=paste0(drop$ID, "-", drop$Year), data=drop, id.method="y",
				ylab="Fecundity without 14655",  xlab = "", main="Fecundity without 14655")
			# In addition to 14655, try subsetting off a few more extreme outliers.
				par(mfrow=c(1,1))
				drop <- drop[-which(drop$ID==12343 | drop$ID==14586 | drop$ID==12397 | drop$ID==14572), ]
				Boxplot((drop$z1)-(drop$z), drop$SiteID, labels=paste0(drop$ID, "-", drop$Year), data=drop, id.method="y",
				ylab="Growth Increment, WITHOUT ALL MAJOR OUTLIERS",  xlab = "", main="Growth Increment, WITHOUT ALL MAJOR OUTLIERS")
				Boxplot(drop$Fec, drop$SiteID, labels=paste0(drop$ID, "-", drop$Year), data=drop, id.method="y",
				ylab="Fecundity WITHOUT ALL MAJOR OUTLIERS",  xlab = "", main="Fecundity WITHOUT ALL MAJOR OUTLIERS")
			# outliers somewhat restricted to limited number of sites 
		
	
	# Cleveland box plots 	
		par(mfrow=c(1,2))
		drop <- drop[order(drop$z1),]
		dotchart(drop$z, col=as.factor(drop$PlotID), pch=19, lcolor = NULL, xlab="Initial Size z", 
		ylab="Data ordered by final size (z1)", main="Cleveland's Dot Plots: col=PlotID")
	
		drop <- drop[order(drop$z),]
		dotchart(drop$z1, col=as.factor(drop$PlotID), pch=19, lcolor = NULL, xlab="final Size z1", 
		ylab="Data ordered by initial size (z)", main="Cleveland's Dot Plots: col=PlotID")
	
		par(mfrow=c(1,2))
		drop <- drop[order(drop$z1),]
		dotchart(drop$Fec, col=as.factor(drop$PlotID), pch=19, lcolor = NULL, xlab="Fec - fecundity", 
		ylab="Data ordered by final size (z1)", main="Fec Cleveland's Dot Plots: col=PlotID")
	
		drop <- drop[order(drop$z),]
		dotchart(drop$Fec, col=as.factor(drop$PlotID), pch=19, lcolor = NULL, xlab="Fec - fecundity", 
		ylab="Data ordered by initial size (z)", main="Fec Cleveland's Dot Plots: col=PlotID")
	
	dev.off(); dev.off()
	# outliers are important for population ecology
	# but might want to make a cutoff point for super giants? that only occur in about 1 or two sties
	

	
#============================================================================================#
# part 3: brief check for homogeneity of variance
# Homogeneity - conditional boxplots (again brief & preliminary)
#============================================================================================#
	## save plots
	setwd(path.fig)
	pdf(file="01_ScatterHist_state.pdf", width=11, height=8.5)
		# growth scatter plot
			par(mfrow=c(1,2))
			plot(drop$z, drop$z1, xlab="Initial Size", ylab="Final Size", bg=as.factor(drop$SiteID), pch=21, main="basic growth")
			abline(a = 0, b = 1, col="red", lty=2)
			# log size
				plot(log(drop$z), log(drop$z1), xlab="Initial Size", ylab="Final Size", bg=as.factor(drop$SiteID), pch=21, main="Log growth")
				abline(a = 0, b = 1, col="red", lty=2)
				
		# fecundity
				par(mfrow=c(1,2))
				plot(log(drop$z1), log(drop$Fec), xlab="Log - final size", ylab="Log - fecundity", bg=as.factor(drop$SiteID), pch=21, main="log(fecundity) given final size log(z1)")
				plot(log(drop$z1), log(drop$Fec), xlab="Log - final size", ylab="Log - fecundity", bg=as.factor(drop$Region), pch=21, main="log(fecundity) given final size log(z1)")

		# quick super basic models to look at residuals 
			g_mod <- lm(drop$z1 ~ drop$z, data=drop)	
			par(mfrow=c(2,2)); plot(g_mod, main="Growth")
			fec_mod <- glm(drop$Fec ~ drop$z1, data=drop, family = "poisson")	
			par(mfrow=c(2,2)); plot(fec_mod, main="Fecundity")

	
	
#============================================================================================#
# part 4: brief check for normality 
# Normality - histograms 
#============================================================================================#
		
			# include with previous
			par(mfrow=c(3,2))
			hist(drop$z, main="raw", col="lightgrey")
			hist(log(drop$z), main="Log trans", col="darkgrey")
			hist(drop$z1, main="raw", col="lightgrey")
			hist(log(drop$z1), main="Log trans", col="darkgrey")
			hist(drop$Fec, main="raw", col="lightgrey")
			hist(log(drop$Fec), main="Log trans", col="darkgrey")
			dev.off(); dev.off()

			
###################################################################
# Outlier plants dropped from the analysis
	# how many outliers were dropped 
		length(levels(as.factor(explor$ID))) - length(levels(as.factor(drop$ID)))		
		# dropped 
		# ID==14655, 12343, 14586, 12397, 14572
	# need to drop these from the 12-13 year yr1213
	yr1213sub <- yr1213sub[-which(yr1213sub$ID==14655 | yr1213sub$ID==12343 | yr1213sub$ID==14586 | yr1213sub$ID==12397 | yr1213sub$ID==14572), ]
	par(mfrow=c(1,2)); hist(log(yr1213sub$z+1)); plot(log(yr1213sub$z+1), log(yr1213sub$Fec+1))

	
#============================================================================================#
# part 5: save a refined data frame with some MAJOR outliers excluded 
# 	SAVE FILE  
#============================================================================================#
		# outliers excluded: 5 super giants IDs= 14655, 12343, 14586, 12397, 14572
		setwd(path.dat)
		write.csv(drop, file="demo_noMajOutlie.csv", row.names = FALSE)
		write.csv(yr1213sub, file="demo_addonFecRepr2012Data.csv", row.names = FALSE)

		