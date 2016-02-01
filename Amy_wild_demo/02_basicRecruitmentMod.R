#============================================================================================#
# ANALYSIS OF WILD CARDINALIS TRANSPLANT DATA 
#============================================================================================#
# April 24, 2015
# SCRIPT 2: Exploratory check of potentially suitable fecundity model 
# file name: 02_basicFecundity.R



# Build a basic fecundity model suitable for both Amy's Demography Dataset & Matt's Lat transplant dataset
# IPM P Kernel = growth kernel (P) + fecundity kernel (F)
# growth kernel P(z, z1) = s(z)G(z,z')

# Fecundity Kernel Layout:
	# In order to reproduce and add new plants to the population parents must survive and
	# probability of survival is a direct function of size s(z). Then individuals must become 
	# reproductive & start investing some energy towards reproduction. This will be a binary 
	# response of reproductive / not reproductive as a function of plant size flr(z'). Notice that 
	# we use z' here because the final size at the second census will be a better predictor of whether or
	# not the plant is flowering. Next the individual must produce fruit and seeds, a process 
	# that is also size dependant fec(z'). Unfortunately we have to assume for this project that
	# all fruits have the same number of seeds. A very small portion of these seeds will disperse to safe 
	# sites, survive the winter, germinate, emerge and survive over the summer until the next census (area of uncertainty). 
	# At the time of the next census, these new individuals will have been growing over the summer 
	# and their size will be a pdf of  C(0,z). 
		# Formula: s(z)Flr(z')Fec(z')[successful recruits]C(0,z)
		
# This script will work towards figuring out the best solution for the "successful recruits - area of uncertainty" given the data at hand. 
	
	# Basic idea of possible models (options?): 
		# A. successful recruits = a constant global portion (# new individuals / # of fruits dropped)
				# most basic - too simplistic? could be averaged across sites or site specific. 
		# B. successful recruits = number fruits in plot + (1|year/region/site/plot)
				# mixed model with random intercept 
		# C. successful recruits = number fruits in plot + (number fruits in plot|year/region/site/plot)
				# mixed model with random intercept & slope
		# D. successful recruits = number fruits in SITE + (number fruits in plot|year/region/site/plot)
				# same as above, but maybe number of fruits in site is a better predictor. 

library(lme4)
library(nlme)
library(glmmML)

				
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
	path.funct = paste(path.root,"/Rcode/Amy_wild_demo/functions", sep="")
	setwd(path.dat); dir()

	
# load in 2010-2012 data (two time step observations) 
# major outliers removed. 
dfclean <- read.csv(file="demo_noMajOutlie.csv")

summary(dfclean)
table(dfclean$SiteID)
table(dfclean$Region)

	# start a pdf to save all plots to 
	setwd(path.fig)
	pdf(file="02_Fecundity_Explor.pdf", width=11, height=8.5)

setwd(path.dat)
#============================================================================================#
# part A: Explore possibility of using a simple fixed constant based on number of fruits dropped
#============================================================================================#

# make variable for new recruit in the t+1 census (0/1 binary)
		dfclean$New[ is.na(dfclean$z) & is.na(dfclean$z1) ] <- NA # keep na if was NA both years 
		dfclean$New[ dfclean$z > 0 & dfclean$z1 > 0 ] <- 0 # it was a survivor, but not a new individual
		dfclean$New[ dfclean$z > 0 & is.na(dfclean$z1)] <- 0 # it died & did not make it to the next year
		dfclean$New[ is.na(dfclean$z) & dfclean$z1 > 0] <- 1 # it was newly found in the second census, but not in the first 
			# a new individual, a new recruit?
		# New individuals in the second year were stripped off the dataframe accidentaly in the previous script will add them bakc on here. 
setwd(path.dat)
getNew <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
	# New individuals: 
		# status in first year = NA
		# status in second year = J or A 
		# New plot in second year = FALSE
		new1011 <- getNew[which(is.na(getNew$Class_10) & getNew$NewPlot_11=="FALSE" & getNew$TotStLn_11 > 0),]
		new1112 <- getNew[which(is.na(getNew$Class_11) & getNew$NewPlot_12=="FALSE" & getNew$TotStLn_12 > 0),]

		# set variables to match existing data
		new1011$Surv <- NA; new1112$Surv <- NA 
		new1011$z <- NA; new1112$z <- NA
		new1011$z1 <- new1011$TotStLn_11; new1112$z1 <- new1112$TotStLn_12
		new1011$Repr <- 0; new1112$Repr <- 0;
		new1011$Repr[new1011$TotFlStLn_11 > 0] <- 1; new1112$Repr[new1112$TotFlStLn_12 > 0] <- 1
		new1011$Fec <- 0; new1112$Fec <- 0;
		new1011$Fec[new1011$TotFr_11 > 0] <- 1; new1112$Fec[new1112$TotFr_12 > 0] <- 1
		new1011$Year <- 1011; new1112$Year <- 1112
		new1011$New <- 1; new1112$New <- 1
		new1011 <- new1011[,c("Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_11", "NewPlot_12", "Year", "New")]
		new1112 <- new1112[,c("Surv", "z", "z1", "Repr", "Fec", "SiteID", "PlotID", "Region", "ID", "NewPlot_11", "NewPlot_12", "Year", "New")]

	# quick checks before merge 
		# summary(new1011); dim(1011)
		# summary(new1112); dim(1112)
		par(mfrow=c(1,2))
		hist(new1011$z1); hist(log(new1011$z1+1))
		hist(new1112$z1); hist(log(new1112$z1+1)); 
		
		# some box plots
		library(car)
			
			Boxplot(new1011$z1, new1011$Region, labels=paste0(new1011$ID, "-", new1011$Year),
				data=new1011, id.method="y", xlab="Region", ylab="size", main="Recruit Outliers in 11")
			Boxplot(log(new1011$z1 + 1), new1011$Region, labels=paste0(new1011$ID, "-", new1011$Year),
				data=new1011, id.method="y", xlab="Region", ylab="size", main="Recruit Outliers in 11")
			
			Boxplot(new1112$z1, new1112$Region, labels=paste0(new1112$ID, "-", new1112$Year),
				data=new1112, id.method="y", xlab="Region", ylab="size", main="Recruit Outliers in 12")
			Boxplot(log(new1112$z1 + 1), new1112$Region, labels=paste0(new1112$ID, "-", new1112$Year),
				data=new1112, id.method="y", xlab="Region", ylab="size", main="Recruit Outliers in 12")
		
		# 	QUICK A FEW IMPOSSIBLY LARGE RECRUITS - WILL NEED TO INVESTIGATE THESE POTENTIALLY 
		#	ERRONEOUS ID'S, MIGHT BE PLOT BOUNDARY MOVEMENT RATHER THAN RECRUIT. 

################################################
		# bind df together 
		dim(new1011); dim(new1112); dim(dfclean)
		newRecruits <- rbind(new1011, new1112)
		dfclean <- rbind(dfclean, newRecruits)
		dim(dfclean)

###############################################
# Is recruitment at the plot, site or region level related to seed rain recorded in the demo study?

	table(dfclean$New)
	# summarize data
	library(plyr)
	
	# have to add on fruit from 2010 for 2011 seedlings 
		setwd(path.dat)
		getNew <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
		add <- getNew[which(getNew$TotFr_10 > 0), ]
	
	# for 2011 
		new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
		test1_11recs <- ddply(new1011,~PlotID,summarise,recruits=sum(New))
		test1_11fruits <- ddply(add,~PlotID,summarise,fruits=sum(TotFr_10))
		test1_11 <- merge(test1_11fruits, test1_11recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_11[is.na(test1_11)] <- 0
	
	# for 2012 
		new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
		test1_12recs <- ddply(new1112,~PlotID,summarise,recruits=sum(New))
		add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
		test1_12fruits <- ddply(add,~PlotID,summarise,fruits=sum(Fec))
		test1_12 <- merge(test1_12fruits, test1_12recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_12[is.na(test1_12)] <- 0
			
			
	# fruits last year in plot & recruits this year in plot 2011
		par(mfrow=c(1,3))
		plot(log(test1_11$fruits+1), test1_11$recruits, xlab="Log fruits in plot previous year", 
		ylab="Recruits in plot for current year", bg="green", main="PlotID 2011-green & 2012-blue", pch=21)
		points(log(test1_12$fruits+1), test1_12$recruits, bg="blue", pch=21)
	#********* pretty much no clear relationship at all!!!!
	# maybe fruits dropped at site in the previous year might be a better predictor
	# For SITES
			# for 2011 
				new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
				test2_11recs <- ddply(new1011,~SiteID,summarise,recruits=sum(New))
				add <- getNew[which(getNew$TotFr_10 > 0), ]
				test2_11fruits <- ddply(add,~SiteID,summarise,fruits=sum(TotFr_10))
				test2_11 <- merge(test2_11fruits, test2_11recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_11[is.na(test2_11)] <- 0
			# for 2012 
				new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
				test2_12recs <- ddply(new1112,~SiteID,summarise,recruits=sum(New))
				add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
				test2_12fruits <- ddply(add,~SiteID,summarise,fruits=sum(Fec))
				test2_12 <- merge(test2_12fruits, test2_12recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_12[is.na(test2_12)] <- 0
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test2_11$fruits+1), test2_11$recruits, xlab="Log fruits in SITE previous year", 
		ylab="Recruits in SITE for current year", bg="green", main="SiteID 2011-green & 2012-blue", pch=21)
		points(log(test2_12$fruits+1), test2_12$recruits, bg="blue", pch=21)
	
	# MAYBE EVEN REGIONS???
					# for 2011 
					new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
					test3_11recs <- ddply(new1011,~Region,summarise,recruits=sum(New))
					add <- getNew[which(getNew$TotFr_10 > 0), ]
					test3_11fruits <- ddply(add,~Region,summarise,fruits=sum(TotFr_10))
					test3_11 <- merge(test3_11fruits, test3_11recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_11[is.na(test3_11)] <- 0
				# for 2012 
					new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
					test3_12recs <- ddply(new1112,~Region,summarise,recruits=sum(New))
					add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
					test3_12fruits <- ddply(add,~Region,summarise,fruits=sum(Fec))
					test3_12 <- merge(test3_12fruits, test3_12recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_12[is.na(test3_12)] <- 0
		# fruits last year in plot & recruits this year in plot 2011
			plot(log(test3_11$fruits+1), test3_11$recruits, xlab="Log fruits in REGION previous year", 
			ylab="Recruits in REGION for current year", bg="green", main="Region 2011-green & 2012-blue", pch=21)
			points(log(test3_12$fruits+1), test3_12$recruits, bg="blue", pch=21)
		

###############################################
# is recruitment consistent across the three scales 	
		par(mfrow=c(1,3))
	# at plot level
		temp_plot <- merge(test1_11, test1_12, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		plot(temp_plot$recruits.x, temp_plot$recruits.y, xlab="2011 recruits", ylab="2012 recruits", main="Site LEVEL", pch=21, bg="darkgrey")

	# at site level
		temp_sites <- merge(test2_11, test2_12, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
		plot(temp_sites$recruits.x, temp_sites$recruits.y, xlab="2011 recruits", ylab="2012 recruits", main="Site LEVEL", pch=21, bg="darkgrey")

	# at region level
			plot(test3_11$recruits, test3_12$recruits, xlab="2011 recruits", ylab="2012 recruits", main="REGION LEVEL", pch=21, bg="darkgrey")
			text(test3_11$recruits, test3_12$recruits,labels=test3_11$Region)
			
###############################################
# is FECUNDITY consistent across the three scales FOR EACH YEAR
		par(mfrow=c(1,3))	
	# at plot level
		temp_plot <- merge(test1_11, test1_12, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		plot(log(temp_plot$fruits.x), log(temp_plot$fruits.y), xlab="2011 fruits", ylab="2012 fruits", main="Site LEVEL", pch=21, bg="darkred")

	# at site level
		temp_sites <- merge(test2_11, test2_12, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
		plot(log(temp_sites$fruits.x), log(temp_sites$fruits.y), xlab="2011 fruits", ylab="2012 fruits", main="Site LEVEL", pch=21, bg="darkred")

	# at region level
			plot(log(test3_11$fruits), log(test3_12$fruits), xlab="2011 fruits", ylab="2012 fruits", main="REGION LEVEL", pch=21, bg="darkred")
			text(log(test3_11$fruits), log(test3_12$fruits),labels=test3_11$Region)
		
	
# close plots for script



#============================================================================================#
# part B: Possible ideas for the discrepancy & heavy noise between number of fruits
# in the previous year & the number of recruits in the current year 
#============================================================================================#
# - 1. Wandering population - strong seed import and export out of plots inhibits our ability to 
							# detect the relationship?
# - 2. Heavy noise - perhaps largely determined by stochastic chaotic processes & deterministic process 
							# that we can explain with our data is ther but very very small ~ XXXL sample 
							# size would be needed. 
# - 3. Recruits from current years seed - if seeds in the current year are dropped in July recruitment happens
							# & individuals are recorded, but in reality these individuals have 
							# no real connection to last years fecundity. There XXXS size means they will 
							# probably not survive the winter & are really just temporay byproducts of the
							# fruit - recruit continuum?
# - 4. Asexual reproduction (colonality) - perhaps reproduction from seed only accounts for ~ 1/2 of new recruits
							# skewing the relationship in the dataset?
# - 5. seeds/fruit is more ecologically relevant - if # of fruit produced is only part of the picture
							# seeds per fruit highly variable within and between sites & is very relevant?
										
		# there is a weak positive relationship with the amount of fruits at the SITE level 
		# in the previous year & recruits in the subsequent year. This would make the most sense if 
		# the "wandering populations" effect is really strong. 
		
		# Sould we exclude XXXS seedlings from the rectuit count as they might be recruits from the 
		# current year rather than recruits from the previous year. 
		
			# use histograms to see distribution if we exclude individuals under a critical size 
				# arbitrary split at 10cm for now
				par(mfrow=c(2,2))
				check1112 <- new1112[which(new1112$z1 < 10), ]
				hist(check1112$z1)
				check1112 <- new1112[which(new1112$z1 > 10 & new1112$z1 < 100), ]
				hist(check1112$z1)
				check1011 <- new1011[which(new1011$z1 < 10), ]
				hist(check1011$z1)
				check1011 <- new1011[which(new1011$z1 > 10 & new1011$z1 < 100), ]
				hist(check1011$z1)
					# *** here they are? the current year recruits?
		
			# what about the XXL recruits, could they be from fragments, typos, moved plots boundaries ect. 
			check1011 <- new1011[which(new1011$z1 > 100), ]
			hist(check1011$z1, breaks = 20)
			# huge drop off after about 200cm, but need to figure out whats going on?
			check1112 <- new1112[which(new1112$z1 > 100), ]
			hist(check1112$z1) # 5000 cm ? tsl no way! 
			
		dev.off()
	
		# Should check to see what's the better predictor of recruits in the current year. 
		# Is it the number of FRUITS in the previous year? or the number of RECRUITS in the previous year
		# Would hopefully expect it to be the number of fruits in the previous year. If number of recruits 
		# in previous year turns out to be the better predictor, then were hooped for make the 
		# fruits-to-recruits connection 
		
			# close any open plotting window
		dev.off()
		# source script called extra plots due to long repetitive code. 
		setwd(path.funct); dir()
		
	# Run some additional plots to explore the relationship between fruit output & recruits 
		#    source("extraFecPlots.R") # make a big figure 
		
	
	
	
#============================================================================================#
# part C: Distribution fitting: Understand size-recruits relationship with distribution fitting. 
#============================================================================================#


	setwd(path.fig)
	pdf(file="02_Fecundity_Explor3_Distributions.pdf", width=11, height=8.5)

		# What is the probability distribution of recruits, what type of curve best describes the 
		# data & what should we be aiming for with our model?
			# Explore recruits size distribution by fitting different curves to data
			# Probability distribution:
				# We know distribution of recruits is a continuous variables with an asymmetrical distribution 
				# outliers are mostly positive. 
					# Possible distributions:
						# - exponential - can't have an infinite number of XXXS individuals 
						# - log normal - 
						# - Gamma
						# - Weibul 
		
		library(fitdistrplus)
		
		# take a quick look at the breaks again
			par(mfrow=c(4,2))
			check1011 <- new1011[which(new1011$z1 > 10 & new1011$z1 < 300), ]
			hist(check1011$z1); hist(log(check1011$z1+1))  # follows an exponential distribution
			check1011 <- new1011[which(new1011$z1 > 0 & new1011$z1 < 5), ]
			hist(check1011$z1); hist(log(check1011$z1+1))  # follows an exponential distribution
		
			check1112 <- new1112[which(new1112$z1 > 10 & new1112$z1 < 300), ]
			hist(check1112$z1); hist(log(check1112$z1+1))  # follows an exponential distribution
			check1112 <- new1112[which(new1112$z1 > 0 & new1112$z1 < 5), ]
			hist(check1112$z1); hist(log(check1112$z1+1))  # follows an exponential distribution
		
		# really looks like an exponential distribution, but cant have infinite number of XXXXXS plants so try 
		# log normal distribution. 
			year1011 <- new1011[which(new1011$z1 < 300), ]
			year1011 <- as.numeric(year1011$z1) 
			year1112 <- new1112[which(new1112$z1 < 300), ]
			year1112 <- as.numeric(year1112$z1) 
			
			
# CONVERT TO LOG SCALE TO FIT WITH REST OF FUNCTIONS 
	year1011 <- log(year1011+1)
	year1112 <- log(year1112+1)
	new1011$z1 <- log(new1011$z1 + 1)
	new1112$z1 <- log(new1112$z1 + 1)

		
			# basic check 
			par(mfrow=c(2,2))
			plotdist(year1011)
			plotdist(year1112)

			# skewness kurtosis plots 
			par(mfrow=c(1,2))
			descdist(year1011, boot = 1000)
			descdist(year1112, boot = 1000)
			# Weibull, gamma and lognormal distributions should be considered 
		
			# fit distributions 
			par(mfrow=c(2,2))
			fw <- fitdist(year1011, "weibull")
			fg <- fitdist(year1011, "gamma")
			fln <- fitdist(year1011, "lnorm")
			fn <- fitdist(year1011, "norm")

			
			# plot distributions 
			plot.legend <- c("Weibull", "lognormal", "gamma", "norm")
			denscomp(list(fw, fln, fg, fn), legendtext = plot.legend)
			qqcomp(list(fw, fln, fg, fn), legendtext = plot.legend)
			cdfcomp(list(fw, fln, fg, fn), legendtext = plot.legend)
			ppcomp(list(fw, fln, fg, fn), legendtext = plot.legend)
			summary(fw)
			summary(fg)
			summary(fln)
			summary(fn)
	
			# try exponential
			fexp <- fitdist(year1011, "exp")
			summary(fexp)
	
			# log normal seems to be the winner for 2010 - 2011, check 2012 data
							fw <- fitdist(year1112, "weibull")
							fg <- fitdist(year1112, "gamma")
							fln <- fitdist(year1112, "lnorm")
							fn <- fitdist(year1112, "norm")

							summary(fw)
							summary(fg)
							summary(fln)
							summary(fn)
							# again log normal here too. 
			
			fln11 <- fitdist(year1011, "lnorm")
			fln12<- fitdist(year1112, "lnorm")
			summary(fln11)
			summary(fln12)
			descdist(year1011, boot = 1000)
			descdist(year1112, boot = 1000)	
		
			# bootstrap estimates
			fln11.B <- bootdist(fln11, niter = 1001)
			fln12.B <- bootdist(fln12, niter = 1001)
			summary(fln11.B)
			summary(fln12.B)
			plot(fln11.B); 
	
			# random sample of distributions 
			samp <- rlnorm(1000, meanlog = fln12$estimate[1], sdlog = fln12$estimate[2])
			par(mfrow=c(1,1))

				
			# Investigate differences across sites & years: 
					my_table11 <- data.frame()
					Region <- levels(new1011$Region)
					for(i in 1:length(Region)){
						 tempy <- Region[i]
						just_reg <- new1011[which(new1011$Region==tempy), ]
						year1011 <- just_reg[which(just_reg$z1 < 300), ]
						year1011 <- as.numeric(year1011$z1) 
						fln11 <- fitdist(year1011, "lnorm")
						samp <- rlnorm(1000, meanlog = fln11$estimate[1], sdlog = fln11$estimate[2])
						d <- density(samp)
						assign(paste("d", Region[i], sep=""), d)
						meanlog <- fln11$estimate[1]; sdlog <- fln11$estimate[2]
						temper <- data.frame(tempy, meanlog, sdlog)
						my_table11 <- rbind(my_table11, temper)
					}
					my_table11
					plot(dS1, col="red", xlim=c(0, 10), main="2011")
					lines(dS2, col="red")
					lines(dC1, col="green"); lines(dC2, col="green"); lines(dC3, col="green");
					lines(dN1, col="blue"); lines(dN2, col="blue")
			
			# Investigate differences across sites & years: 
					my_table12 <- data.frame()
					Region <- levels(new1112$Region)
					for(i in 1:length(Region)){
						 tempy <- Region[i]
						just_reg <- new1112[which(new1112$Region==tempy), ]
						year1112 <- just_reg[which(just_reg$z1 < 300), ]
						year1112 <- as.numeric(year1112$z1) 
						fln11 <- fitdist(year1112, "lnorm")
						samp <- rlnorm(1000, meanlog = fln11$estimate[1], sdlog = fln11$estimate[2])
						d <- density(samp)
						assign(paste("d", Region[i], sep=""), d)
						meanlog <- fln11$estimate[1]; sdlog <- fln11$estimate[2]
						temper <- data.frame(tempy, meanlog, sdlog)
						my_table12 <- rbind(my_table12, temper)
					}
					plot(dS1, col="red", xlim=c(0, 10), main="2012")
					lines(dS2, col="red")
					lines(dC1, col="green"); lines(dC2, col="green"); lines(dC3, col="green");
					lines(dN1, col="blue"); lines(dN2, col="blue")
					my_table12
				dev.off()
		
				my_table12$year <- 2012
				my_table11$year <- 2011
				both_years <- rbind(my_table11, my_table12)
				setwd(path.obj)
				write.csv(both_years, file = "lnormFecKern.csv", row.names = FALSE)
					
		

		
#============================================================================================#
# part D: Explore with mixed effect model framework  
#============================================================================================#

# make dataframe to work with 
	# from plot level (add together both years)
		test1_11$Year <- 2011
		test1_12$Year <- 2012
		test1 <- rbind(test1_11, test1_12)
	# add on site & region ID
		head(dfclean)
		addSitesRegions <- dfclean # make a temporary df
		addSitesRegions$Unique <- paste0(addSitesRegions$PlotID, addSitesRegions$Year) # concatenate year & plot to a unique column 
		addSitesRegions <- addSitesRegions[, c("PlotID", "SiteID", "Region")] # subset columns of interest 
		
		recruits1 <- merge(test1, addSitesRegions, by.x="PlotID", by.y="PlotID", all.x=TRUE)
		dim(recruits1); dim(test1) # should be the same
		#recruits1$Unique <- paste0(recruits1$PlotID, recruits1$Year) # concatenate year & plot to a unique column 
		recruits2 <- recruits1[!duplicated(recruits1), ]
		rec <- recruits2 # rename for easy coding
		dim(rec)

	####################################################
	# begin MEM analysis 
	# will follow Zuur 2009: Mixed Effect Models and Extensions in Ecology
	# Chapter 19: Mixed Effects Modelling Applied on American Foulbrood Affecting Honey Bees Larvae
	# pages: 447 - 
			# plots are nested within sites nested within regions nested within year
			# YEAR/REGION/SITE/Plot
			# 19.2 Data exploration: 
				rec$PlotID <- factor(rec$PlotID) # turn plot ID to factor
				rec$Lrecruits <- log(rec$recruits + 1)
				rec$Lfruits <- log(rec$fruits + 1)
				op<- par(mfrow = c(2, 2), mar = c(3, 4, 1, 1))
				dotchart(rec$recruits, group = rec$PlotID)
				dotchart(rec$Lrecruits, group = rec$PlotID)
				dotchart(rec$fruits, group = rec$PlotID)
				dotchart(rec$Lfruits, group = rec$PlotID)
				rec$Year <- as.factor(rec$Year)
				# Comments: two plots have an absurdly high number of fruits. This will effect problems with homoeneity
				#One option is to use different variance per plot (GLS)
				#but this will results in dozens of extra variance estimates
				#For each sub model. So should use log transformation 
				# to deal with fruit outliers. 

				# should first remove the hyper extreme outliers 
				# around 10,000 & 15,000 fruits per plot
				rec <- rec[which(rec$fruits < 9000), ]
				rec <- rec[complete.cases(rec),] # drop NAs

				# Try re-grouping regions to just three levels, rather than 
				# six. Results identical if regions are grouped with 3 or 6 levels.
				# regional groupings still fail to enter final model. 
				rec$Region <- revalue(rec$Region, c("C1"="C", "C2"="C", "C3"="C", "N1"="N", "N2"="N", "S1"="S", "S2"="S"))

			# 19.3 Analysis of the data:
				# plot a simple linear model of recruits~fruits
					M1 <- lm(Lrecruits ~ Lfruits, data=rec)
					E1 <- rstandard(M1)
					par(mfrow=c(2,2))
					plot(E1 ~ rec$PlotID, xlab="Plot")
					plot(E1 ~ rec$SiteID, xlab="Site")
					plot(E1 ~ rec$Region, xlab="Region")
					plot(E1 ~ as.factor(rec$Year), xlab="Year")
				# Protocol for fitting MEM
					# Step 1: start with model that contains all possible explanatory variables
						# and their interactions. Does MEM out preform linear model?
						library(nlme)
						M2 <- gls(Lrecruits ~ Lfruits, data=rec) # have to use gls to make comparison
						M3 <- lme(Lrecruits ~ Lfruits, random=~1|PlotID, data=rec)
						anova(M2, M3)
							# not really pretty close for plots, but what about SiteID
						M3b <- lme(Lrecruits ~ Lfruits, random=~1|SiteID, data=rec, na.action=na.exclude)
						anova(M2, M3b)
							# yes, strong effect of site
								# but testing on the boundary, still significant?
								0.5 * (1 - pchisq(15.78221, 1)) # yes!
						M3c <- lme(Lrecruits ~ Lfruits, random=~1|Region, data=rec, na.action=na.exclude)
						anova(M2, M3c)
							# no, not for region or year either
					# Then explore site-level random effects
						# is a random slope also significant?
						M4 <- lme(Lrecruits ~ Lfruits, random=~1 + Lfruits|PlotID, data=rec, na.action=na.exclude)
						M4b <- lme(Lrecruits ~ Lfruits, random=~1 + Lfruits|SiteID, data=rec, na.action=na.exclude)
						anova(M2, M3b, M4b)
						# No, no need for random slope at site level or plot level 
						# residual plot 
						plot(M4b)
						plot(M4)
						# some major problems with the residuals for the random intercept only model 
						
				# Compare global nested model: 
					M5 <- lmer(Lrecruits ~ Lfruits + (1|Region/SiteID/PlotID),  na.action=na.omit, data=rec)
					M6 <- lmer(Lrecruits ~ Lfruits + (1|SiteID/PlotID),  na.action=na.omit, data=rec)
					M7 <- lmer(Lrecruits ~ Lfruits + (1|Region/PlotID),  na.action=na.omit, data=rec)
					anova(M5, M6)
					anova(M6, M7)
					anova(M5, M7)
				   M8 <- lmer(Lrecruits ~ Lfruits + (1|SiteID),  na.action=na.omit, data=rec)
					anova(M6, M8) # model 8 is winner
				  # Does regional does adding year as a random effect improve the model?
				   M9 <- lmer(Lrecruits ~ Lfruits + (1|SiteID) + (1|Year),  na.action=na.omit, data=rec)
					anova(M8, M9) # no 
					
				# regional groupings or year as fixed effects?
				  # Does regional does adding year as a random effect improve the model?
				   M10 <- lmer(Lrecruits ~ Lfruits + Year + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				   M11 <- lmer(Lrecruits ~ Lfruits + Region + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				   M12 <- lmer(Lrecruits ~ Lfruits + Year + Region + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
					anova(M10, M11)
					anova(M10, M12)
					anova(M11, M12)
				  M13 <- lmer(Lrecruits ~ Lfruits + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
					anova(M11, M13)
				 M14 <- lmer(Lrecruits ~ Region*Lfruits + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 anova(M11, M14)
				 M15 <- lmer(Lrecruits ~ Region*Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=TRUE)
				 anova(M14, M15)
				 M16 <- lmer(Lrecruits ~ Region + Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 M15 <- lmer(Lrecruits ~ Region*Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 anova(M15, M16)
				 anova(M11, M16)

# TOP MODEL AT PLOT LEVEL
Rec_model <- lmer(Lrecruits ~ Lfruits + Region + (Lfruits|SiteID),  na.action=na.omit, data=rec)
		
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################

# Run model at 'site level' (sum of recruits & fruits in site-year) to see if predictions change/improve?

# make dataframe to work with 
	# from plot level (add together both years)
		test2_11$Year <- 2011
		test2_12$Year <- 2012
		test2 <- rbind(test2_11, test2_12)
	# add on site & region ID
		head(dfclean)
		addSitesRegions <- dfclean # make a temporary df
		addSitesRegions$Unique <- paste0(addSitesRegions$PlotID, addSitesRegions$Year) # concatenate year & plot to a unique column 
		addSitesRegions <- addSitesRegions[, c("SiteID", "Region")] # subset columns of interest 
		
		recruits2 <- merge(test2, addSitesRegions, by.x="SiteID", by.y="SiteID", all.x=TRUE)
		recruits2 <- recruits2[!duplicated(recruits2), ]
		dim(recruits2); dim(test2) # should be the same
		rec <- recruits2 # rename for easy coding
		dim(rec)

				rec$Lrecruits <- log(rec$recruits + 1)
				rec$Lfruits <- log(rec$fruits + 1)
				op<- par(mfrow = c(2, 2), mar = c(3, 4, 1, 1))
				dotchart(rec$recruits, group = rec$SiteID)
				dotchart(rec$Lrecruits, group = rec$SiteID)
				dotchart(rec$fruits, group = rec$SiteID)
				dotchart(rec$Lfruits, group = rec$SiteID)
				rec$Year <- as.factor(rec$Year)
				# Comments: two plots have an absurdly high number of fruits. This will effect problems with homoeneity
				#One option is to use different variance per plot (GLS)
				#but this will results in dozens of extra variance estimates
				#For each sub model. So should use log transformation 
				# to deal with fruit outliers. 

				# should first remove the hyper extreme outliers 
				# around 10,000 & 15,000 fruits per plot
				rec <- rec[which(rec$fruits < 9000), ]
				rec <- rec[complete.cases(rec),] # drop NAs

				# Try re-grouping regions to just three levels, rather than 
				# six. Results identical if regions are grouped with 3 or 6 levels.
				# regional groupings still fail to enter final model. 
				rec$Region <- revalue(rec$Region, c("C1"="C", "C2"="C", "C3"="C", "N1"="N", "N2"="N", "S1"="S", "S2"="S"))

			# 19.3 Analysis of the data:
				# plot a simple linear model of recruits~fruits
					M1 <- lm(Lrecruits ~ Lfruits, data=rec)
					E1 <- rstandard(M1)
					par(mfrow=c(2,2))
					plot(E1 ~ rec$SiteID, xlab="Site")
					plot(E1 ~ rec$Region, xlab="Region")
					plot(E1 ~ as.factor(rec$Year), xlab="Year")
				# Protocol for fitting MEM
					# Step 1: start with model that contains all possible explanatory variables
						# and their interactions. Does MEM out preform linear model?
						library(nlme)
						M2 <- gls(Lrecruits ~ Lfruits, data=rec) # have to use gls to make comparison
							#  what about SiteID
						M3b <- lme(Lrecruits ~ Lfruits, random=~1|SiteID, data=rec, na.action=na.exclude)
						anova(M2, M3b)
							# yes, strong effect of site
								# but testing on the boundary, still significant?
								0.5 * (1 - pchisq(15.78221, 1)) # yes!
						M3c <- lme(Lrecruits ~ Lfruits, random=~1|Region, data=rec, na.action=na.exclude)
						anova(M2, M3c)
							# no, not for region or year either
					# Then explore site-level random effects
						# is a random slope also significant?
						M4 <- lme(Lrecruits ~ Lfruits, random=~1 + Lfruits|SiteID, data=rec, na.action=na.exclude)
						anova(M3b, M4)
						# No, no need for random slope at site level or plot level 
						# residual plot 
						plot(M4b)
						plot(M4)
						# some major problems with the residuals for the random intercept only model 
						
				# Compare global nested model: 
					M5 <- lmer(Lrecruits ~ Lfruits + (1|Region/SiteID),  na.action=na.omit, data=rec)
					M6 <- lmer(Lrecruits ~ Lfruits + (1|SiteID),  na.action=na.omit, data=rec)
					M7 <- lmer(Lrecruits ~ Lfruits + (1|Region),  na.action=na.omit, data=rec)
					anova(M5, M6)
					anova(M6, M7)
					anova(M5, M7)
				   M8 <- lmer(Lrecruits ~ Lfruits + (1|SiteID),  na.action=na.omit, data=rec)
				  # Does adding year as a random effect improve the model?
				   M9 <- lmer(Lrecruits ~ Lfruits + (1|SiteID) + (1|Year),  na.action=na.omit, data=rec)
					anova(M8, M9) # no 
					
				# regional groupings or year as fixed effects?
				  # Does regional does adding year as a random effect improve the model?
				   M10 <- lmer(Lrecruits ~ Lfruits + Year + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				   M11 <- lmer(Lrecruits ~ Lfruits + Region + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				   M12 <- lmer(Lrecruits ~ Lfruits + Year + Region + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
					anova(M10, M11)
					anova(M10, M12)
					anova(M11, M12)
				  M13 <- lmer(Lrecruits ~ Lfruits + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
					anova(M10, M13)
				 M14 <- lmer(Lrecruits ~ Region*Lfruits + (1|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 anova(M11, M14)
				 M15 <- lmer(Lrecruits ~ Region*Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=TRUE)
				 anova(M14, M15)
				 M16 <- lmer(Lrecruits ~ Region + Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 M15 <- lmer(Lrecruits ~ Region*Lfruits + (Lfruits|SiteID),  na.action=na.omit, data=rec, REML=FALSE)
				 anova(M15, M16)
				 anova(M11, M16)

# TOP MODEL AT SITE LEVEL		
	Top_site_model <- lmer(Lrecruits ~ Region*Lfruits + (1|SiteID),  na.action=na.omit, data=rec)
	
# save rec object
setwd(path.obj)
write.csv(rec, file="siteFec.csv", row.names=FALSE)

	
#============================================================================================#
# part E: Save model objects & relevant components
#============================================================================================#

		Rec_model <- lmer(Lrecruits ~ Region*Lfruits + (1|SiteID),  na.action=na.omit, data=rec)
		setwd(path.obj)
		save(Rec_model, file = "Rec_model_site.rda")
		
		##################################
		# model components 
		Rec_model
		slotNames(Rec_model)
		
		methods(class = "merMod")
		fixef(Rec_model)
		confint(Rec_model, level = 0.95)
		coef(Rec_model)
		fixef(Rec_model) / sigma(Rec_model)
		ranef(Rec_model)
		head(Rec_model@frame)
		####################################
		re1 <- ranef(Rec_model, condVar=TRUE) # save the ranef.mer object
		class(re1)
		# conditional variance 
		attr(re1[[1]], which = "postVar")
		re1 <- ranef(MLexamp1, condVar=TRUE, whichel = "SiteID")
		print(re1)		
		dotplot(re1)

	# Use simulation and plots to explore random effects. 
		# remember fruit counts and recruit counts are logged transformed. 
		# need to run with GLMM to see if results consistent - since actually count data. 
		
		
		
#########################################
#########################################
# save plot, site & region level summaries 
setwd(path.obj)

	# plot & region
	plot_fec <- rbind(test1_11, test1_12)
	reg_fec <- rbind(test3_11, test3_12)
	
	write.csv(plot_fec, file="plotFecRec.csv", row.names=FALSE)
	write.csv(reg_fec, file="regionFecRec.csv", row.names=FALSE)