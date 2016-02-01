	setwd(path.fig)
	pdf(file="02_Fecundity_Explor2.pdf", width=11, height=8.5)

###############################################
# TOTAL RECRUIT SIZE VS. NUMBER OF FRUITS DROPPED?

	table(dfclean$New)
	# summarize data
	library(plyr)
	
	# have to add on fruit from 2010 for 2011 seedlings 
		setwd(path.dat)
		getNew <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
		add <- getNew[which(getNew$TotFr_10 > 0), ]
	
	# for 2011 
		new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
		test1_11recs <- ddply(new1011,~PlotID,summarise,recruits=sum(z1))
		test1_11fruits <- ddply(add,~PlotID,summarise,fruits=sum(TotFr_10))
		test1_11 <- merge(test1_11fruits, test1_11recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_11[is.na(test1_11)] <- 0
	
	# for 2012 
		new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
		test1_12recs <- ddply(new1112,~PlotID,summarise,recruits=sum(z1))
		add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
		test1_12fruits <- ddply(add,~PlotID,summarise,fruits=sum(Fec))
		test1_12 <- merge(test1_12fruits, test1_12recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_12[is.na(test1_12)] <- 0
			
			
	# fruits last year in plot & recruits this year in plot 2011
		par(mfrow=c(3,3))
		plot(log(test1_11$fruits+1), log(test1_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log SUM all recruit size @ t", bg="green", main="PlotID 2011-green & 2012-blue", pch=21)
		points(log(test1_12$fruits+1), log(test1_12$recruits+1), bg="blue", pch=21)
	#********* pretty much no clear relationship at all!!!!
	# maybe fruits dropped at site in the previous year might be a better predictor
	# For SITES
			# for 2011 
				new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
				test2_11recs <- ddply(new1011,~SiteID,summarise,recruits=sum(z1))
				add <- getNew[which(getNew$TotFr_10 > 0), ]
				test2_11fruits <- ddply(add,~SiteID,summarise,fruits=sum(TotFr_10))
				test2_11 <- merge(test2_11fruits, test2_11recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_11[is.na(test2_11)] <- 0
			# for 2012 
				new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
				test2_12recs <- ddply(new1112,~SiteID,summarise,recruits=sum(z1))
				add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
				test2_12fruits <- ddply(add,~SiteID,summarise,fruits=sum(Fec))
				test2_12 <- merge(test2_12fruits, test2_12recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_12[is.na(test2_12)] <- 0
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test2_11$fruits+1), log(test2_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log SUM all recruit size @ t", bg="green", main="SiteID 2011-green & 2012-blue", pch=21)
		points(log(test2_12$fruits+1), log(test2_12$recruits+1), bg="blue", pch=21)
	
	# MAYBE EVEN REGIONS???
					# for 2011 
					new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
					test3_11recs <- ddply(new1011,~Region,summarise,recruits=sum(z1))
					add <- getNew[which(getNew$TotFr_10 > 0), ]
					test3_11fruits <- ddply(add,~Region,summarise,fruits=sum(TotFr_10))
					test3_11 <- merge(test3_11fruits, test3_11recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_11[is.na(test3_11)] <- 0
				# for 2012 
					new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
					test3_12recs <- ddply(new1112,~Region,summarise,recruits=sum(z1))
					add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
					test3_12fruits <- ddply(add,~Region,summarise,fruits=sum(Fec))
					test3_12 <- merge(test3_12fruits, test3_12recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_12[is.na(test3_12)] <- 0
		# fruits last year in plot & recruits this year in plot 2011
			plot(log(test3_11$fruits+1), log(test3_11$recruits+1), xlab="Log fruits @ t-1", 
			ylab="log SUM all recruit size @ t", bg="green", main="Region 2011-green & 2012-blue", pch=21)
			points(log(test3_12$fruits+1), log(test3_12$recruits+1), bg="blue", pch=21)
		
		

###############################################
###############################################
###############################################
###############################################
# MEAN RECRUIT SIZE VS. NUMBER OF FRUITS DROPPED?

	table(dfclean$New)
	# summarize data
	library(plyr)
	
	# have to add on fruit from 2010 for 2011 seedlings 
		setwd(path.dat)
		getNew <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
		add <- getNew[which(getNew$TotFr_10 > 0), ]
	
	# for 2011 
		new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
		test1_11recs <- ddply(new1011,~PlotID,summarise,recruits=mean(z1))
		test1_11fruits <- ddply(add,~PlotID,summarise,fruits=sum(TotFr_10))
		test1_11 <- merge(test1_11fruits, test1_11recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_11[is.na(test1_11)] <- 0
	
	# for 2012 
		new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
		test1_12recs <- ddply(new1112,~PlotID,summarise,recruits=mean(z1))
		add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
		test1_12fruits <- ddply(add,~PlotID,summarise,fruits=sum(Fec))
		test1_12 <- merge(test1_12fruits, test1_12recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_12[is.na(test1_12)] <- 0
			
			
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test1_11$fruits+1), log(test1_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log MEAN all recruit size @ t", bg="green", main="PlotID 2011-green & 2012-blue", pch=21)
		points(log(test1_12$fruits+1), log(test1_12$recruits+1), bg="blue", pch=21)
	#********* pretty much no clear relationship at all!!!!
	# maybe fruits dropped at site in the previous year might be a better predictor
	# For SITES
			# for 2011 
				new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
				test2_11recs <- ddply(new1011,~SiteID,summarise,recruits=mean(z1))
				add <- getNew[which(getNew$TotFr_10 > 0), ]
				test2_11fruits <- ddply(add,~SiteID,summarise,fruits=sum(TotFr_10))
				test2_11 <- merge(test2_11fruits, test2_11recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_11[is.na(test2_11)] <- 0
			# for 2012 
				new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
				test2_12recs <- ddply(new1112,~SiteID,summarise,recruits=mean(z1))
				add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
				test2_12fruits <- ddply(add,~SiteID,summarise,fruits=sum(Fec))
				test2_12 <- merge(test2_12fruits, test2_12recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_12[is.na(test2_12)] <- 0
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test2_11$fruits+1), log(test2_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log MEAN all recruit size @ t", bg="green", main="SiteID 2011-green & 2012-blue", pch=21)
		points(log(test2_12$fruits+1), log(test2_12$recruits+1), bg="blue", pch=21)
	
	# MAYBE EVEN REGIONS???
					# for 2011 
					new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
					test3_11recs <- ddply(new1011,~Region,summarise,recruits=mean(z1))
					add <- getNew[which(getNew$TotFr_10 > 0), ]
					test3_11fruits <- ddply(add,~Region,summarise,fruits=sum(TotFr_10))
					test3_11 <- merge(test3_11fruits, test3_11recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_11[is.na(test3_11)] <- 0
				# for 2012 
					new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
					test3_12recs <- ddply(new1112,~Region,summarise,recruits=mean(z1))
					add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
					test3_12fruits <- ddply(add,~Region,summarise,fruits=sum(Fec))
					test3_12 <- merge(test3_12fruits, test3_12recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_12[is.na(test3_12)] <- 0
		# fruits last year in plot & recruits this year in plot 2011
			plot(log(test3_11$fruits+1), log(test3_11$recruits+1), xlab="Log fruits @ t-1", 
			ylab="log MEAN all recruit size @ t", bg="green", main="Region 2011-green & 2012-blue", pch=21)
			points(log(test3_12$fruits+1), log(test3_12$recruits+1), bg="blue", pch=21)
		
		
		
		
		
		
	####################################
	####################################
	####################################
	####################################
	# MEDIAN RECRUIT SIZE VS. NUMBER OF FRUITS DROPPED?

	table(dfclean$New)
	# summarize data
	library(plyr)
	
	# have to add on fruit from 2010 for 2011 seedlings 
		setwd(path.dat)
		getNew <- read.csv(file="Demog.data.2010-2012.forMatt.csv")
		add <- getNew[which(getNew$TotFr_10 > 0), ]
	
	# for 2011 
		new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
		test1_11recs <- ddply(new1011,~PlotID,summarise,recruits=median(z1))
		test1_11fruits <- ddply(add,~PlotID,summarise,fruits=sum(TotFr_10))
		test1_11 <- merge(test1_11fruits, test1_11recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_11[is.na(test1_11)] <- 0
	
	# for 2012 
		new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
		test1_12recs <- ddply(new1112,~PlotID,summarise,recruits=median(z1))
		add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
		test1_12fruits <- ddply(add,~PlotID,summarise,fruits=sum(Fec))
		test1_12 <- merge(test1_12fruits, test1_12recs, by.x = "PlotID", by.y = "PlotID", all.x = TRUE, all.y = TRUE)
		test1_12[is.na(test1_12)] <- 0
			
			
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test1_11$fruits+1), log(test1_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log MEDIAN all recruit size @ t", bg="green", main="PlotID 2011-green & 2012-blue", pch=21)
		points(log(test1_12$fruits+1), log(test1_12$recruits+1), bg="blue", pch=21)
	#********* pretty much no clear relationship at all!!!!
	# maybe fruits dropped at site in the previous year might be a better predictor
	# For SITES
			# for 2011 
				new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
				test2_11recs <- ddply(new1011,~SiteID,summarise,recruits=median(z1))
				add <- getNew[which(getNew$TotFr_10 > 0), ]
				test2_11fruits <- ddply(add,~SiteID,summarise,fruits=sum(TotFr_10))
				test2_11 <- merge(test2_11fruits, test2_11recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_11[is.na(test2_11)] <- 0
			# for 2012 
				new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
				test2_12recs <- ddply(new1112,~SiteID,summarise,recruits=median(z1))
				add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
				test2_12fruits <- ddply(add,~SiteID,summarise,fruits=sum(Fec))
				test2_12 <- merge(test2_12fruits, test2_12recs, by.x = "SiteID", by.y = "SiteID", all.x = TRUE, all.y = TRUE)
				test2_12[is.na(test2_12)] <- 0
	# fruits last year in plot & recruits this year in plot 2011
		plot(log(test2_11$fruits+1), log(test2_11$recruits+1), xlab="Log fruits @ t-1", 
		ylab="log MEDIAN all recruit size @ t", bg="green", main="SiteID 2011-green & 2012-blue", pch=21)
		points(log(test2_12$fruits+1), log(test2_12$recruits+1), bg="blue", pch=21)
	
	# MAYBE EVEN REGIONS???
					# for 2011 
					new1011 <- dfclean[which(dfclean$New==1 & dfclean$Year==1011), ]
					test3_11recs <- ddply(new1011,~Region,summarise,recruits=median(z1))
					add <- getNew[which(getNew$TotFr_10 > 0), ]
					test3_11fruits <- ddply(add,~Region,summarise,fruits=sum(TotFr_10))
					test3_11 <- merge(test3_11fruits, test3_11recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_11[is.na(test3_11)] <- 0
				# for 2012 
					new1112 <- dfclean[which(dfclean$New==1 & dfclean$Year==1112), ]
					test3_12recs <- ddply(new1112,~Region,summarise,recruits=median(z1))
					add <- dfclean[which(dfclean$Fec > 0 & dfclean$Year==1112), ]
					test3_12fruits <- ddply(add,~Region,summarise,fruits=sum(Fec))
					test3_12 <- merge(test3_12fruits, test3_12recs, by.x = "Region", by.y = "Region", all.x = TRUE, all.y = TRUE)
					test3_12[is.na(test3_12)] <- 0
		# fruits last year in plot & recruits this year in plot 2011
			plot(log(test3_11$fruits+1), log(test3_11$recruits+1), xlab="Log fruits @ t-1", 
			ylab="log MEDIAN all recruit size @ t", bg="green", main="Region 2011-green & 2012-blue", pch=21)
			points(log(test3_12$fruits+1), log(test3_12$recruits+1), bg="blue", pch=21)
			
			
			dev.off()
		