#---
#title: "1.Data_structure"
#author: "Matthew Bayly"
#date: "August 26, 2015"
#---

#**EXPLORE DATA STRUCTURE**

# This script is for exploratory data analysis. The purpose 
# of this script it to generate a series of tables
# that look for anomalies in the data structure 
# & values prior to beginning the fulldata analysis.
# This data structure exploratory step is run in
# addition to precautions taken during data entry. 

## **INDEX**

# 1. SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS
# 2. Quick plots for size-specific transitions in the global dataset 
# 3. Revise dataframe change variable names for easier coding, make binary factors for survivorship
# 4. Custom individual growth tracking plots
# 5. Check data structure & do log transformation (where necessary).
# 6. Size/stage starting state variable (what should be used as start size from greenhouse?). 
# 7. Attach plot level variables 
# 8. Attach SITE-level variables 
# 9. Save new dataframe for other use with other scripts


#============================================================================================#
# 1. START: SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS
#============================================================================================#

### _Set directories for computer_ ###
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"	
setwd(path.set)
source("00_SetDirectories.R") # directory script (edit for your own computer). 
setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)

#Open 2015 plant datafile 
setwd(path.dat.raw)
plantdat <- read.csv(file="plant_data2015.csv")
dim(plantdat); #colnames(plantdat)

# Check whole data structure
#str(plantdat)

# **Check level **: are there any wierd typos or other anomalies
# are all values plausible? 
summary(plantdat$year) #OK
levels(plantdat$plot) # ''(empty)
levels(plantdat$planted) # OK
levels(plantdat$census.posttrans) # OK
levels(plantdat$census.fall) # OK
summary(plantdat$y.cord); summary(plantdat$x.cord) # OK
levels(plantdat$site) # NEED TO DELETE EMPTY ROWS AT BOTTOM OF DATAFRAME
levels(plantdat$pot) # Problems: 'X', 'MARKER', ''(empty), 'seed'.
levels(plantdat$source) # Problems: 'RR', ''(empty).
#plantdat$start.height
summary(plantdat$start.bsd) #Problems: "?" 
levels(plantdat$tag.marker)
levels(plantdat$post.trans.condition) #Problems: "?" 
levels(plantdat$mid.condition) #Problems:  ~ COMMENTS
levels(plantdat$fall.condition) #Problems: "D", "DRY", ~ COMMENTS
#summary(plantdat)
#FALL.FR.full max = 37!?!? too high!


#**Check quantities**: Are plants, plots, pots & sources represented?
# Table of # Plants at each site
table(plantdat$site) # seems ok (but, Hunter slightly too high!?)

# Table of # Plots at each site
table(plantdat$plot) # across whole study
  temp <- as.matrix(table(plantdat$plot)); temp <- temp/9; temp # per site

# Table of # of R & C at each site
with(plantdat, table(plantdat$source, plantdat$site))

# Table of # of pot sizes at each site & source
with(plantdat, table(plantdat$pot, plantdat$site))
with(plantdat, table(plantdat$pot, plantdat$source))

# Table of plants in plot per site
  # Most useful!
  with(plantdat, table(plantdat$plot, plantdat$site))
  # save table 
  setwd(path.obj); 
  write.csv(with(plantdat, table(plantdat$plot, plantdat$site)), file="plants_per_plot.csv")

### 3. Quick data structure plots

#***Quick data structure plot to visualize***:
# Plots of above for each site 
setwd(paste0(path.fig, "/01_DataStruct"))
pdf(file="01_DataStruct_global_2015.pdf", width=11, height=8.5)
plot(as.numeric(plantdat$start.height), log(plantdat$start.bsd + 0.5)) # correlation extremely weak as expected 
boxplot(plantdat$start.height~plantdat$pot, ylab="Start Height") # some very suspect values
boxplot(plantdat$start.bsd~plantdat$pot, ylab="start BSD") # some very suspect values (including zeros & NAs?)
# medium pots also tended to have more branching
boxplot(plantdat$start.height~plantdat$source,  ylab="Start Height") # seems reasonable
### for site & pot size

# Additional tables
with(plantdat, table(plantdat$source, plantdat$plot))
with(plantdat, table(plantdat$pot, plantdat$plot))
with(plantdat, table(plantdat$plot, plantdat$site))
dev.off()
# END: SET DIRECTORIES, LOAD LIBRARIES, LOAD PLANT DATA FILE, EXAMIN LEVELS


#============================================================================================#
# 2. START: Quick plots for size-specific transitions in the global dataset 
#============================================================================================#

#***SUMMER TO FALL GROWTH***
	MID.total.height <- plantdat$spring.total.height # might need to change this
	FALL.total.height <- plantdat$fall.total.height # might need to change this
	plot(MID.total.height, FALL.total.height, xlab="July Total Stem Length", ylab="FALL Total Stem Length", col=as.factor(plantdat$site), main="col of Sites global")
	abline(0, 1)

# BY SITE - split down by site for easier viewing (LOOP)
	site <- levels(plantdat$site); site <- as.factor(site)
	#op <- par(ask=FALSE)
	setwd(paste0(path.fig, "/01_DataStruct"))
	pdf(file="01_DataStruct_sites_2015.pdf", width=11, height=8.5)
	for (i in 1:length(site)){
	   current <- plantdat[plantdat$site==site[i], ]
	  plot(current$spring.total.height, current$fall.total.height, xlab="July tot stm", ylab="FALL tot stm", 
			  col=as.factor(current$pot), main=site[i], ylim=c(0, 390), xlim=c(0,250), pch=19, 
			  cex=0.7); abline(0, 1)
	  }
	#par(op)

# BY POT SIZE TYPE - split down by pots for easier viewing (LOOP)
pot <- levels(plantdat$pot); pot <- as.factor(pot)
	#op <- par(ask=TRUE)
	for (i in 1:length(pot)){
	   current <- plantdat[plantdat$pot==pot[i], ]
	  plot(current$spring.total.height, current$fall.total.height, xlab="July tot stm", ylab="FALL tot stm", 
			  col=as.factor(current$site), main=pot[i], ylim=c(0, 390), xlim=c(0,250), pch=19, 
			  cex=0.7); abline(0, 1)
	  }
	#par(op)
	#op <- par(ask=FALSE)
	dev.off(); dev.off()

# END: Quick plots for size-specific transitions in the global dataset 


#============================================================================================#
# 3. START: Revise dataframe change variable names for easier coding, make binary facotrs for survivorship
#============================================================================================#
# MAKE REVISED DATAFRAME TO WORK WITH
# subset down useful variables for current analysis
revised_data <- plantdat[ ,c("ID", "site", "plot", "pot", 
		"source", "start.height", "post.trans.condition", 
	 "fr.total", "fr.fl.total", "spring.total.height",
	 "spring.stem.number", "spring.max.stem", "fall.total.height", 
	 "fall.stem.number", "fall.max.stem")]

	summary(revised_data) # check it out
#make a unique plot variable
	revised_data$Uplot <- paste(revised_data$site, revised_data$plot, sep="")
# remove SEED rows (not included for this analysis
	revised_data <- revised_data[ !grepl("SEED", revised_data$pot) , ] 
	# to drop plants with heavy herbivory at Looking Glass Creek
	#revised_data <- revised_data[ !grepl("HERB", revised_data$mid.condition) , ]

###############################################################################
# RENAME VARIABLES 
	#size3 - early spring growth size in May 2015 after planting in May 2014
	#size4 - final size t2, for growth to the end of the summer OR for whole season growth 
		#***used in final analysis
	#surv3 - survivorship to May 2015 (over winter survivorship). 
	#surv4 - survivorship to the end of the summer (late July 2015). 
		#***used in final analysis
	# (old) surv_end15 - survivorship from May planting to Sept 2014

##########################################################
# Rename variables to be more appropriate and easier to work with.
	names(revised_data)[names(revised_data) == 'spring.total.height'] <- 'size3' #size3= May 2015 total stem length
	names(revised_data)[names(revised_data) == 'fall.total.height'] <- 'size4' #size4= late July 2015 end of summer. 
	revised_data$size3_ln <- log(revised_data$size3 + 1)
	revised_data$size4_ln <- log(revised_data$size4 + 1)
# change data & substitute values for survival 
##############################################
	revised_data$surv0 <- revised_data$post.trans.condition # Post transplant survivorship
	revised_data$surv0 <- sub("A", "1", revised_data$surv0) # Alive = 1
	revised_data$surv0 <- sub("D", "0", revised_data$surv0) # Dead = 0
	revised_data$surv0 <- as.numeric(revised_data$surv0)
#########################################
# Survival for whole summer period = surv_end15 (June 2014 - The final census July 2015)
	# total end of season survival
	revised_data$surv_end15 <- revised_data$size4 # assuming only AUGUST census
	revised_data$surv_end15[revised_data$surv_end15 > 0 ] <- 1 # Make 0/1
	# exclude post transplant mortality from estimate
	revised_data$surv_end15[ revised_data$surv0 == 0] <- NA 

#########################################
# Merge data frames 2014 + 2015
# Will need to merge dataframes to get surv0, surv1, surv2, start, size1, size2 from 2014 data. 
	setwd(path.dat); dir()
	Data2014 <- read.csv(file="Data_2014.csv")
	Data2014 <- Data2014[ ,c("ID.x", "start", "size1", "size1_ln", "size2", "size2_ln", 
	"surv0", "surv1", "surv2", "pFlower", "fec")]
	# make sure not to mix up 2014 & 2015 fecundity data 
	names(Data2014)[names(Data2014)=="fec"] <- "fec2014"
	names(Data2014)[names(Data2014)=="pFlower"] <- "pFlower2014"

	revised_data <- merge(revised_data, Data2014, by.x="ID", by.y="ID.x")
	revised_data <- revised_data[, !(colnames(revised_data) %in% c("surv0.x",
	"post.trans.condition", "fr.total", "spring.stem.number", 
	"spring.max.stem", "fall.stem.number", "fall.max.stem"))]
	names(revised_data)[names(revised_data)=="surv0.y"] <- "surv0"
#########################################

# To look only at survivorship seasonally (don't keep zeros in inter-seasonal dataset)
# Surv 3 = overwinter survivorship
	revised_data$surv3[ revised_data$surv2 == 0 & revised_data$size3 == 0 ] <- NA
	revised_data$surv3[ revised_data$surv2 == 1 & revised_data$size3 == 0 ] <- 0
	revised_data$surv3[ revised_data$surv2 == 1 & revised_data$size3 >= 1 ] <- 1
# Surv4 - looking only at the May to July data
	revised_data$surv4[ revised_data$surv3 == 0 & revised_data$size4 == 0 ] <- NA
	revised_data$surv4[ revised_data$surv3 == 1 & revised_data$size4 == 0 ] <- 0
	revised_data$surv4[ revised_data$surv3 == 1 & revised_data$size4 >= 1 ] <- 1
# Surv_end15 - whether or not the plant survived from Sept 2014 - Final Census (July 2015)
	# (Survival from t1 to t+1, annual time step)
	revised_data$surv_end15 <- revised_data$size4 # start with final census 
	revised_data$surv_end15[ revised_data$surv2 == 0] <- NA  # remove all plants that didn't survive to the fall of 2014
	revised_data$surv_end15[is.na(revised_data$surv2)] <- NA  # and all plants that died before the fall of 2014
	revised_data$surv_end15[ revised_data$surv2 == 0 & revised_data$size4 == 0 ] <- NA
	revised_data$surv_end15[ revised_data$surv2 == 1 & revised_data$size4 == 0 ] <- 0
	revised_data$surv_end15[ revised_data$surv2 == 1 & revised_data$size4 >= 1 ] <- 1	

# Check values - do they add up?
	with(revised_data, table(revised_data$surv0, revised_data$site))
	with(revised_data, table(revised_data$surv1, revised_data$site))
	with(revised_data, table(revised_data$surv2, revised_data$site))
	with(revised_data, table(revised_data$surv3, revised_data$site))
	with(revised_data, table(revised_data$surv4, revised_data$site))
	with(revised_data, table(revised_data$surv_end15, revised_data$site))

#################################################
# Make growth conditional on survival
	revised_data$size3[ revised_data$size3 == 0 ] <- NA
	revised_data$size3_ln[ revised_data$size3_ln == 0 ] <- NA
	revised_data$size4[ revised_data$size4 == 0 ] <- NA
	revised_data$size4_ln[ revised_data$size4_ln == 0 ] <- NA
#################################################
# Change pot sizes to numerical values (replace values)
############################################
############################################
############################################
############################################

# START: final plot to verify revisions.
	d <- revised_data
	site <- levels(d$site); site <- as.factor(site) # study sites
	setwd(paste0(path.fig, "/01_DataStruct"))

#### GROWTH ACROSS SITES
	pdf(file="Growth_2015.pdf", width=11, height=8.5)
	par(mfrow=c(2,4),mar=c(4,4,2,1))
	for (i in 1:length(site)){
	   current <- d[d$site==site[i], ]
	   current$moot <- current$size2_ln - current$size1_ln; current$moot <- (current$moot + 2.5)/1.5
	  plot(current$size2_ln, current$size3_ln, xlab="2014 fall size", ylab="May - Spring size 2015", 
			  col=as.factor(current$Uplot), main=site[i], pch=19, 
			  cex=current$moot); abline(0, 1)
	  }
	for (i in 1:length(site)){
	   current <- d[d$site==site[i], ]
	 current$moot <- current$size2_ln - current$size1_ln; current$moot <- (current$moot + 2.5)/1.5
	plot(current$size2_ln, current$size4_ln, xlab="2014 fall size", ylab="July size 2015", 
			  col=as.factor(current$Uplot), main=site[i], pch=19, 
			  cex=current$moot); abline(0, 1)
	  }	  
	for (i in 1:length(site)){
	   current <- d[d$site==site[i], ]
	   current$moot <- current$size2_ln - current$size1_ln; current$moot <- (current$moot + 2.5)/1.5
	 plot(current$size3_ln, current$size4_ln, xlab="May - Spring 2015", ylab="July - summer 2015", 
			  col=as.factor(current$Uplot), main=site[i], pch=19, xlim=c(0,6), ylim=c(0,6), 
			  cex=current$moot); abline(0, 1)
	  }
	for (i in 1:length(site)){
	   current <- d[d$site==site[i], ]
	   current$increment1 <- current$size2_ln - current$size1_ln # 2014 growth
	   	current$increment2 <- current$size4_ln - current$size3_ln # 2015 growth
	 plot(current$increment1, current$increment2, xlab="2014 Growth Increment (ABS)", ylab="2015 Growth Increment (ABS)", 
			  col=as.factor(current$Uplot), main=site[i], pch=19, xlim=c(0,3.5), ylim=c(0,3.5), 
			  cex=1.5); abline(0, 1)
	  }  
	for (i in 1:length(site)){
	   current <- d[d$site==site[i], ]
	   current$increment1 <- current$size2_ln/current$size1_ln # 2014 growth
	   	current$increment2 <- current$size4_ln/current$size3_ln # 2015 growth
	 plot(current$increment1, current$increment2, xlab="2014 Growth Increment (REL)", ylab="2015 Growth Increment (REL)", 
			  col=as.factor(current$Uplot), main=site[i], pch=19, xlim=c(0.2,3.2), ylim=c(0.6,3.2), 
			  cex=1.5); abline(0, 1)
	  }  
	   
	  
dev.off(); dev.off()

  
#============================================================================================#
# 4. START: Custom individual growth tracking plots
#============================================================================================#

# STAR INDIVIDUAL TRACKING PLOT
#library(reshape)
#id <- seq(1:nrow(d))# unique plant id
#d <- cbind(d, id)
#d$begin <-  runif(nrow(d), 0, 4)
# order sites based on location rather than alphabetical
#d$site <- factor(d$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
#site <- levels(d$site); site <- as.factor(site) # slit sites
# melt to planting & monitoring time intervals as line
#mdata <- melt(d, id=c("ID", "site", "plot", "pot", "source", "start.height", 
# "post.trans.condition", "mid.condition", 
# "fall.condition", "fr.total", "fr.fl.total", "size3", 
 #"MID.stem.number", "MID.max.stem", "size4", "FALL.stem.number", 
# "FALL.max.stem", "Uplot", "size3_ln", "size4_ln", "surv0", "surv_end15",
# "surv3", "surv4"))
#dead_planting <- subset(mdata, mdata$post.trans.condition=="D")  # plants dead ~ 1 week after transplanting
#dead_planting$value <- dead_planting$value - 10 # just to put it below the line
# set ploting frame & id column for quick plant reference in plotting 
#par(mfrow=c(2,4),mar=c(4,4,2,1)) 
#id <- as.factor(mdata$id)

# SUPER PLOT FOR LARGE PLANTS (INDIVIDUALS)
#for (i in 1:length(site)){
#   current <- mdata[mdata$site==site[i], ]
#	plot(current$size, current$value, type="n", xlab="Greenhouse", ylab="GROWTH", main=site[i], ylim=c(-15, 390), xlim=c(0, 390)) 
#		points(dead_planting$size, dead_planting$value, pch=21, bg="lightgrey", cex = 1.4) # dead transplants
#		#each plant has its own line
#		for (j in 1:length(id)) { 
#			plant <- current[current$id==id[j], ]
#			plant2 <- subset(plant, plant$variable=="MID.total.height" | plant$variable=="begin")  # RED just for spring growth
#				lines(plant2$size, plant2$value, col="red", type="b", lwd=3) 
#				lines(plant$size, plant$value, type="b", lwd=1.5, 
#						bg=as.factor(plant$variable), pch=21, cex = 1.4) 				
#
# }
#  }

# END: Custom individual growth tracking plots
 


#============================================================================================#
# 5. START: Check data structure & do log transformation (where necessary).
#============================================================================================#

d <- revised_data
d$start.height_ln <- log(d$start.height + 1)
colnames(d) 
site <- levels(d$site); site <- as.factor(site) # study sites
d$fec2015 <- d$fr.fl.total
d$fec2015[d$fec2015==0]<- NA
# binary response variable for flowering probability, given survival to fall
d$pFlower2015[ d$surv4 == 0 ] <- NA
d$pFlower2015[ d$surv4 == 1 & d$fr.fl.total == 0 ] <- 0
d$pFlower2015[ d$surv4 == 1 & d$fr.fl.total >= 1 ] <- 1
#summary & str
summary(d)
str(d)

# END: Check data structure & do log transformation (where necessary).


#============================================================================================#
# 6. START: Size/stage starting state variable (what should be used as start size from greenhouse?). 
#============================================================================================#

# From data exploration script.
d$start <- 0.01*(d$start.height_ln) + 0.05*(d$pot) +  0.1*(d$start.height_ln*d$pot)# multiplicative

#  END: Size/stage starting state variable (what should be used as start size from greenhouse?). 


#============================================================================================#
# 7. START: Attach plot level variables 
#============================================================================================#
setwd(path.dat.raw) # where environmental data is stored
enviro <- read.csv(file="environmental_variables.csv") # plot level enviornmental variables
colnames(enviro) # seems ok
levels(enviro$SITE)[levels(enviro$SITE)=="HUNT"] <- "HUNTER"
levels(d$site); levels(enviro$SITE) 
levels(d$plot); levels(enviro$PLOT) 
enviro$Un.plot <- paste(enviro$SITE, enviro$PLOT, sep="") # unique plot variable to mere with 
###### will need to merge these data frames
head(enviro$Un.plot); head(d$Uplot)

##### MERGE
d <- merge(d, enviro, by.x = "Uplot", by.y = "Un.plot", all.x = TRUE, all.y = FALSE)

# END: Attach plot-level variables 


#============================================================================================#
# 8. START: Attach SITE-level variables 
#============================================================================================#

setwd(path.dat.raw) # where site data is stored
SDMsite <- read.csv(file="occ_site_preds_sept2014.csv") # site level enviornmental variables
colnames(SDMsite) # seems ok
levels(d$site); levels(SDMsite$MERGE) # HUNTER needs to be renamed to match
levels(d$plot); levels(enviro$PLOT) # HUNTER needs to be renamed to match

##### MERGE
d <- merge(d, SDMsite, by.x = "site", by.y = "MERGE", all.x = TRUE, all.y = FALSE)
d$site_type <- d$site
	library(plyr)
	revalue(d$site_type, c("CALAPOOIA" = "beyond")) -> d$site_type
	revalue(d$site_type, c("WILEY" = "beyond")) -> d$site_type
	revalue(d$site_type, c("HUNTER" = "beyond")) -> d$site_type
	revalue(d$site_type, c("THOMAS" = "beyond")) -> d$site_type
	revalue(d$site_type, c("ROCK" = "occupied")) -> d$site_type
	revalue(d$site_type, c("COAST" = "occupied")) -> d$site_type
	revalue(d$site_type, c("LOOK" = "unoccupied")) -> d$site_type
	revalue(d$site_type, c("MOSBY" = "unoccupied")) -> d$site_type
levels(d$site_type)

# END: Attach SITE-level variables

	
#============================================================================================#
# 9. Save new dataframe for other use with other scripts
#============================================================================================#

###########################################
# Standardize Environmental Variables
d$start_sd <- ((d$start - mean(d$start, na.rm=TRUE))/sd(d$start, na.rm=TRUE))
d$ENSEMBLE_sd <- ((d$ENSEMBLE - mean(d$ENSEMBLE, na.rm=TRUE))/sd(d$ENSEMBLE, na.rm=TRUE))
d$moist_score_sd <- ((d$moist_score - mean(d$moist_score, na.rm=TRUE))/sd(d$moist_score, na.rm=TRUE))
###########################################


# trim down DF
names(d)[names(d)=="ID.x"] <- "ID"
d <- d[,c("ID", "site", "Uplot","surv0", "surv1", "surv2", "surv3", "surv4", "surv_end15", "start", 
	"size1", "size2", "size3", "size4", "size1_ln", "size2_ln", "size3_ln",
	"size4_ln", "pFlower2014", "fec2014", "pFlower2015", "fec2015")]
#########################################

# Save & check out 
setwd(path.dat)
dir()
write.csv(d, file="Data_2015.csv", row.names = FALSE)
# END
