#title: "MAJOR SOURCE OF VARIATION SPATIAL & TEMPORAL"
#author: "Matthew Bayly"
#modified: "Saturday, January 17, 2015"

# PURPOSE OF SCRIPT: 
# Quick overview of major source of variability for each 
# vital rate in terms of spatial (within sites vs. between
# sites) and temporal (May - July Vs. July - August vs. 
# May - August?) coverage of data. Prior to diving into 
# exploration of all combinations of environmental variables, 
# a quick overview of variability is useful to justify the 
# data analysis plan. 

# Goal: to produce final figures similar to Figure 2 in:
	#Ramula (2014) Linking vital rates to invasiveness of a perennial herb
	#Oecologia, 174:1255–1264.
		
# INDEX: 
	# 0. Load libraries & set directories to computer

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


######################################
# TEMPORARY STANDARD DEVIATION VALUES 


#STANDARD DEVIATION:
#ratios plot to site

#Survival
#	May - July
plot <- 1.803/(1.803+1.474)
site <- 1.474/(1.803+1.474)
level <- "S.May_July"
S.May_July <- data.frame(level, plot, site)
#		Within sites 1.803
#		Across sites 1.474
#	July - Aug
plot <- 2.438/(2.438+2.230)
site <- 2.230/(2.438+2.230)
level <- "S.July_Sept"
S.July_Sept <- data.frame(level, plot, site)
#		Within sites 2.438
#		Across sites 2.230
#	May - Aug
plot <- 2.539/(2.539+2.126)
site <- 2.126/(2.539+2.126)
level <- "S.May_Sept"
S.May_Sept <- data.frame(level, plot, site)
#		Within sites 2.539 
#		Across sites 2.126

#Growth 
#	May - July
plot <- 0.5685559/(0.5685559+0.3750124)
site <- 0.3750124/(0.5685559+0.3750124)
level <- "G.May_July"
G.May_July <- data.frame(level, plot, site)
#		Within sites 0.5685559
#		Across sites 0.3750124
#		resid 0.6435375
#	July - Aug
plot <- 0.3743681/(0.3743681+0.02571898)
site <- 0.02571898/(0.3743681+0.02571898)
level <- "G.July_Sept"
G.July_Sept <- data.frame(level, plot, site)
#		Within sites 0.3743681 
#		Across sites 0.02571898
#		resid 0.4520948
#	May - Aug
plot <- 0.836431/(0.836431+0.2965563)
site <- 0.2965563/(0.836431+0.2965563)
level <- "G.May_Sept"
G.May_Sept <- data.frame(level, plot, site)
#		Within sites 0.836431 
#		Across sites 0.2965563	
#		resid 0.6527844

#pFlower
#	Aug, given size
plot <- 1.330/(1.330+1.321)
site <- 1.321/(1.330+1.321)
level <- "pFlower"
pFlower <- data.frame(level, plot, site)
#		Within sites 1.330
#		Across sites 1.321	 

#Fecundity
#	Fruit Aug, given size 
plot <- 0.4994/(0.4994+0.3991)
site <- 0.3991/(0.499+0.3991)
level <- "Fec"
Fec <- data.frame(level, plot, site)

Surv <- S.May_Sept; Growth <- G.May_Sept; Fecundity <- Fec
df <- rbind(Surv, Growth, pFlower, Fecundity)

library(reshape)
mydata <- melt(df, id=c("level"))

library(ggplot2)
p <- ggplot(mydata, aes(x = level, y = value, fill=variable)) +
    geom_bar(stat='identity')

#eliminates background, gridlines, and chart border
p + theme_classic() + xlab("") + ylab("Source of additional variability")

	