#title: "Site grouping tests, performance within and 
		#beyond range limit and across occupied vs. unoccupied sites"
#author: "Matthew Bayly"
#modified: "Saturday, January 20, 2015"

# PURPOSE OF SCRIPT: 
# To compare the performance of transplants:
	# i.) Across the range limit (within and beyond)
	# ii.) Across occupied vs. unoccupied sites.

# Assuming no climatic differences across sites, did
# sites within the range limit out-preform those beyond 
# the range limit (as would be expected from fitness limitation 
# limitation)? 
# Assume occupied sites within the range limit would outperform
# unoccupied sites. 

# For each vital rate run a mixed effect model with siteType as a varible
# and sites nested within siteType 
# e.g. lmer( y ~ x + siteType + (1|siteType/site/plot) ) 
# siteType should be a fixed categorical factor. 

# INDEX: 
	# 0. load libraries & set directories to computer
	# 1. SITETYPE GROUP TESTS: 
		# 1.1 MEM with sitType as "within vs. beyond"
			# 1.1.1 SURVIVORSHIP
			# 1.1.2 GROWTH
			# 1.1.3 pFLOWER
			# 1.1.4 FECUNDITY
		# 1.2 MEM with sitType as "occupied vs. unoccupied"
		 #* exclude all sites beyond RL.
			# 1.2.1 SURVIVORSHIP
			# 1.2.2 GROWTH
			# 1.2.3 pFLOWER
			# 1.2.4 FECUNDITY
		# 1.3 MEM with site type as "occupied vs. unoccWithin vs. unoccBeyond"
			# 1.3.1 SURVIVORSHIP
			# 1.3.2 GROWTH
			# 1.3.3 pFLOWER
			# 1.3.4 FECUNDITY
	
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