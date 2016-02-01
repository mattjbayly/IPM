# File: 00_SetDirectories.R 
######################################################
# SETUP DIRECTORIES & FILE PATHWAYS
######################################################
# Initiate folder directories & layout so later scripts  
# can quickly source, write & save files
######################################################
# path.root = main master folder for project with all other files contained(root) 
# path.dat = raw data files & main data files for project
# path.code = extra source code & functions developed by others for specific models ect. (not yet in specific r packages) 
# path.obj = R output & files produced from R
			
# SET ACCORDING TO YOUR SPECIFIC COMPUTER
path.root="C:/Users/DW/Desktop/transplant_analysis/" 
	
# SHOULD REMIAN UNCHANGED
path.dat=paste(path.root, "/Data", sep="") # MAIN DATA FILES
path.dat.raw=paste(path.dat, "/raw_data", sep="") # Don't write here!
path.code=paste(path.root, "/Rcode", sep="") # Source Code
path.obj=paste(path.root, "/Robjects", sep="") # To store temporary objects
path.fig=paste(path.root, "/Figures", sep="") # Polished & exploratory figures
setwd(path.dat); setwd(path.dat.raw); setwd(path.code); setwd(path.fig); setwd(path.obj)
path.funct = paste(path.root,"/Rcode/functions", sep="")

###############################################
# ALL OTHER SUBSEQUENT SCRIPTS WILL SOURCE THIS FILE
# via ... 
#path.root="C:/Users/DW/Desktop/temporalAUC/" 
#setwd(path.root)
#source("1_Setup_directoris.R")

setwd(path.dat) # set working directories
