###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   Create pretty map from model predictions
##   Project to region
## last update:  March 25, 2015 
###################################################################################

# load libraries now if desired; otherwise loaded below as needed
library(raster)
library(maptools)
library(gam)
library(randomForest)
library(dismo)
library(gbm)
library(rgdal)


#http://gis.stackexchange.com/questions/63577/joining-polygons
#http://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon

## set pathnames - Amy
path.root="/Users/amyangert/Desktop/Amy work August 2014" 
path.dat = paste(path.root, "/data files", sep="")
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.dat, "/ecoregions.shp", sep="")
path.bio = paste(path.dat, "/wc0.5", sep="")
path.sta = paste(path.dat, "/gz_2010_us_040_00_500k", sep="")
path.fig = paste(path.root, "/figures", sep="")

## set pathnames - Matthew
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.obj, "/ecoregions.shp", sep="")
path.bio = paste(path.obj, "/wc0.5", sep="")
path.cod=paste(path.root, "/R code", sep="")
#path.fig=paste(path.root, "/figures", sep="")
path.sta = paste(path.obj, "/gz_2010_us_040_00_500k", sep="")


################################################################################
######## START INITIALIZATION OF DATA STRUCTURES

## read in replicate training datasets
library(sp)
library(rgdal)
setwd(path.dat) 
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

## read in occupancy dataset
library(sp)
library(rgdal)
setwd(path.dat)
all = read.csv("all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
occ = all[all$DATASET=="occ",] #pull out occupancy dataset
dim(occ) #check size
occ$bio3 = log(occ$bio3+0.5) #make needed ln-transforms of predictors
occ$bio10 = log(occ$bio10+0.5)
occ$bio12 = log(occ$bio12+0.5)
occ$bio14 = log(occ$bio14+0.5)
coordinates(occ) = ~Longitude + Latitude
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
proj4string(occ) = CRS(prj.wgs)
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
occ.lcc = spTransform(occ, CRS=CRS(prj.lcc))

## split by elevation 
L.thresh = 400
H.thresh = 1200
high <- occ[occ$Elevation > H.thresh,]
low <- occ[occ$Elevation < L.thresh,]
mid <- occ[occ$Elevation >= L.thresh & occ$Elevation <= H.thresh,]
high.aea = occ.aea[occ.aea$Elevation > H.thresh,]
low.aea <- occ.aea[occ.aea$Elevation < L.thresh,]
mid.aea <- occ.aea[occ.aea$Elevation >= L.thresh & occ.aea$Elevation <= H.thresh,]

high.lcc = occ.lcc[occ.lcc$Elevation > H.thresh,]
low.lcc <- occ.lcc[occ.lcc$Elevation < L.thresh,]
mid.lcc <- occ.lcc[occ.lcc$Elevation >= L.thresh & occ.lcc$Elevation <= H.thresh,]


## read in occupancy dataset split by ecoregions
setwd(path.dat)
north = read.csv("north.csv")
center = read.csv("center.csv")
south = read.csv("south.csv")
prj.wgs = "+proj=longlat +ellps=WGS84"

coordinates(north) = ~Longitude + Latitude
proj4string(north) = CRS(prj.wgs)
coordinates(center) = ~Longitude + Latitude
proj4string(center) = CRS(prj.wgs)
coordinates(south) = ~Longitude + Latitude
proj4string(south) = CRS(prj.wgs)

## reproject to albers equal area
library(rgdal)
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
north.aea = spTransform(north, CRS=CRS(prj.aea))
center.aea = spTransform(center, CRS=CRS(prj.aea))
south.aea = spTransform(south, CRS=CRS(prj.aea))
north.lcc = spTransform(north, CRS=CRS(prj.lcc))
center.lcc = spTransform(center, CRS=CRS(prj.lcc))
south.lcc = spTransform(south, CRS=CRS(prj.lcc))

library(raster) 

###################
###################
###################
###################
###################
# 
# crop to buffered bounding box

occ2 = spTransform(occ, CRS=CRS(prj.lcc)) #ADDED
ymin = min(occ$Latitude) - 1.5
ymax = max(occ$Latitude) + 1.5
xmin = min(occ$Longitude) - 1.5
xmax = max(occ$Longitude) + 1.5
e = extent(xmin, xmax, ymin, ymax)

###################
###################
###################

setwd(path.eco)
# see 'help' for this function, you  need to set the directory to one folder level above where the destination file is saved 
ecoreg = readOGR(dsn="C:/Users/DW/Desktop/temp.sept.30/R objects/ecoregions.shp", layer="us_eco_l3_no_st") #ADDED

# project to match bioclim
library(rgdal)
prj.wgs = "+proj=longlat +ellps=WGS84"
ecoreg2 = spTransform(ecoreg, CRS(prj.lcc))
ecoreg.wgs = spTransform(ecoreg, CRS(prj.wgs))
ecoreg.aea = spTransform(ecoreg, CRS(prj.aea))


# crop polygon to buffered bounding box
# note: this loses the region names, so it doesn't work for subsetting data
library(rgeos)
bbox = as(e, "SpatialPolygons")
proj4string(bbox) = CRS(prj.wgs)
bbox.lcc = spTransform(bbox, CRS=CRS(prj.lcc))
bbox.aea = spTransform(bbox, CRS=CRS(prj.aea))
#ecoreg.crop.lcc = gIntersection(ecoreg2, bbox, byid=TRUE)
ecoreg.crop.lcc <- crop(ecoreg2, bbox.lcc) # ADDED

# pull out relevant ecoregions
eco.names = levels(ecoreg2$US_L3NAME)
poly.list = c(10,12,13,14,19,25,36,43,64,66,70,71,84)
for (i in 1:length(poly.list)) {
	poly = ecoreg2[ecoreg2$US_L3NAME==eco.names[poly.list[i]],]
	assign(paste("poly",poly.list[i], sep=""), poly)
	poly.aea = spTransform(poly, CRS=CRS(prj.aea))
	poly.lcc = spTransform(poly, CRS=CRS(prj.lcc))
	assign(paste("poly",poly.list[i],".lcc", sep=""), poly.lcc)
	assign(paste("poly",poly.list[i],".aea", sep=""), poly.aea)
	}

## world map polygons
library(rworldmap); data(countriesLow)
library(rworldxtra); data(countriesHigh)
library(rgdal)
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
countriesLow.aea = spTransform(countriesLow, CRS=CRS(prj.aea))
countriesHigh.aea = spTransform(countriesHigh, CRS=CRS(prj.aea))
countriesLow.lcc = spTransform(countriesLow, CRS=CRS(prj.lcc))
countriesHigh.lcc = spTransform(countriesHigh, CRS=CRS(prj.lcc))

## state polygons
setwd(path.sta)
sta = readShapePoly("gz_2010_us_040_00_500k.shp")
projection(sta) = CRS(prj.wgs)
sta.aea = spTransform(sta, CRS=CRS(prj.aea))
sta.lcc = spTransform(sta, CRS=CRS(prj.lcc))

## gridlines
library(sp)
library(rgdal)#
#create unprojected gridlines
ecoreg.wgs = spTransform(ecoreg, CRS=CRS(prj.wgs))
grd.wgs = gridlines(ecoreg.wgs, ndiscr=200)
grd.wgs2 = gridlines(ecoreg.wgs, ndiscr=2000)


#project gridlines
grd.aea = spTransform(grd.wgs, CRS=CRS(prj.aea))
grd.lcc = spTransform(grd.wgs, CRS=CRS(prj.lcc))

#prepare labels for gridlines in unprojected space
gridat <- gridat(ecoreg.wgs, side="EN")
#project labels for gridlines
gridat.aea=spTransform(gridat, CRS=CRS(prj.aea), side="EN")
gridat.lcc=spTransform(gridat, CRS=CRS(prj.lcc), side="EN")

#slim down labels to fit on plots
lab.all = parse(text=as.character(gridat.aea$labels))
lab.cull = lab.all[c(2,4,6:9)]
coord.all = coordinates(gridat.aea)
coord.cull = coord.all[c(2,4,6:9),]

lab.all = parse(text=as.character(gridat.lcc$labels))
lab.cull = lab.all[c(2,4,6:9)]
coord.all = coordinates(gridat.lcc)
coord.cull = coord.all[c(2,4,6:9),]

######## END INITIALIZATION OF DATA STRUCTURES
################################################################################
  
#Clip function to trim probability rasters by bbox.lcc (extent obj)
clip<-function(raster,shape) {
          a1_crop<-crop(raster,shape)
          step1<-rasterize(shape,a1_crop)
          a1_crop*step1}

################################################################################
################################################################################
################################################################################
################################################################################
######## START PRETTY MAPS

#Make base plotting frame with gridlines
frame <- crop(countriesHigh.lcc, bbox.lcc); 
frame.aea <- crop(countriesHigh.aea, bbox.aea); 
frame.wgs = spTransform(frame, CRS=CRS(prj.wgs))
frame.grd = gridlines(frame.wgs, ndiscr=100)
frame.grd.lcc = spTransform(frame.grd, CRS=CRS(prj.lcc))
frame.grd.aea = spTransform(frame.grd, CRS=CRS(prj.aea))
plot(frame); plot(frame.grd.lcc, add=T) # check it out..OK?
plot(frame.aea); plot(frame.grd.aea, add=T) # check it out..OK?

## replicate herb + pseudo training sets
setwd(path.dat)
for (i in 1:10) {
	dat = read.csv(paste("dat",i,"c.csv", sep=""))
	dat = dat[,c(4:8,59:61,67:69,71:72)] #column indices changed 9/4/14
	assign(paste("dat",i, sep=""), dat)
	coordinates(dat) = ~Longitude + Latitude
	projection(dat) = prj.wgs
	dat = spTransform(dat, CRS=CRS(prj.aea))
	assign(paste("dat",i,".aea", sep=""), dat)
	}
rm(dat)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

path.fig <- "F:/THESIS_mobile/Figures & table plan/code and objects"
setwd(path.fig)
dir()
pdf(file = "2.1.MAP_herb.replicates_and_pseudos.pdf", width=(8.5 - (1.25 + 0.75)), height=(3), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(1,3))

	pres = all[all$PRESABS==1,]
	pres = pres[pres$DATASET=='herb',]
	coordinates(pres) = ~Longitude + Latitude
	projection(pres) = prj.wgs
	pres = spTransform(pres, CRS=CRS(prj.aea))
	plot(frame.aea, col="lightgrey", main="Herbarium Presence Records")
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(pres, pch=21, col="lightgrey", bg="black")
	backup <- pres
	#legend("bottomleft", legend=c("pres","abs"), pch=21, col="lightgrey", pt.bg=c("black","white"), bg="white", box.col="white")
	text(coord.cull, labels=lab.cull, pos=c(1,1,1,1,1,1,4,4,4), offset=0, col="black") #this works for single panels, but not multi-panel
	
for (i in 1:1) {
	dat = get(paste("dat",i,".aea", sep=""))
	pres = dat[dat$PRESABS==1,]
	abs = dat[dat$PRESABS==0,]
	plot(frame.aea, col="lightgrey", main="Sample Thinned Replicate")
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(abs, pch=21, col="lightgrey", bg="white")
	points(pres, pch=21, col="lightgrey", bg="black")
	legend("bottomleft", legend=c("pres","abs"), pch=21, col="lightgrey", pt.bg=c("black","white"), bg="white", box.col="white")
	text(coord.cull, labels=lab.cull, pos=c(1,1,1,1,1,1,4,4,4), offset=0, col="black") #this works for single panels, but not multi-panel
	}

	pres = all[all$PRESABS==1,]
	pres = pres[pres$DATASET=='occ',]
	coordinates(pres) = ~Longitude + Latitude
	projection(pres) = prj.wgs
	abs = all[all$PRESABS==0,]
	abs = abs[abs$DATASET=='occ',]
	coordinates(abs) = ~Longitude + Latitude
	projection(abs) = prj.wgs
	pres = spTransform(pres, CRS=CRS(prj.aea))
	abs = spTransform(abs, CRS=CRS(prj.aea))
	plot(frame.aea, col="lightgrey", main="Independent Testing Data")
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(abs, pch=21, col="lightgrey", bg="white")
	points(pres, pch=21, col="lightgrey", bg="black")
	legend("bottomleft", legend=c("pres","abs"), pch=21, col="lightgrey", pt.bg=c("black","white"), bg="white", box.col="white")
	text(coord.cull, labels=lab.cull, pos=c(1,1,1,1,1,1,4,4,4), offset=0, col="black") #this works for single panels, but not multi-panel
		
dev.off()

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



par(mfrow=c(1,1))
	setwd(path.dat)
	preds = read.csv("site_preds.csv")
	preds = preds[preds$ID1=='trans',]
	preds <- preds[order(-preds$Latitude),] 
	coordinates(preds) = ~Longitude + Latitude
	projection(preds) = prj.wgs
							  # -2077015 
	test = as(extent(-2267634, -2040000, 2530000, 2770000), "SpatialPolygons")					# convert extent box to shapefile (rectangle)
	proj4string(test) = prj.aea			# assign spatial projection to extent object
	test.aea = spTransform(test, CRS=CRS(prj.aea))

# Crop down DEM
	library(raster)
	dem <- raster("C:/Users/DW/Desktop/temp.sept.30/path.gis/PRISM_us_dem_800m_bil/PRISM_us_dem_800m_bil.bil")
	frame.wgs <- spTransform(frame.aea, CRS=CRS(prj.wgs))
	dem <- crop(dem, frame.wgs)
	test.wgs = spTransform(test.aea, CRS=CRS(prj.wgs))
	dem_zoom <- crop(dem, test.wgs)
	dem_zoom.aea <- projectRaster(dem_zoom, crs=CRS(prj.aea))
	dem_zoom.aea <- crop(dem_zoom.aea, test.aea)
	dem_zoom.aea[dem_zoom.aea < 1] <- NA
	
	dem <- aggregate(dem, fact=3, fun=mean)
	dem[dem < 1] <- NA
	dem.aea <- projectRaster(dem, crs=CRS(prj.aea))

	#plot(dem_zoom.aea, col=palette(gray(seq(0.6,0.98,len=150))))
	
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	# TRANSPLANT SITE MAP 
	setwd(path.fig)
	pdf(file = "2.2.TransPlant_SiteMap.pdf", width=(8.5 - (1.25 + 0.75)), height=(4.5), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/

	par(mfrow=c(1,2))
	# TRIM DOWN AND REPROJECT TO REGION 
	test.wgs = spTransform(test.aea, CRS=CRS(prj.wgs))
	plot(test.aea, border="purple")
	ecoreg.aea2 <- crop(ecoreg.aea, test.aea) 
	plot(dem_zoom.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE)
	plot(ecoreg.aea2, add=TRUE, lwd=1.7)
	preds = spTransform(preds, CRS=CRS(prj.aea))
	plot(preds, add=TRUE, pch=24, col="lightgrey", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	backup2 <- crop(backup, test.aea)
	points(backup2, pch=21, col="lightgrey", bg="black", cex=0.8)
	verts = gridlines(ecoreg.wgs, easts=c(seq(-124.75, -121.75, by=0.5)), norths=c(seq(42.75, 45.75, by=0.5)))
	gridlittle.aea = spTransform(verts, CRS=CRS(prj.aea))
	gridlittle.aea <- crop(gridlittle.aea, test.aea)
	plot(gridlittle.aea, add=TRUE, lty=2, col="darkgrey")
	points(preds, cex=1.5, pch=24, col="black", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	u <- data.frame(name=c("one", "two"), long=c(-2229695, -2229695), lat=c(2570000, 2690000))
	coordinates(u)= ~long+lat
	projection(u) = prj.aea
	myline <- Line(coordinates(u))
	lines(myline, lwd=4)
	text(-2252695, 2630000, "120 km", srt=90)
	plot(test.aea, border="purple", add=TRUE, lwd=1.7)
	# distance = 119km 
	par(xpd=TRUE)
	legend("bottom", inset=-0.072, legend=c("within","beyond", "wild population"), pch=c(24, 24, 21), pt.bg=c("red","blue", "black"), bg="white", border="black")

	pres = all[all$PRESABS==1,]
	pres = pres[pres$DATASET=='herb',]
	coordinates(pres) = ~Longitude + Latitude
	projection(pres) = prj.wgs
	pres = spTransform(pres, CRS=CRS(prj.aea))
	plot(frame.aea)
	plot(dem.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE,)
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(pres, pch=21, col="lightgrey", bg="black", cex=0.8)
	backup <- pres
	#legend("bottomleft", legend=c("pres","abs"), pch=21, col="lightgrey", pt.bg=c("black","white"), bg="white", box.col="white")
	text(coord.cull, labels=lab.cull, pos=c(1,1,1,1,1,1,4,4,4), offset=0, col="black") #this works for single panels, but not multi-panel
	plot(preds, add=TRUE, pch=24, col="black", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	plot(test.aea, add=TRUE, border="purple", lwd=1.5)
	dev.off()
	
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



modMap <- c("C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/final_rasters/ALL_LR.tif",
	"C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/final_rasters/ALL_GAM.tif",
	"C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/final_rasters/ALL_RF.tif",
	"C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/final_rasters/ALL_BRT.tif",
	"C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/final_rasters/ALL_MAX.tif")

cuty <- c("LR.mod2.extaccs.pseudo11.rda",
"GAM.mod4.extaccs.pseudo11.rda",
"RF.mod1.extaccs.pseudo11.rda",
"BRT.mod4.extaccs.pseudo11.rda",	
"MAX.mod1.extaccs.pseudo11.rda")

moddy <- c("GLM", "GAM", "RF", "BRT", "MAX")
	
	setwd(path.fig)
	pdf(file = "2.3.Map_Probs_ENMS.pdf", width=(8.5 - (1.25 + 0.75)), height=(9 - (0.75 + 0.75)), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(3,2))

for(i in 1:5){
	
# RUN the next 22 lines for the first time only 
# TRIM DOWN AND REPROJECT TO REGION 
#	test.wgs = spTransform(test.aea, CRS=CRS(prj.wgs))
#	setwd(path.obj)
#	tt <- raster(modMap[i])
#	big <- tt
#	big <- aggregate(big,fact=10,fun=max)
#	#plot(big)
#	big <- crop(big, frame.wgs)
#	big.aea <- projectRaster(big, crs=CRS(prj.aea))
#	big.aea[big.aea < 0.0001] <- NA
#	ttc <- crop(tt, test.wgs)
#	ttc[ttc < 0.0001] <- NA
#	tt2 <- aggregate(ttc,fact=4,fun=mean)
#	pr1 <- projectRaster(tt2, crs=CRS(prj.aea))
#	pr1 <- pr1
#	
#	# Save models to avoid repeating process
#	setwd("C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/TransplantFig")
#	writeRaster(pr1, filename=paste0(moddy[i], '.zoom.tif'))
#	writeRaster(big.aea, filename=paste0(moddy[i], '.full.tif'))
#	setwd(path.fig)

	# load in saved rasters 
	setwd("C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/TransplantFig")
	big.aea <- raster(paste0(moddy[i], ".full.tif"))
	pr1 <- raster(paste0(moddy[i], ".zoom.tif"))
	setwd(path.fig)


# GET MODEL CUTS 
	setwd(path.obj)
	cuts <- get(load(cuty[i]))
	cuts <- cuts[cuts$thresh=="SensSpec",]
	cuts <- mean(cuts$threshold)

# COLOR BREAK POINTS 
	breakUnd <- c(seq(0, cuts*1000))
	breakAbove <- c(seq(cuts*1000, 1000))
	colorUnd =colorRampPalette(c("blue", "green", "yellow"))
	colorAbove =colorRampPalette(c("orange", "red"))
	plotcolUnd <- colorUnd(length(breakUnd))
	plotcolAbove <- colorAbove(length(breakAbove))
	plotcolAbove <- add.alpha(plotcolAbove, alpha=0.6)
	plotcolUnd <- add.alpha(plotcolUnd, alpha=0.6)

##################################
# MAKE FINAL PLOT
	plot(test.aea, border="black", main=paste0(moddy[i], " transplant sites"))
	# add on DEM under
	#plot(dem_zoom.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE)
	
# PLOT ON RASTER 
	Und <- pr1
	Und[Und > (cuts*1000)] <- NA
	Und[Und < 1] <- NA
	Und <- crop(Und, test.aea)
	plot(Und, add=TRUE, breaks=breakUnd, legend=FALSE, col=plotcolUnd)

	above <- pr1
	above[above < (cuts*1000)] <- NA
	above <- crop(above, test.aea)
	plot(above, add=TRUE, breaks=breakAbove, legend=FALSE, col=plotcolAbove)
	
	# other stuff
	preds = spTransform(preds, CRS=CRS(prj.aea))
	plot(preds, add=TRUE, pch=24, col="lightgrey", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	backup2 <- crop(backup, test.aea)
	points(backup2, pch=21, col="lightgrey", bg="black", cex=0.8)

	# gridlines 
	verts = gridlines(ecoreg.wgs, easts=c(seq(-124.75, -121.75, by=0.5)), norths=c(seq(42.75, 45.75, by=0.5)))
	gridlittle.aea = spTransform(verts, CRS=CRS(prj.aea))
	gridlittle.aea <- crop(gridlittle.aea, test.aea)
	plot(gridlittle.aea, add=TRUE, lty=2, col="darkgrey")

	# ecoregions
	ecoreg.aea2 <- crop(ecoreg.aea, test.aea) 
	plot(ecoreg.aea2, add=TRUE, lwd=1.7)
	points(preds, cex=1.5, pch=24, col="black", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))

	# scale bar known distance = 119km 
	u <- data.frame(name=c("one", "two"), long=c(-2229695, -2229695), lat=c(2570000, 2690000))
	coordinates(u)= ~long+lat
	projection(u) = prj.aea
	myline <- Line(coordinates(u))
	lines(myline, lwd=2)
	text(-2252695, 2630000, "120 km", srt=90)
	# distance = 119km 

	# add on legend
	plot(pr1, legend.only=TRUE, legend.width=2, legend.shrink=0.75, breaks=c(seq(0,1000)), col=c(plotcolUnd, plotcolAbove), 
	axis.args=list(at=seq(0, 1000, by=100), labels=seq(0, 1, by=0.1), cex.axis=0.8),
	 legend.args=list(text='Predicted Suitability', side=4, font=2, line=2.5, cex=0.8))
	# labels to match original model output 0 -1 (rather than 0 - 1000)
	# add the purple frame for reference
	plot(test.aea, add=TRUE, border="purple", lwd=1.5)
		
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# MAKE BIG PLOT 
	Und <- big.aea
	Und[Und > (cuts*1000)] <- NA
	Und[Und < 1] <- NA
	above <- big.aea
	above[above < (cuts*1000)] <- NA
	
	# START WITH FRAME 
	plot(frame.aea, col="white", main=paste0(moddy[i], " global"))
	#plot(dem.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE,)

	plot(Und, add=TRUE, breaks=breakUnd, legend=FALSE, col=plotcolUnd)
	plot(above, add=TRUE, breaks=breakAbove, legend=FALSE, col=plotcolAbove)

	# ADD STATES, GRIDS AND BACKGROUND PRES 
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(backup, pch=19, col="black", cex=0.2)
	plot(test.aea, add=TRUE, border="purple", lwd=1.5)
	plot(preds, add=TRUE, pch=17, cex=0.5, col=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	
	}
	
	setwd(path.fig)
	dev.off()
	
	
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
# MAKE FIGURE OF ENSEMBLE MAP OF PREDICTIONS FOR MAIN TEXT

for(i in 1:5){
	# load in saved rasters 
	setwd("C:/Users/DW/Desktop/temp.sept.30/R objects/model_predictions/TransplantFig")
	big.aea <- raster(paste0(moddy[i], ".full.tif"))
	pr1 <- raster(paste0(moddy[i], ".zoom.tif"))
	assign(paste0(moddy[i], "zoom"), pr1)
	assign(paste0(moddy[i], "big"), big.aea)
	rm(pr1); rm(big.aea)
	setwd(path.fig)
}
	big.aea <- mean(GLMbig, GAMbig, RFbig,  BRTbig, MAXbig, na.rm = TRUE)
	pr1 <- mean(GLMzoom, GAMzoom, RFzoom,  BRTzoom, MAXzoom, na.rm = TRUE)
	
	# get thereshold 
	myocc <- read.csv(file="C:/Users/DW/Desktop/temp.sept.30/data files/site_preds.csv")
	myocc <- myocc[,c("ID2", "LRavg", "GAMavg", "RFavg", "BRTavg", "MAXavg")]
	for(i in 1:dim(myocc)[1]){
		myocc$EN[i] <- mean(c(myocc$LRavg[i], myocc$GAMavg[i], myocc$RFavg[i], myocc$BRTavg[i], myocc$MAXavg[i]), na.rm=TRUE)
		}
	occ4 <- merge(occ, myocc, by.x="ID", by.y="ID2", all.x=TRUE, all.y=FALSE)
	library(SDMTools)
	thresh <- optim.thresh(occ4$PRESABS, occ4$EN, threshold = 101)
	cuts <- as.numeric(thresh[5])
		#$`max.sensitivity+specificity`
		#[1] 0.56

# COLOR BREAK POINTS 
	breakUnd <- c(seq(0, cuts*1000))
	breakAbove <- c(seq(cuts*1000, 1000))
	colorUnd =colorRampPalette(c("blue", "green", "yellow"))
	colorAbove =colorRampPalette(c("orange", "red"))
	plotcolUnd <- colorUnd(length(breakUnd))
	plotcolAbove <- colorAbove(length(breakAbove))
	plotcolAbove <- add.alpha(plotcolAbove, alpha=0.6)
	plotcolUnd <- add.alpha(plotcolUnd, alpha=0.6)

###########################################################################	
# Make ensemble map 
	setwd(path.fig); par(new = TRUE)
	pdf(file = "2.4.ENM_ensemble maps_main_text.pdf", width=(8.5 - (1.25 + 0.75)), height=(4.5), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	
	par(mfrow=c(1,2))
	# TRIM DOWN AND REPROJECT TO REGION 
	plot(test.aea, border="purple")	
# plot on raster under 
	Und <- pr1
	Und[Und > (cuts*1000)] <- NA
	Und[Und < 1] <- NA
	Und <- crop(Und, test.aea)
	plot(Und, add=TRUE, breaks=breakUnd, legend=FALSE, col=plotcolUnd)
# plot on raster over
	above <- pr1
	above[above < (cuts*1000)] <- NA
	above <- crop(above, test.aea)
	plot(above, add=TRUE, breaks=breakAbove, legend=FALSE, col=plotcolAbove)
	
	#plot(dem_zoom.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE)
	plot(ecoreg.aea2, add=TRUE, lwd=1.7)
	plot(preds, add=TRUE, pch=24, col="lightgrey", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	points(backup2, pch=21, col="lightgrey", bg="black", cex=0.8)
	verts = gridlines(ecoreg.wgs, easts=c(seq(-124.75, -121.75, by=0.5)), norths=c(seq(42.75, 45.75, by=0.5)))
	plot(gridlittle.aea, add=TRUE, lty=2, col="darkgrey")
	points(preds, cex=1.5, pch=24, col="black", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	lines(myline, lwd=4)
	text(-2252695, 2630000, "120 km", srt=90)
	plot(test.aea, border="purple", add=TRUE, lwd=1.7)
	# distance = 119km 
	par(xpd=TRUE)
	#legend("bottom", inset=-0.072, legend=c("within","beyond", "wild population"), pch=c(24, 24, 21), pt.bg=c("red","blue", "black"), bg="white", border="black")
	
	# add on legend
	#plot(pr1, legend.only=TRUE, legend.width=2, legend.shrink=0.75,
	#	breaks=c(seq(0,1000)), col=c(plotcolUnd, plotcolAbove), 
	#	axis.args=list(at=seq(0, 1000, by=100), labels=seq(0, 1, by=0.1), cex.axis=0.8),
	#	legend.args=list(text='Predicted Suitability', side=4, font=2, line=1.5, cex=0.8))	
	plot(test.aea, add=TRUE, border="purple", lwd=1.5)
	
# MAKE BIG MAP
	plot(frame.aea)
	Und <- big.aea
	Und[Und > (cuts*1000)] <- NA
	Und[Und < 1] <- NA
	above <- big.aea
	above[above < (cuts*1000)] <- NA
	plot(Und, add=TRUE, breaks=breakUnd, legend=FALSE, col=plotcolUnd)
	plot(above, add=TRUE, breaks=breakAbove, legend=FALSE, col=plotcolAbove)
	#plot(dem.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE,)
	plot(sta.aea2, add=T, border="black", lwd=0.2)
	plot(grd.aea2, add=T, lty="dashed", lwd=1)
	points(pres, pch=19, col="black", cex=0.3)
	backup <- pres
	text(coord.cull, labels=lab.cull, pos=c(1,1,1,1,1,1,4,4,4), offset=0, col="black") #this works for single panels, but not multi-panel
	plot(preds, add=TRUE, pch=24, cex=0.7, col="black", bg=c("blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"))
	plot(test.aea, add=TRUE, border="purple", lwd=1.5)
	
	# add on legend
	plot(big.aea, legend.only=TRUE, legend.width=1, legend.shrink=0.75,
		breaks=c(seq(0,1000)), col=c(plotcolUnd, plotcolAbove), 
		smallplot=c(0,.09, .3,.75),
		axis.args=list(at=seq(0, 1000, by=100), labels=seq(0, 1, by=0.1), cex.axis=0.8),
		legend.args=list(text='Predicted Suitability', side=4, font=2, line=1.8, cex=0.8))
	# labels to match original model output 0 -1 (rather th
	dev.off()
	

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	# RESIDUAL SCORES FOR MODELS
	WhichAvg <- c("LRavg", "GAMavg", "RFavg", "BRTavg", "MAXavg")
	# RESIDUAL SCORES FOR MODELS
	setwd(path.dat)
	# get values 
	allpres <- read.csv("site_preds.csv")
	allpres = allpres[allpres$ID1=='herb' | allpres$ID1 =='occ',]
	# add on ensemble scores
	justThese <- allpres[,c("LRavg", "GAMavg", "RFavg", "BRTavg", "MAXavg")]
	EN <- rowMeans(justThese, na.rm = TRUE)
	allpres$EN <- EN
	
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	setwd(path.fig)
	setwd(path.fig); par(new = TRUE)
	pdf(file = "2.5.RESIDUALS_from_predictions.pdf", width=(8.5 - (1.25 + 0.75)), height=(9 - (0.75 + 0.75)), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(1,3,2,3), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(3,2))

	for(i in 1:5){
	# GET MODEL CUTS 
	setwd(path.obj)
	cuts <- get(load(cuty[i]))
	cuts <- cuts[cuts$thresh=="SensSpec",]
	cuts <- mean(cuts$threshold)

	# COLOR BREAK POINTS 
	breakUnd <- c(seq(0, cuts, by=0.001))
	breakAbove <- c(seq(cuts, 1, by=0.001))
	colorUnd =colorRampPalette(c("blue", "green", "yellow"))
	colorAbove =colorRampPalette(c("orange", "red"))
	plotcolUnd <- colorUnd(length(breakUnd))
	plotcolAbove <- colorAbove(length(breakAbove))
	
	allpresO <- allpres[which(allpres[,WhichAvg[i]] > cuts), ]
	allpresU <- allpres[which(allpres[,WhichAvg[i]] < cuts), ]
	allpresU$Col <- colorUnd(length(breakUnd))[as.numeric(cut(allpresU[,WhichAvg[i]], breaks = breakUnd))]
	allpresO$Col <- colorAbove(length(breakAbove))[as.numeric(cut(allpresO[,WhichAvg[i]], breaks = breakAbove))]

	# make spatial & plot
	coordinates(allpresO) = ~Longitude + Latitude
	projection(allpresO) = prj.wgs
	allpres.aeaO = spTransform(allpresO, CRS=CRS(prj.aea))
	coordinates(allpresU) = ~Longitude + Latitude
	projection(allpresU) = prj.wgs
	allpres.aeaU = spTransform(allpresU, CRS=CRS(prj.aea))

	plot(frame.aea, col="lightgrey", main=paste0(moddy[i], " residuals"))
	plot(dem.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE,)
	sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
	grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
	points(allpres.aeaU, pch=21, col="lightgrey", bg=allpres.aeaU$Col, cex=1.2)
	points(allpres.aeaO, pch=21, col="lightgrey", bg=allpres.aeaO$Col, cex=1.2)

	# add on legend
	plot(pr1, legend.only=TRUE, legend.width=2, legend.shrink=0.75, breaks=c(seq(0,1000)), col=c(plotcolUnd, plotcolAbove), 
	smallplot=c(0.70,0.75, 0.3,.8),
	axis.args=list(at=seq(0, 1000, by=100), labels=seq(0, 1, by=0.1), cex.axis=0.8),
	 legend.args=list(text='Predicted Suitability', side=4, font=2, line=2.5, cex=0.8))

	}
	# ENSEMBLE 
	cuts <- 0.56
		# COLOR BREAK POINTS 
		breakUnd <- c(seq(0, cuts, by=0.001))
		breakAbove <- c(seq(cuts, 1, by=0.001))
		colorUnd =colorRampPalette(c("blue", "green", "yellow"))
		colorAbove =colorRampPalette(c("orange", "red"))
		plotcolUnd <- colorUnd(length(breakUnd))
		plotcolAbove <- colorAbove(length(breakAbove))
		allpresO <- allpres[which(allpres$EN > cuts), ]
		allpresU <- allpres[which(allpres$EN < cuts), ]
		allpresU$Col <- colorUnd(length(breakUnd))[as.numeric(cut(allpresU[,c("EN")], breaks = breakUnd))]
		allpresO$Col <- colorAbove(length(breakAbove))[as.numeric(cut(allpresO[,c("EN")], breaks = breakAbove))]
		# make spatial & plot
		coordinates(allpresO) = ~Longitude + Latitude
		projection(allpresO) = prj.wgs
		allpres.aeaO = spTransform(allpresO, CRS=CRS(prj.aea))
		coordinates(allpresU) = ~Longitude + Latitude
		projection(allpresU) = prj.wgs
		allpres.aeaU = spTransform(allpresU, CRS=CRS(prj.aea))
		plot(frame.aea, col="lightgrey", main=paste0("ENSEMBLE", " residuals"))
		plot(dem.aea, col=palette(gray(seq(0.6,0.98,len=150))), add=TRUE, legend=FALSE,)
		sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
		grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
		points(allpres.aeaU, pch=21, col="lightgrey", bg=allpres.aeaU$Col, cex=1.2)
		points(allpres.aeaO, pch=21, col="lightgrey", bg=allpres.aeaO$Col, cex=1.2)
	
		# add on legend
		plot(pr1, legend.only=TRUE, legend.width=2, legend.shrink=0.75, breaks=c(seq(0,1000)), col=c(plotcolUnd, plotcolAbove), 
		smallplot=c(0.70,0.75, 0.3,.8),
		axis.args=list(at=seq(0, 1000, by=100), labels=seq(0, 1, by=0.1), cex.axis=0.8),
		 legend.args=list(text='Predicted Suitability', side=4, font=2, line=2.5, cex=0.8))
		# labels to m
	
	setwd(path.fig)
	dev.off()
	###
	