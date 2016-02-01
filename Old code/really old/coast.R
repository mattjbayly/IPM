#title: "Polygons multiple points to coast"
#author: "Matthew Bayly"
#date: "Thursday, December 4, 2014"
#

library(sp)
library(rgdal)
library(dismo)
library(rgeos)
library(raster)

# define basic spatial projections (to maybe use later)
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# generate random coordinates for weather stations (but replace with real coordinates
# once everything is up and running
Latitude <- runif(200, 42, 45) # random points in PNW
Longitude <- runif(200, 120, 123) # random points in PNW
temp_id <- seq(1, 200, by=1)
stations <- data.frame(temp_id)
stations <- cbind(stations, Latitude, Longitude)
# replace the last 5 lines with your own data frame
# stations <- read.csv(file.choose()) # remove hashtag

# make weather stations into spatial points
coordinates(stations) = ~Longitude + Latitude
proj4string(stations) = CRS(prj.wgs)
# just a demonstration of how to reproject spatial data in R 
stations.aea = spTransform(stations, CRS=CRS(prj.aea)) # not used just demo 

# import the coast line polygon
?readOGR
# dsn = the root directory of the shape file
# layer = the layer name with no extensions
coastline = readOGR(dsn="C:/Users/DW/Desktop/temp.sept.30/R objects/ecoregions.shp", layer="us_eco_l3_no_st") 




###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
