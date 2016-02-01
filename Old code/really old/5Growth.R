#title: "Merow_tutorial"
#author: "Matthew Bayly"
#date: "Tuesday, November 25, 2014"
#**A simple IPM for a long-lived perennial plant**

# SOURCE CODE:
#Merow C, Dalgren JP, Metcalf CJE, Childs DZ, Evans MEK, Jongejans E, Record S, Rees M, Salguero-Gómez R, McMahon S (2014) Advancing population ecology with integral projection models: a practical guide. Methods in Ecology and Evolution 5(2): 99-110. http://dx.doi.org/10.1111/2041-210X.12146 

### _Set directories for computer_ ###
########################### 
## set directories
path.set="C:/Users/DW/Desktop/transplant_analysis/Rcode"
	
setwd(path.set)
source("01_SetDirectories.R")
# LIBRARIES
library(IPMpack)
########################### 
#Open 2014 plant datafile 
setwd(path.dat.raw)
plantdat <- read.csv(file="plantdata.csv")
dim(plantdat);
#colnames(plantdat)

###########################
###########################
# INDEX
#### GROWTH PLOTS				~ basic global plots for project across sites
#### SURVIVORSHIP PLOTS			~ pool plots & check plant survival across sites
#### SITE GROWTH PLOTS 
#### FANCY GROTH PLOTS (SITES)
		### LARGE PLANTS
		### SMALL PLANTS
		### LARGE PLANTS (PLOTS)
		### SMALL PLANTS (PLOTS)

#### BOXPLOTS PER PLOT / PER SITE
		### SPRING GROWTH
		### FALL GROWTH
		
#### GROWTH & SURVIOR CURVES 
		### VITAL RATE REGRESSION



###########################
###########################
###########################
### START DATA INITIALIZATION

# Set up Merow plant parameters according to Merow
# remove seeds rows  
plantdat <- plantdat[ !grepl("SEED", plantdat$pot) , ]
#plantdat <- plantdat[ !grepl("HERB", plantdat$mid.condition) , ]
summary(plantdat$pot) 

d <- plantdat # dataframe
d$size <- d$start.height # start height is hight in greenhouse for now
d$sizeNext <- d$FALL.total.height # mid summer height is 
d$sizeNext[d$sizeNext<1] <- NA # Make NA if die

d$surv <- d$FALL.total.height # assuming only july census
d$surv[d$surv>0] <- 1 # Make 0/1
d$fec.seed <- d$FALL_fruit

#fec.seed <- rowSums(d[, c("FALL.FR.devo.", "FALL.FR.full.", "FALL.FR.dehis.")]) # adds fall fruit columns
fec.flower <- rowSums(d[, c("FALL.plus.devo", "FALL.flower.fail")]) # adds fall flower columns

d <- cbind(d, fec.flower) # bring them together
site <- levels(d$site); site <- as.factor(site) # slit sites
# filter down dataframe
d <- d[ ,c('size', 'sizeNext', 'surv', 'fec.seed', 'fec.flower', 'site', 'plot', 'source','MID.total.height', 'post.trans.condition', 'pot')]
d$MID.total.height[d$MID.total.height<1] <- NA # Make NA if die in spring too 

### END DATA INITIALIZATION

###########################
###########################
###########################

### START BASIC LIFE HIST PLOTS
# plotting frame
par(mfrow=c(2,2),mar=c(4,4,2,1))

plot(d$size,jitter(d$surv),xlab="Greenhouse size", ylab="Survival to FALL", col=as.factor(d$site)) # SURVIVAL
plot(d$size,d$sizeNext,xlab="Greenhouse size",ylab="Size FALL", col=as.factor(d$site)) # GROWTH
plot(d$size,d$fec.seed,xlab="Greenhouse size",ylab="Fruits", col=as.factor(d$site)) # FECUNDITY
#hist(d$sizeNext[is.na(d$size)]) # seedlings (not possible with data)
############

#### GROWTH TO FALL ACROSS SITES
par(mfrow=c(2,4),mar=c(4,4,2,1))
for (i in 1:length(site)){
   current <- d[d$site==site[i], ]
  plot(current$size, current$sizeNext, xlab="Greenhouse", ylab="FALL size", 
          col=as.factor(current$source), main=site[i], ylim=c(0, 390), xlim=c(0,250), pch=19, 
          cex=1.2); abline(0, 1)
  }

#### GROWTH TO SPRING ACROSS SITES
par(mfrow=c(2,4),mar=c(4,4,2,1))
for (i in 1:length(site)){
   current <- d[d$site==site[i], ]
  plot(current$size, current$MID.total.height, xlab="Greenhouse", ylab="SPRING size", 
          col=as.factor(current$source), main=site[i], ylim=c(0, 390), xlim=c(0,250), pch=19, 
          cex=1.2); abline(0, 1)
  }
### END BASIC LIFE HISTORY PLOTS

###########################
###########################
###########################

### STAR INDIVIDUAL TRACKING PLOT
library(reshape)
id <- seq(1:nrow(d))# unique plant id
d <- cbind(d, id)
d$begin <-  runif(nrow(d), 0, 4)
# order sites based on location rather than alphabetical
d$site <- factor(d$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
site <- levels(d$site); site <- as.factor(site) # slit sites
# melt to planitng & monitoring time intervals as line
mdata <- melt(d, id=c("size", "surv", "fec.seed", "fec.flower", "plot", "site", "source", "id", "post.trans.condition", "pot"))
dead_planting <- subset(mdata, mdata$post.trans.condition=="D")  # plants dead ~ 1 week after transplanting
dead_planting$value <- dead_planting$value - 10 # just to put it below the line
# set ploting frame & id column for quick plant reference in plotting 
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
id <- as.factor(mdata$id)

# SUPER PLOT FOR LARGE PLANTS (INDIVIDUALS)
#for (i in 1:length(site)){
   current <- mdata[mdata$site==site[i], ]
	plot(current$size, current$value, type="n", xlab="Greenhouse", ylab="GROWTH", main=site[i], ylim=c(-15, 390), xlim=c(0, 390)) 
		points(dead_planting$size, dead_planting$value, pch=21, bg="lightgrey", cex = 1.4) # dead transplants
		#each plant has its own line
		for (j in 1:length(id)) { 
			plant <- current[current$id==id[j], ]
			plant2 <- subset(plant, plant$variable=="MID.total.height" | plant$variable=="begin")  # RED just for spring growth
				lines(plant2$size, plant2$value, col="red", type="b", lwd=3) 
				lines(plant$size, plant$value, type="b", lwd=1.5, 
						bg=as.factor(plant$variable), pch=21, cex = 1.4) 				

  }
  }
  
####
# SUPER PLOT FOR SMALL PLANTS (INDIVIDUALS)
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
#for (i in 1:length(site)){
   current <- mdata[mdata$site==site[i], ]
	plot(current$size, current$value, type="n", xlab="Greenhouse", ylab="GROWTH", main=site[i], ylim=c(-15, 100), xlim=c(0, 100)) 
		points(dead_planting$size, dead_planting$value, pch=21, bg="lightgrey", cex = 1.4) # dead transplants
		#each plant has its own line
		for (j in 1:length(id)) { 
			plant <- current[current$id==id[j], ]
			plant2 <- subset(plant, plant$variable=="MID.total.height" | plant$variable=="begin")  # RED just for spring growth
				lines(plant2$size, plant2$value, col="red", type="b", lwd=3) 
				lines(plant$size, plant$value, type="b", lwd=1.5, 
						bg=as.factor(plant$variable), pch=21, cex = 1.4) 				

  }
  }
#### END INDIVIDUAL TRACKING PLOT

####
# SUPER PLOT FOR LARGE PLANTS (PLOTS)
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
#for (i in 1:length(site)){
   current <- mdata[mdata$site==site[i], ]
	plot(current$size, current$value, type="n", xlab="Greenhouse", ylab="GROWTH", main=site[i], ylim=c(-15, 390), xlim=c(0, 390)) 
		points(dead_planting$size, dead_planting$value, pch=21, bg=as.factor(dead_planting$plot), cex = 1.4) # dead transplants
		#each plant has its own line
		for (j in 1:length(id)) { 
			plant <- current[current$id==id[j], ]
			plant2 <- subset(plant, plant$variable=="MID.total.height" | plant$variable=="begin")  # RED just for spring growth
				lines(plant2$size, plant2$value, col="red", type="b", lwd=3) 
				lines(plant$size, plant$value, type="b", lwd=1.5, 
						bg=as.factor(plant$plot), pch=21, cex = 1.4) 				

  }
  }

####
# SUPER PLOT FOR SMALL PLANTS (PLOTS)
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
#for (i in 1:length(site)){
   current <- mdata[mdata$site==site[i], ]
	plot(current$size, current$value, type="n", xlab="Greenhouse", ylab="GROWTH", main=site[i], ylim=c(-15, 100), xlim=c(0, 100)) 
		points(dead_planting$size, dead_planting$value, pch=21, bg=as.factor(dead_planting$plot), cex = 1.4) # dead transplants
		#each plant has its own line
		for (j in 1:length(id)) { 
			plant <- current[current$id==id[j], ]
			plant2 <- subset(plant, plant$variable=="MID.total.height" | plant$variable=="begin")  # RED just for spring growth
				lines(plant2$size, plant2$value, col="red", type="b", lwd=3) 
				lines(plant$size, plant$value, type="b", lwd=1.5, 
						bg=as.factor(plant$plot), pch=21, cex = 1.4) 				

  }
  }

  #### END INDIVIDUAL TRACKING PLOT

###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################

#### START BARPLOT OF GROWTH
# Which starting variable correlates best with growth & survival ( ~ a possible starting potential?)
par(mfrow=c(1,2))
plot(d$size, d$sizeNext, pch=19, cex=0.5, ylim=c(0, 390), xlim=c(0, 390), xlab="Greenhouse Height", ylab="Fall Height")
d$pot <- factor(d$pot, levels = c("XS", "S", "M", "L", "XL"))
boxplot(d$sizeNext~d$pot, col="darkgrey", xlab="Greenhouse AGE ~ POT", ylab="Fall Height")
###########################
# make model of growth given starting potential. 
# replace values
d$pot <- as.character(d$pot)
d$pot[d$pot=="XL"]<-5
d$pot[d$pot=="XS"]<-1
d$pot[d$pot=="S"]<-2
d$pot[d$pot=="L"]<-4 #LARGE & MEDIUM MIHT HAVE TO SWITCH
d$pot[d$pot=="M"]<-3
d$pot <- as.numeric(d$pot) # BACK TO NUMERIC
summary(d$pot) 
#######
# quick normality check 
d$ln_size <- log(d$size + 1)
hist(d$ln_size) # better 
z1 <- lm(d$MID.total.height ~ d$ln_size + d$pot)  # no interaction 
z2 <- lm(d$MID.total.height ~ d$ln_size + d$pot + d$ln_size:d$pot)  # interaction term present
z3 <- lm(d$MID.total.height ~ d$size + d$pot + d$size:d$pot) # untransformed
plot(d$MID.total.height ~ d$ln_size * d$pot)
plot(d$MID.total.height ~ d$size)
### EXAMIN
anova(z1) # no interaction
anova(z2) # interaction
####
summary(z1) # no interaction
summary(z2) # interaction
summary(z3) # untransformed
####
AIC(z1)  # no interaction
AIC(z2)  # interaction
AIC(z3)  # interaction

preddom <- d[ ,c('size', 'pot')]
potential <- predict(z2, preddom) # choose z2 because z3 violates normality assumptions
dim(d)
length(potential)
d2 <- cbind(d, potential)
plot(d2$MID.total.height ~ d2$potential)
d2$growth_given_potential_FALL <- d2$sizeNext/d2$potential # FALL growth given greenhouse starting potential
d2$growth_given_potential_SPRING <- d2$MID.total.height/d2$potential # SPRING growth given greenhouse starting potential

############### examin & log transform
par(mfrow=c(2,2),mar=c(4,4,2,1))
	hist(d2$growth_given_potential_SPRING); hist(d2$growth_given_potential_FALL)
d2$ln_growth_given_potential_FALL <- log(d2$growth_given_potential_FALL + 1) # FALL growth given greenhouse starting potential
d2$ln_growth_given_potential_SPRING <- log(d2$growth_given_potential_SPRING + 1) # SPRING growth given greenhouse starting potential
	hist(d2$ln_growth_given_potential_SPRING); hist(d2$ln_growth_given_potential_FALL)
###
###
#### AVERAGE DOWN TO PLOT LEVEL
# average down to plot level for sites
frame <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(frame) <- c("ln_growth_given_potential_FALL", "site")
for(i in 1:length(site)){
  current <- d2[which(d2$site==paste(site[i])), ]
  current2 <-  current[!is.na(current$ln_growth_given_potential_FALL),]
  final <- tapply(current2$ln_growth_given_potential_FALL, current2$plot, mean)
  final <- data.frame(final)
  final$siteR <- paste(site[i])
  colnames(final) <- c("ln_growth_given_potential_FALL", "site")
  frame <- rbind(frame, final); rm(final) 
}
frame_FALL <- frame
#### repeat for spring
frame <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(frame) <- c("ln_growth_given_potential_SPRING", "site")
for(i in 1:length(site)){
  current <- d2[which(d2$site==paste(site[i])), ]
  current2 <-  current[!is.na(current$ln_growth_given_potential_SPRING),]
  final <- tapply(current2$ln_growth_given_potential_SPRING, current2$plot, mean)
  final <- data.frame(final)
  final$siteR <- paste(site[i])
  colnames(final) <- c("ln_growth_given_potential_SPRING", "site")
  frame <- rbind(frame, final); rm(final) 
}
frame_SPRING <- frame

# BIND FRAMES TOGETHER
colnames(frame_SPRING) <- c("ln_growth_given_potential", "site")
frame_SPRING$season <- "SPRING"

colnames(frame_FALL) <- c("ln_growth_given_potential", "site")
frame_FALL$season <- "FALL"

frame <- rbind(frame_SPRING, frame_FALL)

#### GROWTH TO SPRING ACROSS SITES
par(mfrow=c(2,4),mar=c(4,4,2,1))
for (i in 1:length(site)){
      current <- frame[frame$site==site[i], ]
		boxplot(current$ln_growth_given_potential~current$season, data=current, 
			col=(c("azure4","coral1")), ylab="Relative growth", main=site[i], ylim=c(0, 3))
			
 stripchart(current$ln_growth_given_potential~current$season, data=current, 
            vertical = TRUE, method = "jitter", 
            pch = 21, col = "black", bg = "lightgrey", 
            add = TRUE, cex=1.2) 	
  }



###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
