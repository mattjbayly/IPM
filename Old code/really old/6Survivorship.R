#title: "SURVIVORSHIP"
#author: "Matthew Bayly"
#date: "Tuesday, November 25, 2014"

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

d$surv2 <- d$FALL.total.height # assuming only AUGUST census
d$surv2[d$surv2>0] <- 1 # Make 0/1

d$surv1 <- d$MID.total.height # assuming only JULY census
d$surv1[d$surv1>0] <- 1 # Make 0/1

d$surv0 <- d$post.trans.condition # assuming only JULY census
d$surv0 <- sub("A", "1", d$surv0)
d$surv0 <- sub("D", "0", d$surv0)
d$surv0 <- as.numeric(d$surv0)


d$fec.seed <- d$FALL_fruit # sum of fruits

#fec.seed <- rowSums(d[, c("FALL.FR.devo.", "FALL.FR.full.", "FALL.FR.dehis.")]) # adds fall fruit columns
fec.flower <- rowSums(d[, c("FALL.plus.devo", "FALL.flower.fail")]) # adds fall flower columns

d <- cbind(d, fec.flower) # bring them together
site <- levels(d$site); site <- as.factor(site) # slit sites
# filter down dataframe
d <- d[ ,c('size', 'sizeNext', 'surv2', 'surv1', 'surv0', 'fec.seed', 'fec.flower', 'site', 'plot', 'source','MID.total.height', 'post.trans.condition', 'pot')]
d$MID.total.height[d$MID.total.height<1] <- NA # Make NA if die in spring too 

### END DATA INITIALIZATION

###########################
###########################
###########################
###########################
###########################
###########################

### STAR INDIVIDUAL TRACKING PLOT
library(reshape)
id <- seq(1:nrow(d))# unique plant id
d <- cbind(d, id)
d$site <- factor(d$site, levels = c("THOMAS", "HUNTER", "WILEY", "CALAPOOIA", "MOSBY", "COAST", "ROCK", "LOOK"))
site <- levels(d$site); site <- as.factor(site) # slit sites
# melt to planitng & monitoring time intervals as line
mdata <- melt(d, id=c("size", "sizeNext", "fec.seed", "fec.flower", "site", "plot", "source", "MID.total.height", "post.trans.condition", "pot", "id")) 

# set ploting frame & id column for quick plant reference in plotting 
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
id <- as.factor(mdata$id)

###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
###########################
##################

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

##########################################
##########################################
##########################################
##########################################

params=data.frame(
surv.int=NA, # Intercept from logistic regression of survival
surv.slope=NA, # Slope from logistic regression of survival
growth.int=NA, # Intercept from linear regression of growth
growth.slope=NA, # Slope from linear regression of growth
growth.sd=NA, # Residual sd from the linear regression of growth
seed.int=NA, # Intercept from Poisson regression of seed number
seed.slope=NA, # Slope from Poisson regression of seed number
recruit.size.mean=NA, # Mean recruit size
recruit.size.sd=NA, # Standard deviation of recruit size
establishment.prob=NA # Probability of establishment
)

# make sure none are missing
d2 <- d2[complete.cases(d2[,c("potential")]),]
d2 <- d2[order(d2$potential),]
d2_post <- d2[complete.cases(d2[,c("surv0")]),]

# START PROBABILITY SURVIVE GIVEN STARTING POTENTIAL
setwd(path.funct)
source("modforms.R")
library("MASS")

# 1. logistic regression survivorship (GLOBAL)
par(mfrow=c(2,2))
surv.reg0 = glm(d2_post$surv0 ~ d2_post$potential, family=binomial(logit), data=d2_post)
surv.reg1 = glm(d2$surv1 ~ d2$potential, family=binomial(logit), data=d2)
surv.reg2 = glm(d2$surv2 ~ d2$potential, family=binomial(logit), data=d2)
plot(jitter(d2$surv2, amount=0.05) ~ d2$potential, data=d2, xlab="Start potential", ylab="Survival to FALL")
lines(d2$potential, surv.reg2$fitted, type="l", col="red")
lines(d2$potential, surv.reg1$fitted, type="l", col="blue")
lines(d2_post$potential, surv.reg0$fitted, type="l", col="green")
params$surv.int=coefficients(surv.reg2)[1]
params$surv.slope=coefficients(surv.reg2)[2]

# 2. logistic regression growth (GLOBAL)
d2_sN <- d2[complete.cases(d2[,c("sizeNext")]),]
d2_sN <- d2_sN[order(d2_sN$potential),]
growth.reg=lm(sizeNext~potential,data=d2_sN)
growth.reg1=lm(MID.total.height~size,data=d2)
params$growth.int=coefficients(growth.reg)[1]
params$growth.slope=coefficients(growth.reg)[2]
params$growth.sd=sd(resid(growth.reg))

# 3. seeds: Poisson regression
seed.reg=glm(fec.seed~potential,data=d2,family=poisson())
params$seed.int=coefficients(seed.reg)[1]
params$seed.slope=coefficients(seed.reg)[2]

#### PLOTING
# Growth plot
plot(d2_sN$potential,d2_sN$sizeNext,xlab="Start potential", ylab="FALL Size")
lines(d2_sN$potential, growth.reg$fitted, type="l", col="red")
# Fecundity plot
plot(d2$potential,d2$fec.seed,xlab="Start potential",ylab="Fall Fruits")
lines(d2$potential, seed.reg$fitted, type="l", col="red")

##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

# SITE SPECIFIC REGRESSIONS FOR SURVIVORSHIP; GROWTH & FECUNDITY
# make sure none are missing
d2 <- d2[complete.cases(d2[,c("potential")]),]
d2 <- d2[order(d2$potential),]
d2_post <- d2[complete.cases(d2[,c("surv0")]),]
# START PROBABILITY SURVIVE GIVEN STARTING POTENTIAL ~ for each site
par(mfrow=c(2,2))

# SURVIVORSHIP PLOT
plot(jitter(d2$surv2, amount=0.05) ~ d2$potential, data=d2, xlab="Start potential", ylab="Survival", xlim=c(-10, 120))
lines(d2$potential, surv.reg2$fitted, type="l", col="black")
for (i in 1:length(site)){
   current <- d2[d2$site==site[i], ]
   current_post <- d2_post[d2_post$site==site[i], ]
		surv.reg0 = glm(current_post$surv0 ~ current_post$potential, family=binomial(logit), data=current_post)
		surv.reg1 = glm(current$surv1 ~ current$potential, family=binomial(logit), data=current)
		surv.reg2 = glm(current$surv2 ~ current$potential, family=binomial(logit), data=current)
			lines(current$potential, surv.reg2$fitted, type="l", col="red")
			#lines(current$potential, surv.reg1$fitted, type="l", col="blue")
			#lines(current_post$potential, surv.reg0$fitted, type="l", col="green")
		text(x=-5, y=logit(coefficients(surv.reg2)[1]), site[i])
		assign(paste("surv.reg0",site[i], sep=""), surv.reg0)
		assign(paste("surv.reg1",site[i], sep=""), surv.reg1)
		assign(paste("surv.reg2",site[i], sep=""), surv.reg2)			
  }

# FALL GROWTH LINE FOR EACH SITE
d2_sN <- d2[complete.cases(d2[,c("sizeNext")]),]
d2_sN <- d2_sN[order(d2_sN$potential),]
plot(d2_sN$potential,d2_sN$sizeNext,xlab="Start potential", ylab="FALL Size")
lines(d2_sN$potential, growth.reg$fitted, type="l", col="black")
for (i in 1:length(site)){
   current <- d2_sN[d2_sN$site==site[i], ]
		growth.reg=lm(sizeNext~potential,data=current)
			lines(current$potential, growth.reg$fitted, type="l", col="red")
		assign(paste("growth.reg",site[i], sep=""), growth.reg)
}

# FECUNDITY LINE FOR EACH SITE
plot(d2$potential,d2$fec.seed,xlab="Start potential",ylab="Fall Fruits")
lines(d2$potential, seed.reg$fitted, type="l", col="black")
for (i in 1:length(site)){
   current <- d2[d2$site==site[i], ]
		seed.reg=glm(fec.seed~potential,data=current,family=poisson())
			lines(current$potential, seed.reg$fitted, type="l", col="red")
		assign(paste("seed.reg",site[i], sep=""), seed.reg)
}


##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

###### INDIVIDUAL SITE FECUNDITY PLOTS	
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
for (i in 1:length(site)){
	current <- d2[d2$site==site[i], ]
	current_post <- d2_post[d2_post$site==site[i], ]
		plot(current$fec.seed ~ current$potential, data=current, 
		xlab="Start potential", ylab="Fall Fruits", main=paste(site[i], "FALL", sep="-"),
		pch=21, bg=as.factor(current$plot), xlim=c(0, 120), ylim=c(0,30))
		######
		 #3. fruits: Poisson regression
			seed.reg=glm(current$fec.seed~current$potential,data=current,family=poisson())
			lines(current$potential, seed.reg$fitted, type="l", col="red")
				for (j in 1:length(site)) { 
					seed.regTEMP = get((paste("seed.reg",site[j], sep="")))
					temp <- d2[d2$site==site[j], ]
					lines(temp$potential, seed.regTEMP$fitted, type="l", col="lightgrey", lwd=2)
				}
	#### MAIN SITE LINE
		lines(current$potential, seed.reg$fitted, type="l", col="red", lwd=2)
	}
	
########################################################################################
##### INDIVIDUAL SITE FECUNDITY	
par(mfrow=c(2,4),mar=c(4,4,2,1)) 
for (i in 1:length(site)){
	current <- d2[d2$site==site[i], ]
	current_post <- d2_post[d2_post$site==site[i], ]
		plot(jitter(current_post$surv0, amount=0.05) ~ current_post$potential, data=current_post, 
		xlab="Start potential", ylab="Survival", main=paste(site[i], "MAY", sep="-"),
		pch=21, bg=as.factor(current$plot), xlim=c(0, 120), ylim=c(-0.05, 1.1))
		######
		surv.reg0 = glm(current_post$surv0 ~ current_post$potential, family=binomial(logit), data=current_post)
		surv.reg1 = glm(current$surv1 ~ current$potential, family=binomial(logit), data=current)
		surv.reg2 = glm(current$surv2 ~ current$potential, family=binomial(logit), data=current)
			#lines(current$potential, surv.reg2$fitted, type="l", col="lightgrey")
			#lines(current$potential, surv.reg1$fitted, type="l", col="blue")
			lines(current_post$potential, surv.reg0$fitted, type="l", col="green")
				for (j in 1:length(site)) { 
					surv.regTEMP = get((paste("surv.reg0",site[j], sep="")))
					temp_pt <- d2_post[d2_post$site==site[j], ]
					lines(temp_pt$potential, surv.regTEMP$fitted, type="l", col="lightgrey", lwd=2)
				}
	#### MAIN SITE LINE
		lines(current_post$potential, surv.reg0$fitted, type="l", col="green", lwd=2)
	}
	
	