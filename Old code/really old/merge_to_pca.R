
# Final Growth, Survivorship, Flower and Fecundity Model
Gr <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|Uplot), 
				REML=FALSE, data=d, na.action=na.omit)
				
Fec <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + (1|Uplot),
				family = poisson, data=d, na.action=na.omit)

Flr <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + (1|Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
				
Sur <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|Uplot), 
				family=binomial(link="logit"),
				na.action=na.omit, data=d)				
###################################################################

fixef(Sur)
#confint(Sur, level = 0.99)

moo <- ranef(Sur); moo <- moo$Uplot
moo <- ranef(Gr); moo <- moo$Uplot
moo <- ranef(Flr); moo <- moo$Uplot
moo <- ranef(Fec); moo <- moo$Uplot




enviro$Uplot2 <- paste(enviro$SITE, enviro$PLOT, sep="")





temp <- read.csv(file.choose())

head(temp)
test <- merge(enviro, temp, by.x="Uplot2", by.y="PLOT", all.x=TRUE)
dim(test); 
head(test)

# PCA OF CORRELATION MATRIX 

fit <- princomp(test[, 5:14], cor=TRUE) # correlation matrix rather than covariance matrix (variables different units)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
par(mfrow=c(1,1))
#plot(fit,type="lines") # scree plot 
#fit$scores # the principal components
#biplot(fit, expand=0.5)

# load custom PCA plotting function
setwd(path.funct)
# thanks to Boris!
source("my_Biplot.R") # from Boris https://stat.ethz.ch/pipermail/r-help//2014-April/374098.html
library(MASS)
fit <- prcomp(test[, 5:14], scale. = TRUE) # correlation matrix rather than covariance matrix (variables different units)
# scale. scales variables to have unit variance (mean = 0, stdev =1)
var_exp <- (fit$sdev)^2 / sum(fit$sdev^2) # Percent variance explained by each axis

myBiplot(fit, choices=1:2, type = "t", pch=20, col="red", col.arrows = "#FF6600")

# color by moisture score
r <- range(test[,15]); n <- 16; cv <- rev(topo.colors(n)) # make color range
c <- cv[cut(test[,15],n)]  
myBiplot(fit, choices=1:2, type = "p", pch=20, col=c, col.arrows = "#FF6600")

