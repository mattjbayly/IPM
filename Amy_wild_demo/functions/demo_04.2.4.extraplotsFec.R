# Year Specific
			firstyr <- Fec.plot.data[which(Fec.plot.data$Year == "1011"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataF <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = firstyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataF)
			Fec.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,3))
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range, xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="Year = 1011")
			points(z, firstyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
			
		secondyr <- Fec.plot.data[which(Fec.plot.data$Year == "1112"), ]
			z <- secondyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataS <- within(secondyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = secondyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataS)
			Fec.ps # (actual values for plotting vs. predicted) 
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="Year = 1112")
			points(z, secondyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = secondyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
	
################
# 2012 FRUIT DATA
		secondyr <- Fec.plot.data[which(Fec.plot.data$Year == "1213"), ]
			z <- secondyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataS <- within(secondyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = secondyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataS)
			Fec.ps # (actual values for plotting vs. predicted) 
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="Year = 2012")
			points(z, secondyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = secondyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
			
	
			
# Region Specific
			firstyr <- Fec.plot.data[which(Fec.plot.data$Region == "N"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = firstyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataR)
			Fec.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,3))
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="North")
			points(z, firstyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
	##############################
			firstyr <- Fec.plot.data[which(Fec.plot.data$Region == "C"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = firstyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataR)
			Fec.ps # (actual values for plotting vs. predicted) 
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="Center")
			points(z, firstyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
			
				##############################
			firstyr <- Fec.plot.data[which(Fec.plot.data$Region == "S"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Fec.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Fecival model
			mod.Fec.glm <- glm(Fec ~ z  , family = poisson, data = firstyr)
			summary(mod.Fec.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Fec.ps <- summaryBy(z + Fec ~ z.classes, data = Fec.plot.dataR)
			Fec.ps # (actual values for plotting vs. predicted) 
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main="South")
			points(z, firstyr$Fec, col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Fec.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 40, paste0("Int ", as.character(round(coef(mod.Fec.glm), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.Fec.glm), 2))[2]))
			

##############################
##############################
##############################
##############################
# site level plots

par(mfrow=c(3,3))

z <- Fec.plot.data$z
sites <- levels(Fec.plot.data$SiteID)

	for(i in 1:length(sites)){
		firstyr <- Fec.plot.data[which(Fec.plot.data$SiteID == sites[i]), ]
			plot(Fec.mean ~ z.mean, data = Fec.ps, pch = 19, cex = 0.7, xlim = plot.range,  xlab  =expression("Size , "*italic(z)), ylab = "Number of fruits", main=sites[i], type="n")
			points(firstyr$z, firstyr$Fec, col="black", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			mod.global <- glm(Fec ~ z  , family = poisson, data = Fec.plot.data)
			lines(fitted(mod.global) ~ z, data = Fec.plot.data, col = "lightgrey", lwd=3)
			
			# adjust colour palette for other sites
			palette(c("green", "blue", "red"))
			palette(c(adjustcolor(palette(), alpha.f = 0.2)))

			
			# plot all other sites in big plot 
				for(j in 1:length(sites)){
						temp <- Fec.plot.data[which(Fec.plot.data$SiteID == sites[j]), ]
						mod.temp <- glm(Fec ~ z  , family = poisson, data = temp)
						lines(fitted(mod.temp) ~ z, data = temp, col = as.factor(temp$Region), lwd=1.5)
					}
			
			palette(c("green", "blue", "red")) # change back colors 
			mod.site <- glm(Fec ~ z  , family = poisson, data = firstyr)
			lines(fitted(mod.site) ~ z, data = firstyr, col = as.factor(firstyr$Region), lwd=4)

			legend('topright', legend = levels(Fec.plot.data$Region), col = 1:3, cex = 20, pch = 1)
			text(2.5, 60, paste0("Int ", as.character(round(coef(mod.site), 2))[1]))
			text(2.5, 20, paste0("Slope ", as.character(round(coef(mod.site), 2))[2]))
			text(2.5, 100, paste0("N= ", length(mod.site$fitted.values)))

	}
	
setwd(path.fig)
dev.off()
