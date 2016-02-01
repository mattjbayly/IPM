# Year Specific
			firstyr <- Repr.plot.data[which(Repr.plot.data$Year == "1011"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataF <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = firstyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataF)
			Repr.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,3))
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="Year = 1011")
			points(z, jitter(firstyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
			
		secondyr <- Repr.plot.data[which(Repr.plot.data$Year == "1112"), ]
			z <- secondyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataS <- within(secondyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = secondyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataS)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="Year = 1112")
			points(z, jitter(secondyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = secondyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
			
				
		secondyr <- Repr.plot.data[which(Repr.plot.data$Year == "1213"), ]
			z <- secondyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataS <- within(secondyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = secondyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataS)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="Year = 1213")
			points(z, jitter(secondyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = secondyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
			
					
			
			
# Region Specific
			firstyr <- Repr.plot.data[which(Repr.plot.data$Region == "N"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = firstyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataR)
			Repr.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,3))
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="North")
			points(z, jitter(firstyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
	##############################
			firstyr <- Repr.plot.data[which(Repr.plot.data$Region == "C"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = firstyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataR)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="Center")
			points(z, jitter(firstyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
			
				##############################
			firstyr <- Repr.plot.data[which(Repr.plot.data$Region == "S"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					Repr.plot.dataR <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global Reprival model
			mod.Repr.glm <- glm(Repr ~ z  , family = binomial, data = firstyr)
			summary(mod.Repr.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			Repr.ps <- summaryBy(z + Repr ~ z.classes, data = Repr.plot.dataR)
			Repr.ps # (actual values for plotting vs. predicted) 
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main="South")
			points(z, jitter(firstyr$Repr, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.Repr.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.Repr.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.Repr.glm), 2))[2]))
			

##############################
##############################
##############################
##############################
# site level plots

par(mfrow=c(3,3))

z <- Repr.plot.data$z
sites <- levels(Repr.plot.data$SiteID)

	for(i in 1:length(sites)){
		firstyr <- Repr.plot.data[which(Repr.plot.data$SiteID == sites[i]), ]
			plot(Repr.mean ~ z.mean, data = Repr.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Prob Reproductive", main=sites[i], type="n")
			points(firstyr$z, jitter(firstyr$Repr, amount=0.02), col="black", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			mod.global <- glm(Repr ~ z  , family = binomial, data = Repr.plot.data)
			lines(fitted(mod.global) ~ z, data = Repr.plot.data, col = "lightgrey", lwd=3)
			
			# adjust colour palette for other sites
			palette(c("green", "blue", "red"))
			palette(c(adjustcolor(palette(), alpha.f = 0.2)))

			
			# plot all other sites in big plot 
				for(j in 1:length(sites)){
						temp <- Repr.plot.data[which(Repr.plot.data$SiteID == sites[j]), ]
						mod.temp <- glm(Repr ~ z  , family = binomial, data = temp)
						lines(fitted(mod.temp) ~ z, data = temp, col = as.factor(temp$Region), lwd=1.5)
					}
			
			palette(c("green", "blue", "red")) # change back colors 
			mod.site <- glm(Repr ~ z  , family = binomial, data = firstyr)
			lines(fitted(mod.site) ~ z, data = firstyr, col = as.factor(firstyr$Region), lwd=4)

			legend('topright', legend = levels(Repr.plot.data$Region), col = 1:3, cex = 0.8, pch = 1)
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.site), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.site), 2))[2]))
			text(2.5, 0.7, paste0("N= ", length(mod.site$fitted.values)))

	}
	
setwd(path.fig)
dev.off()
