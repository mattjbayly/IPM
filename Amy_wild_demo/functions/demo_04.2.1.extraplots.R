# Year Specific
			firstyr <- dfclean[which(dfclean$Year == "1011"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					surv.plot.data <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = firstyr)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,2))
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="Year = 1011")
			points(z, jitter(firstyr$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.surv.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.surv.glm), 2))[2]))
			
		secondyr <- dfclean[which(dfclean$Year == "1112"), ]
			z <- secondyr$z # stand along object 
				# Plotting frame 
					surv.plot.data <- within(secondyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = secondyr)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="Year = 1112")
			points(z, jitter(secondyr$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = secondyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.surv.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.surv.glm), 2))[2]))
			
			
# Region Specific
			firstyr <- dfclean[which(dfclean$Region == "N"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					surv.plot.data <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = firstyr)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			par(mfrow=c(1,3))
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="North")
			points(z, jitter(firstyr$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.surv.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.surv.glm), 2))[2]))
	##############################
			firstyr <- dfclean[which(dfclean$Region == "C"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					surv.plot.data <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = firstyr)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="Center")
			points(z, jitter(firstyr$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.surv.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.surv.glm), 2))[2]))
			
				##############################
			firstyr <- dfclean[which(dfclean$Region == "S"), ]
			z <- firstyr$z # stand along object 
				# Plotting frame 
					surv.plot.data <- within(firstyr, {
						z.quantiles <- quantile(z, probs=seq(0,1,length=16))# make up vals to plot qq of size range
						z.quantiles[1] <- z.quantiles[1]-0.4 			# make sure we include the smallest
						z.quantiles <- z.quantiles[!duplicated(z.quantiles)]
						z.classes <- cut(z, z.quantiles)					# add column to df with quantiles
						rm(z.quantiles)
					})
			# basic global survival model
			mod.surv.glm <- glm(Surv ~ z  , family = binomial, data = firstyr)
			summary(mod.surv.glm)

			# calculate size specific groupwise summary statistics (estimates across size range - 16 classes) 
			surv.ps <- summaryBy(z + Surv ~ z.classes, data = surv.plot.data)
			surv.ps # (actual values for plotting vs. predicted) 
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main="South")
			points(z, jitter(firstyr$Surv, amount=0.02), col="#0000FF0A", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			lines(fitted(mod.surv.glm) ~ z, data = firstyr, col = "red")
			add_panel_label(ltype="a") # A.) B.) C.) ... 
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.surv.glm), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.surv.glm), 2))[2]))
			

##############################
##############################
##############################
##############################
# site level plots

par(mfrow=c(3,3))

z <- dfclean$z
sites <- levels(dfclean$SiteID)

	for(i in 1:length(sites)){
		firstyr <- dfclean[which(dfclean$SiteID == sites[i]), ]
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main=sites[i], type="n")
			points(firstyr$z, jitter(firstyr$Surv, amount=0.02), col="black", cex=0.6)
			
			# Could be fitted line from MEM, but will do that later.
			mod.global <- glm(Surv ~ z  , family = binomial, data = dfclean)
			lines(fitted(mod.global) ~ z, data = dfclean, col = "lightgrey", lwd=3)
			
			# adjust colour palette for other sites
			palette(c("green", "blue", "red"))
			palette(c(adjustcolor(palette(), alpha.f = 0.2)))

			
			# plot all other sites in big plot 
				for(j in 1:length(sites)){
						temp <- dfclean[which(dfclean$SiteID == sites[j]), ]
						mod.temp <- glm(Surv ~ z  , family = binomial, data = temp)
						lines(fitted(mod.temp) ~ z, data = temp, col = as.factor(temp$Region), lwd=1.5)
					}
			
			palette(c("green", "blue", "red")) # change back colors 
			mod.site <- glm(Surv ~ z  , family = binomial, data = firstyr)
			lines(fitted(mod.site) ~ z, data = firstyr, col = as.factor(firstyr$Region), lwd=4)

			legend('topright', legend = levels(dfclean$Region), col = 1:3, cex = 0.8, pch = 1)
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.site), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.site), 2))[2]))
			text(2.5, 0.7, paste0("N= ", length(mod.site$fitted.values)))

	}
	
setwd(path.fig)
dev.off()
