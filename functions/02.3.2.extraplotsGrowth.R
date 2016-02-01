myColours = c("green", "blue", "red")
palette(myColours) 
library(RColorBrewer)
palette(brewer.pal(9, "Set1"))
mycols <- adjustcolor(palette(), alpha.f = 0.7)
palette(mycols) 


# Region Specific
			firstyr <- dfclean[which(dfclean$Region == "Beyond"), ]
			z <- firstyr$z # stand along object 
			grow.plot.data <- subset(firstyr, !is.na(z1))
			# Growth regression, simple linear - will do mixed model  with regions later
			mod.grow <- lm(z1 ~ z, data = grow.plot.data)
				plot(z1 ~ z,
					 data = grow.plot.data,
					 xlim = plot.range, ylim = plot.range, pch = 19, cex = 0.7, col=1,
					 xlab = expression("Initial size, "*italic(z)),
					 ylab = expression("Final size, "*italic(z)*"'"),
					 main="Region specific")
			# fitted line from growth regression 
			text(8, 2, paste0("Int ", as.character(round(coef(mod.grow), 2))[1]), col=1)
			text(8, 1, paste0("Slope ", as.character(round(coef(mod.grow), 2))[2]), col=1)
			abline(mod.grow, col=1, lwd=2)
			abline(0, 1, lty=2) # line of zero growth (not positive or negative), useful for negative growth at large size classes. 
	####################################
		secondyr <- dfclean[which(dfclean$Region == "Within"), ]
			z <- secondyr$z # stand along object 
			grow.plot.data <- subset(secondyr, !is.na(z1))
			# Growth regression, simple linear - will do mixed model  with regions later
			mod.grow <- lm(z1 ~ z, data = grow.plot.data)
				points(z1 ~ z, data = grow.plot.data, pch = 19, cex = 0.7, col=2)
			# fitted line from growth regression 
			abline(mod.grow, col=2, lwd=2)
			abline(0, 1, lty=2) # line of zero growth (not positive or negative), useful for negative growth at large size classes. 
		text(2, 8, paste0("Int ", as.character(round(coef(mod.grow), 2))[1]), col=2)
		text(2, 7, paste0("Slope ", as.character(round(coef(mod.grow), 2))[2]), col=2)
		legend('topright', legend = levels(dfclean$Region), col = 1:2, cex = 0.8, bg="white", pch = 19)


		
	











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
			z <- firstyr$z # stand along object 
			grow.plot.data <- subset(firstyr, !is.na(z1))
			grow.plot.data.all <- subset(dfclean, !is.na(z1))

			# Growth regression, simple linear - will do mixed model  with regions later
			mod.grow <- lm(z1 ~ z, data = grow.plot.data)
				plot(z1 ~ z,
					 data = grow.plot.data,
					 xlim = plot.range, ylim = plot.range, pch = 19, cex = 0.7,
					 col=as.factor(grow.plot.data$PlotID),
					 xlab = expression("Initial size, "*italic(z)),
					 ylab = expression("Final size, "*italic(z)*"'"),
					 main=sites[i])
			######
			# Could be fitted line from MEM, but will do that later.
			mod.global <- lm(z1 ~ z, data = grow.plot.data.all)
			lines(fitted(mod.global) ~ z, data = grow.plot.data.all, col = "lightgrey", lwd=3)
			
			# adjust colour palette for other sites
			palette(brewer.pal(9, "Set1"))
			palette(c(adjustcolor(palette(), alpha.f = 0.2)))

			# plot all other sites in big plot 
				for(j in 1:length(sites)){
						temp <- grow.plot.data.all[which(grow.plot.data.all$SiteID == sites[j]), ]
						mod.temp <- lm(z1 ~ z, data = temp)
						lines(fitted(mod.temp) ~ z, data = temp, col = as.factor(temp$Region), lwd=1.5)
					}
			
			palette(brewer.pal(9, "Set1"))
			mod.site <- lm(z1 ~ z  , data = grow.plot.data)
			lines(fitted(mod.site) ~ z, data = grow.plot.data, col = as.factor(grow.plot.data$Region), lwd=4)

			legend('topright', legend = levels(dfclean$Region), col = 1:3, cex = 0.8, pch = 1)
			text(2, 8, paste0("Int ", as.character(round(coef(mod.site), 2))[1]))
			text(2, 6, paste0("Slope ", as.character(round(coef(mod.site), 2))[2]))
			text(2, 4, paste0("N= ", length(mod.site$fitted.values)))

	}
	
setwd(path.fig)
dev.off()
