# EXTRA PLOTS TRANS _ SURVIVAL 

# SITE LEVEL PLOTS
par(mfrow=c(3,3)) # 9 SITES 
z <- dfclean$z # 
sites <- levels(dfclean$SiteID); sites

	for(i in 1:length(sites)){
		mydataframe <- dfclean[which(dfclean$SiteID == sites[i]), ]
			plot(Surv.mean ~ z.mean, data = surv.ps, pch = 19, cex = 0.7, xlim = plot.range, ylim = c(0,1), xlab  =expression("Size , "*italic(z)), ylab = "Probability of surviving", main=sites[i], type="n")
			points(mydataframe$z, jitter(mydataframe$Surv, amount=0.02), col="black", cex=0.6)
			
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
			mod.site <- glm(Surv ~ z  , family = binomial, data = mydataframe)
			lines(fitted(mod.site) ~ z, data = mydataframe, col = as.factor(mydataframe$Region), lwd=4)

			legend('topright', legend = levels(dfclean$Region), col = 1:3, cex = 0.8, pch = 1)
			text(2.5, 0.9, paste0("Int ", as.character(round(coef(mod.site), 2))[1]))
			text(2.5, 0.8, paste0("Slope ", as.character(round(coef(mod.site), 2))[2]))
			text(2.5, 0.7, paste0("N= ", length(mod.site$fitted.values)))

	}