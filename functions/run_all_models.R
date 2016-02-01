# run_all_models 
# Quick background script for running all models 

# SURVIVAL 
SurB <- glmmML(surv_end ~ start, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBP <- glmmML(surv_end ~ start + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBS <- glmmML(surv_end ~ start + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBW <- glmmML(surv_end ~ start + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBWP <- glmmML(surv_end ~ start + moist_score + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBWS <- glmmML(surv_end ~ start + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBSP <- glmmML(surv_end ~ start + moist_score + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
SurBPSW <- glmmML(surv_end ~ start + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	

# growth
GrB <- lmer(size2_ln ~ start + (1|Uplot), na.action=na.omit, REML=FALSE, data=d)
GrBP <- lmer(size2_ln ~ start + moist_score + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBS <- lmer(size2_ln ~ start + ENSEMBLE + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBW <- lmer(size2_ln ~ start + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWP <- lmer(size2_ln ~ start + moist_score + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWS <- lmer(size2_ln ~ start + ENSEMBLE + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBSP <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBPSW <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + site_type + (1|Uplot), REML=FALSE, na.action=na.omit, data=d)	

# pFlower
FlrB <- glmmML(pFlower ~ size2_ln, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBP <- glmmML(pFlower ~ size2_ln + moist_score, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBS <- glmmML(pFlower ~ size2_ln + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBW <- glmmML(pFlower ~ size2_ln + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBWP <- glmmML(pFlower ~ size2_ln + moist_score + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBWS <- glmmML(pFlower ~ size2_ln + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBSP <- glmmML(pFlower ~ size2_ln + moist_score + ENSEMBLE, cluster=Uplot, family=binomial, na.action=na.omit, data=d)
FlrBPSW <- glmmML(pFlower ~ size2_ln + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=binomial, na.action=na.omit, data=d)	

# FECUNDITy
FecB <- glmmML(fec ~ size2_ln, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBP <- glmmML(fec ~ size2_ln + moist_score, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBS <- glmmML(fec ~ size2_ln + ENSEMBLE, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBW <- glmmML(fec ~ size2_ln + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBWP <- glmmML(fec ~ size2_ln + moist_score + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBWS <- glmmML(fec ~ size2_ln + ENSEMBLE + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBSP <- glmmML(fec ~ size2_ln + moist_score + ENSEMBLE, cluster=Uplot, family=poisson, na.action=na.omit, data=d)
FecBPSW <- glmmML(fec ~ size2_ln + moist_score + ENSEMBLE + site_type, cluster=Uplot, family=poisson, na.action=na.omit, data=d)	

###########################################################################################

# SURVIVAL 
SurB2 <- glmer(surv_end ~ start + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBP2 <- glmer(surv_end ~ start + moist_score + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBS2 <- glmer(surv_end ~ start + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBW2 <- glmer(surv_end ~ start + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBWP2 <- glmer(surv_end ~ start + moist_score + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBWS2 <- glmer(surv_end ~ start + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBSP2 <- glmer(surv_end ~ start + moist_score + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
SurBPSW2 <- glmer(surv_end ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	

# pFLOWER 
FlrB2 <- glmer(pFlower ~ size2_ln + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBP2 <- glmer(pFlower ~ size2_ln + moist_score + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBS2 <- glmer(pFlower ~ size2_ln + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBW2 <- glmer(pFlower ~ size2_ln + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBWP2 <- glmer(pFlower ~ size2_ln + moist_score + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBWS2 <- glmer(pFlower ~ size2_ln + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBSP2 <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)
FlrBPSW2 <- glmer(pFlower ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=binomial(link="logit"), na.action=na.omit, data=d)	

# Fecundity 
FecB2 <- glmer(fec ~ size2_ln + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBP2 <- glmer(fec ~ size2_ln + moist_score + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBS2 <- glmer(fec ~ size2_ln + ENSEMBLE + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBW2 <- glmer(fec ~ size2_ln + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBWP2 <- glmer(fec ~ size2_ln + moist_score + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBWS2 <- glmer(fec ~ size2_ln + ENSEMBLE + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBSP2 <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)
FecBPSW2 <- glmer(fec ~ size2_ln + moist_score + ENSEMBLE + site_type + (1|site/Uplot), family=poisson, na.action=na.omit, data=d)	

# growth
GrB2 <- lmer(size2_ln ~ start + (1|site/Uplot), na.action=na.omit, REML=FALSE, data=d)
GrBP2 <- lmer(size2_ln ~ start + moist_score + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBS2 <- lmer(size2_ln ~ start + ENSEMBLE + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBW2 <- lmer(size2_ln ~ start + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWP2 <- lmer(size2_ln ~ start + moist_score + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBWS2 <- lmer(size2_ln ~ start + ENSEMBLE + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBSP2 <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)
GrBPSW2 <- lmer(size2_ln ~ start + moist_score + ENSEMBLE + site_type + (1|site/Uplot), REML=FALSE, na.action=na.omit, data=d)	
