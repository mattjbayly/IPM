# Function list from R code supplimentary material in:
			# Merow, C., Latimer, A. M., Wilson, A. M., McMahon, S. M., Rebelo,
			# A. G. and Silander, J. A. (2014), On using integral projection
			# models to generate demographically driven predictions of species'
			# distributions: development and validation using sparse data. Ecography,
			# 37: 1167–1183. doi: 10.1111/ecog.00839


# stepDIC:
# function for backward stepwise selection with MCMCglmm objects. 
# if quadratic term is included all lower terms are included. 
# but linear terms can be selected without quadratic			
				stepDIC=function(form.best,data,keep=NULL,DIC.diff=3,...){
				# keep specifies some terms to keep for sure to speed it up
				cand=attr(terms(form.best),'term.labels')
				model.best=MCMCglmm( form.best, data=data, verbose=FALSE,...)
				DIC.best=model.best$DIC
				cat('full model dic: ',DIC.best[1],'############################\n')
				print(summary(model.best))
				delta.DIC=-Inf
				step.counter=1
				resp=as.character(model.best$Fixed$formula)[2]
				stored.summaries=list()
				# list the linear terms to help formulas below
				linear=cand[-grep('I',cand)]
				while(all(delta.DIC<DIC.diff)){
				cat('-------------------step ',step.counter,' -----------------------\n')
				cat('best DIC so far: ', DIC.best , ' \n')
				cand=attr(terms(form.best),'term.labels')
				DICs.cand=vector()
				models=list()
				if(is.null(keep)) which.to.test=(1:length(cand))
				if(!is.null(keep)) which.to.test=(1:length(cand))[-match(keep,cand)]
				stored.summaries[[step.counter]]=list()
				# ensure higher order terms are only in if the linear predcitor is
				temp.form=list()
				for(i in 1:length(which.to.test)){
				toss=grep(cand[which.to.test[i]],cand,fixed=T)
				temp.form[[i]]=as.formula(paste(resp,'~', paste(cand[-toss], collapse='+')))
				}
				for(i in 1:length(temp.form)){
				models[[i]]=MCMCglmm( temp.form[[i]], data=data, verbose=FALSE,...)
				DICs.cand[i]=models[[i]]$DIC
				cat('###########################','model',rep(i,20))
				print(summary(models[[i]]))
				stored.summaries[[step.counter]][[i]]=summary(models[[i]])
				for(i in 1:length(temp.form)){
				models[[i]]=MCMCglmm( temp.form[[i]], data=data, verbose=FALSE,...)
				DICs.cand[i]=models[[i]]$DIC
				cat('###########################','model',rep(i,20))
				print(summary(models[[i]]))
				stored.summaries[[step.counter]][[i]]=summary(models[[i]])
				}
				delta.DIC=apply((outer(DICs.cand,DIC.best,'-')),2,min)
				if(all(delta.DIC<DIC.diff)){
				which.best=which.min(DICs.cand- DIC.best[step.counter])
				form.best=models[[which.best]]$Fixed$formula
				DIC.best[step.counter+1]=DICs.cand[which.best]
				model.best=models[[which.best]]
				}
				step.counter=step.counter+1
				}
				cat('+++++++++++++++++++++++++++++ best model +++++++++++++++++++++++')
				print(summary(model.best))
				return(list(model.best=model.best,stored.summaries=stored.summaries))
				}
				}

				

# shrink.matrix: 				
# A function to reduce the size of large matrices to make storage more practical.				
				shrink.matrix=function(m,min=1e-5){
				m[m<min]=0
				m=as(m,'sparseMatrix')
				}


				
# predict.MCMCglmm.cm: 
#== doesn't work well with factors. only binary factors coded as 0,1 work.
#== note that random effects in new data need to be named as their columnname.#
#== if you want to predict with Random effects and use new data, you've got to
#== specify the values of the Random effects in the newdata, which is a little annoying.
#== currently only works for 1 Random effect
#== unfortunately, if you want to do transformations of variables, they need
#== seperate columns, as this function does not make them for you like most predict functions do.

			predict.MCMCglmm.cm=function (object, newdata = NULL, marginal = object$Random$formula,
			type = "response", interval = "none", level = 0.95,
			return.post.pred=FALSE,...)
			{
			#== doesn't work well with factors. only binary factors coded as 0,1 work.
			#== note that random effects in new data need to be named as their columnname.#
			#== if you want to predict with Random effects and use new data, you've got to
			#== specify the values of the Random effects in the newdata, which is a little annoying.
			#== currently only works for 1 Random effect
			#== unfortunately, if you want to do transformations of variables, they need
			#== seperate columns, as this function does not make them for you like most predict functions do.
			if (is.null(object$Random$nfl) == FALSE) {
			rcomponents <- (as.character(object$Random$formula)[2])
			mcomponents <- (as.character(marginal)[2])
			marginal <- rep(rep(as.numeric(rcomponents %in% mcomponents),
			object$Random$nrt), object$Random$nfl)
			st <- c(1, cumsum(rep(object$Random$nrl, object$Random$nfl)) + 1)
			st <- st[-length(st)]
			end <- cumsum(rep(object$Random$nrl, object$Random$nfl))
			comp <- rep(1:length(object$Random$nfl), object$Random$nfl)
			keep <- unlist(mapply(st[which(marginal == 0)], end[which(marginal ==
			0)], FUN = ":"))
			} else {
			keep <- NULL
			rand = NULL
			marginal=-1
			}
			object$Sol <- object$Sol[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), drop = FALSE]
			if ( is.null(newdata) ){
			W <- cBind(object$X, object$Z)
			W <- W[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), drop = FALSE]
			} else {
			X=matrix(NA,nrow(newdata),object$Fixed$nfl)
			temp.form=as.formula(paste('~',as.character(object$Fixed$formula)[3]))
			for(i in 1:nrow(newdata)){
			X[i,]=model.matrix(temp.form, data=newdata[i,])
			}
			if(marginal==0) {
			rand= paste(labels(terms(object$Random$formula)),
			1:object$Random$nrl,sep='.')
			} else {
			rand=NULL
			}
			W=cbind(X,as.matrix(newdata[,rand]))
			W=as(W,"dgCMatrix")
			}
			# only for prediction interval, and probably grabbing variance
			if ((type == "response" & any(object$family != "gaussian" &
			object$family != "cengaussian")) | interval == "prediction") {
			if(is.null(newdata)){
			vpred <- matrix(0, dim(object$X)[1],
			sum(object$Random$nfl[which(marginal ==
			1)]^2, na.rm = T) + sum(object$Residual$nfl^2))
			cnt <- 0
			if (any(marginal == 1)) {
			st <- st[which(marginal == 1)]
			end <- end[which(marginal == 1)]
			comp <- comp[which(marginal == 1)]
			for (i in 1:length(st)) {
			for (j in 1:length(st)) {
			if (comp[i] == comp[j]) {
			cnt <- cnt + 1
			vpred[, cnt] <- diag(object$Z[, st[i]:end[i]] %*%
			t(object$Z[, st[j]:end[j]]))
			}
			}
			}
			}
			comp <- rep(1:length(object$Residual$nfl), object$Residual$nfl)
			for (i in 1:length(comp)) {
			for (j in 1:length(comp)) {
			if (comp[i] == comp[j]) {
			cnt <- cnt + 1

			vpred[, cnt][which(object$error.term == i &
			object$error.term == j)] <- 1
			}
			}
			}
			} else {
			vpred=matrix(1,nrow(newdata),1)
			}
			if (is.null(object$Random$nfl) == FALSE) {
			keep <- which(rep(rep(as.numeric(rcomponents %in% mcomponents),
			object$Random$nrt), object$Random$nfl^2) == 1)
			} else {
			keep=NULL
			}
			keep <- c(keep, which(rep(rep(rep(1, length(object$Residual$nrt)),
			object$Residual$nrt), object$Residual$nfl^2) == 1))
			postvar <- t(apply(object$VCV[, keep, drop = FALSE],
			1, function(x) { (vpred %*% x)}))
			}
			# predictions made here
			post.pred <- t(apply(object$Sol, 1, function(x) { (W %*% x)@x }))
			# probably samples post pred with right variance to make prediction interval
			if (interval == "prediction") {
			post.pred <- matrix(rnorm(prod(dim(post.pred)), post.pred,
			sqrt(postvar)), dim(post.pred)[1], dim(post.pred)[2])
			}
			# changes prediction to response scale
			if (type == "response") {
			if (any(object$family %in% c("poisson", "cenpoisson",
			"multinomial", "categorical", "gaussian", "cengaussian",
			"ordinal") == FALSE)) {
			stop("sorry - prediction on data scale not implemented for this family")
			}
			if (any(object$family %in% c("poisson", "cenpoisson"))) {
			if(is.null(newdata)){
			keep <- which(object$family %in% c("poisson", "cenpoisson"))
			} else {
			keep <- 1:nrow(newdata)
			}
			if (interval == "prediction") {
			post.pred[, keep] <- exp(post.pred[, keep])
			} else {
			post.pred[, keep] <- exp(post.pred[, keep] +
			0.5 * postvar[,keep])
			}
			}
			if (any(object$family %in% c("multinomial", "categorical"))) {
			c2 <- (16 * sqrt(3)/(15 * pi))^2
			if(is.null(newdata)){
			keep <- which(object$family %in% c("multinomial", "categorical"))

			} else {
			keep <- 1:nrow(newdata)
			}
			if (interval == "prediction") {
			post.pred[, keep] <- plogis(post.pred[, keep])
			} else {
			post.pred[, keep] <- plogis(post.pred[, keep]/sqrt(1 +
			c2 * postvar[, keep]))
			}
			}
			if (any(object$family %in% c("ordinal"))) {
			if(is.null(newdata)){
			keep <- which(object$family %in% c("ordinal"))
			} else {
			keep <- 1:nrow(newdata)
			}
			CP <- cbind(-Inf, 0, object$CP, Inf)
			q <- matrix(0, dim(post.pred)[1], length(keep))
			if (interval == "prediction") {
			for (i in 2:(dim(CP)[2] - 1)) {
			q <- q + (pnorm(CP[, i + 1] - post.pred[, keep]) -
			pnorm(CP[, i] - post.pred[, keep])) * (i -
			1)
			}
			} else {
			for (i in 2:(dim(CP)[2] - 1)) {
			q <- q + (pnorm(CP[, i + 1] - post.pred[, keep],
			0, sqrt(postvar[, keep] + 1)) - pnorm(CP[,
			i] - post.pred[, keep], 0, sqrt(postvar[,
			keep] + 1))) * (i - 1)
			}
			}
			post.pred[, keep] <- q
			rm(q)
			}
			}
			# process for final output
			pred <- matrix(colMeans(post.pred), dim(post.pred)[2], 1)
			if (!interval == "none") {
			pred <- cbind(pred, coda::HPDinterval(mcmc(post.pred), prob = level))
			colnames(pred) <- c("fit", "lwr", "upr")
			}
			rownames(pred) <- 1:dim(pred)[1]
			if(return.post.pred) return(post.pred)
			if(!return.post.pred) return(pred)
			}


				
				
				