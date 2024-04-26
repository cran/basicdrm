

#' Selects a best-fitting Hill model given defaults
#'
#' Using the function [fitHillModel()], this function fits four Hill models
#' with minimal and maximal effects either varying or fixed at the given
#' default values; it then selects the best fitting model based on the Bayesian
#' information criterio or Akaike information criterion, and returns a Hill fit
#' object with information from all fits included.
#'
#' @param formula Either an object of class `formula` such as would be provided
#' to a modeling function like [lm()], or a numeric vector of concentration
#' values (including 0 or Inf)
#' @param data If `forumula` is a symbolic formula, a data frame containing the
#' specified values. If `formula` is a numeric vector of concentrations, a
#' numeric vector of response values
#' @param defaults A two value numeric vector containing the default minimal
#' effect and the default maximal effect, in that order
#' @param weights A vector of weights (between 0 and 1) the same length as
#' `conc` and `act` which determines the weight with which each measurement
#' will impact the the sum of squared errors.  Weights will be multiplied by
#' errors *before* squaring.  If `NULL` (the default) all weights will be set
#' to 1. Can be a numeric vector, or the name of a column in `data` if `formula`
#' is a symbolic formula
#' @param start A vector of four starting values for the Hill model to be fit.
#' Any values not being fit will be fixed at these starting values.  If left as
#' `NULL`, a starting vector will be estimated from the data.
#' @param direction Determines the possible directionality of the dose response
#' model.  If 0 (the default) no additional constraints are placed on the
#' parameters.  If greater than 0, the fitting will require that the maximal
#' effect is *greater* than the minimal effect.  If less than 0, the fitting
#' wll require tha the maximal effect is *less* than the minimal effect.
#' @param lower A length-four vector of lower bounds on the Hill parameter
#' values.  Any parameters for which you do not wish to specify a bound can be
#' set to `NA`.
#' @param upper A vector of upper bounds on the Hill parameter values.  Works
#' the same as parameter `lower`.
#' @param useBIC Determines the information criterion to be used.  If `TRUE`
#' (the default), uses the Bayesian information criterion.  If `FALSE`, uses
#' the Akaike information criterion
#'
#' @return An object of class `hillrm`.  Contains all of the values found in
#' any `hillrm` object (see [fitHillModel()]), as well as `allfits`, a named
#' list of lists containing the `coefficients` and `par`vectors for each of the
#' individual fits, as well as the Bayesian information criterion (`bic`) and
#' Akaike informtion criterion (`aic`) values for each fit.
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#'
#' hfit <- findBestHillModel(conc,response,defaults=c(0,100))

findBestHillModel <- function(formula,data,defaults,weights=NULL,start=NULL,
							  direction=0,lower=NULL,upper=NULL,useBIC=TRUE) UseMethod("findBestHillModel")

#' @export
findBestHillModel.formula <- function(formula,data,defaults,weights=NULL,start=NULL,
									  direction=0,lower=NULL,upper=NULL,useBIC=TRUE) {
	mf <- stats::model.frame(formula=formula, data=data)
	conc <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(conc,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { conc <- conc[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	hfit <- findBestHillModel.default(conc,act,defaults,weights,start,
									  direction,lower,upper,useBIC)
	hfit$call <- match.call()
	return(hfit)
}

#' @export
findBestHillModel.default <- function(formula,data,defaults,weights=NULL,start=NULL,
									  direction=0,lower=NULL,upper=NULL,useBIC=TRUE) {
	conc <- formula
	act <- data

	prel <- conc>0 & is.finite(conc)

	if (is.null(weights)) { weights <- rep(1,length(conc)) }
	else if (length(weights)==1) { weights <- rep(weights,length(conc)) }

	# Pick appropriate parameter bounds
	tlower <- c(exp(1.5*log(min(conc[prel]))-0.5*log(max(conc[prel]))), 0.1, -Inf, -Inf)
	tupper <- c(exp(1.5*log(max(conc[prel]))-0.5*log(min(conc[prel]))),  10,  Inf,  Inf)
	if (!is.null(lower)) {
		if (length(lower)==4) { tlower[!is.na(lower)] <- lower[!is.na(lower)] }
		else { stop("Length of parameter 'lower' must equal 4 or the number of free parameters specified by 'model'.") }
	}
	if (!is.null(upper)) {
		if (length(upper)==4) { tupper[!is.na(upper)] <- upper[!is.na(upper)] }
		else { stop("Length of parameter 'upper' must equal 4 or the number of free parameters specified by 'model'.") }
	}
	pbounds <- rbind(tlower,tupper)

	# Pick appropriate starting values
	if (is.null(start)) {
		lmf <- mean(log(conc[prel])*act[prel])-mean(log(conc[prel]))*mean(act[prel])
		if (lmf>0) { start <- c(exp(mean(log(conc[prel]))), 1, defaults) }
		else       { start <- c(exp(mean(log(conc[prel]))), 1, defaults) }

		start <- pmin(pmax(start,pbounds[,1]),pbounds[,2])
	} else {
		start[3:4] <- defaults
		if (is.null(lower)) { pbounds[1,] <- pmin(pbounds[1,],start) }
		if (is.null(upper)) { pbounds[2,] <- pmax(pbounds[2,],start) }
	}
	lower <- pbounds[1,]
	upper <- pbounds[2,]


	models <- list(c(1,2),c(1,2,3),c(1,2,4),c(1,2,3,4))
	allmodels <- list()
	allfits <- list()
	for (i in 1:5) {
		if (i<5) { allmodels[[i]] <- fitHillModel(conc,act,models[[i]],weights,start,direction,lower,upper) }
		else {
			mnv <- mean(act)
			allmodels[[i]] <- list(coefficients=c(start[1:2],mnv,mnv),par=mnv,conc=conc,act=act,residuals=act-mnv)
		}
		allfits[[i]] <- list(coefficients=allmodels[[i]]$coefficients,
							 par=allmodels[[i]]$par,
							 AIC=estimateAIC(weights*allmodels[[i]]$residuals,length(allmodels[[i]]$par)),
							 BIC=estimateBIC(weights*allmodels[[i]]$residuals,length(allmodels[[i]]$par)))
	}
	names(allfits) <- c("m2p","m3puc","m3plc","m4p","const")

	if (useBIC) { icv <- sapply(allfits,function(m) m$BIC) }
	else { icv <- sapply(allfits,function(m) m$AIC) }
	parv <- sapply(allfits,function(m) length(m$par))
	biv <- modelSelect(icv,parv)

	structure(
		list(conc=conc,act=act,weights=weights,par=allmodels[[biv]]$par,coefficients=allmodels[[biv]]$coefficients,
				 fitted.values=allmodels[[biv]]$fitted.values,residuals=allmodels[[biv]]$residuals,
				 mname=allmodels[[biv]]$mname,model=allmodels[[biv]]$model,start=start,direction=direction,
				 pbounds=rbind(lower,upper),allfits=allfits),
		class = "hillrm"
	)
}

#' Estimate Akaike Information Criterion
#'
#' @param resid A vector of residuals from a given fit
#' @param npar The number of paramters in the model
#'
#' @return The Akaike informtion criterion (AIC) value for the fit
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#'
#' hfit4p <- fitHillModel(conc,response,c(1,2,3,4),start=c(0.5,1,0,100))
#' hfit3p <- fitHillModel(conc,response,c(1,2,4),start=c(0.5,1,0,100))
#'
#' aic4p <- estimateAIC(residuals(hfit4p),4)
#' aic3p <- estimateAIC(residuals(hfit3p),3)
estimateAIC <- function(resid,npar) {
	df <- length(resid)
	ssErr <- sum(resid^2)
	lLik <- -(df/2) * (log(2*pi)+log(ssErr)-log(df)+1)
	k <- npar
	aic <- 2*k - 2*lLik + 2*k*(k+1)/(df-k-1)
	return(aic)
}

#' Estimate Bayesian Information Criterion
#'
#' @param resid A vector of residuals from a given fit
#' @param npar The number of paramters in the model
#'
#' @return The Bayesian informtion criterion (BIC) value for the fit
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#'
#' hfit4p <- fitHillModel(conc,response,c(1,2,3,4),start=c(0.5,1,0,100))
#' hfit3p <- fitHillModel(conc,response,c(1,2,4),start=c(0.5,1,0,100))
#'
#' aic4p <- estimateBIC(residuals(hfit4p),4)
#' aic3p <- estimateBIC(residuals(hfit3p),3)
estimateBIC <- function(resid,npar) {
	df <- length(resid)
	ssErr <- sum(resid^2)
	lLik <- -(df/2) * (log(2*pi)+log(ssErr)-log(df)+1)
	k <- npar
	bic <- log(df)*k - 2*lLik
	return(bic)
}

modelSelect <- function(modAIC, modK, ik=2) {
	evRatioCut <- 15
	bmi <- which.min(ifelse(modK==ik,modAIC,NA))
	while(any(exp(0.5*(modAIC[bmi]-modAIC[which(modK>modK[bmi])]))>evRatioCut)) {
		bmi <- which.min(ifelse(modK==modK[bmi]+1,modAIC,NA))
	}
	return(bmi)
}
