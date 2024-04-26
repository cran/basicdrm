
#' Fit a Hill dose response model to data
#'
#' This function uses the R `stats` function `optim` to fit a Hill dose
#' response model to a given set of dose and response values.  Four different
#' model settings are allowed, in which the minimal and maximal effects are
#' either fixed at a provided value or allowed to be fit to the data.
#'
#' @param formula Either an object of class `formula` such as would be provided
#' to a modeling function like [lm()], or a numeric vector of concentration
#' values (including 0 or Inf)
#' @param data If `forumula` is a symbolic formula, a data frame containing the
#' specified values. If `formula` is a numeric vector of concentrations, a
#' numeric vector of response values
#' @param model A vector of values between 1 and 4, specifying the precise
#' model to be fit. The values correspond to the four parameters of the Hill
#' model: dose of median effect, Hill slope, minimal effect, and maximal effect
#' (see [evalHillModel()]).  The first of these two are always fit, so `model`
#' must contain at least `1` and `2`.  The presence of `3` or `4` will determine
#' if those parameters are also fit, or fixed at the given starting value. So
#' `c(1,2,4)` will fit the dose of median effect, the Hill slope, and the
#' maximal effect, but will leave the minimal effect fixed at the starting
#' value.
#' @param weights A vector of weights (between 0 and 1) the same length as
#' `conc` and `act` which determines the weight with which each measurement
#' will impact the the sum of squared errors.  Weights will be multiplied by
#' errors *before* squaring.  If `NULL` (the default) all weights will be set
#' to 1. Can be a numeric vector, or the name of a column in `data` if `formula`
#' is a symbolic formula
#' @param start A vector of four starting values for the Hill model to be fit.
#' Any values not being fit will be fixed at these starting values.  If left as
#' `NULL`, a starting vector will be estimated from the data, but it will almost
#' always be better to provide an explicit staring model.
#' @param direction Determines the possible directionality of the dose response
#' model.  If 0 (the default) no additional constraints are placed on the
#' parameters.  If greater than 0, the fitting will require that the maximal
#' effect is *greater* than the minimal effect.  If less than 0, the fitting
#' wll require tha the maximal effect is *less* than the minimal effect.
#' @param lower A vector of lower bounds on the Hill parameter values.  Can be
#' the same length as `model` (in which case the bounds will be applied to the
#' corresponding fit parameters) or the full length of 4.  Any parameters for
#' which you do not wish to specify a bound can be set to `NA`.
#' @param upper A vector of upper bounds on the Hill parameter values.  Works
#' the same as parameter `lower`.
#'
#' @return An object of class `hillrm`, containing the following values:
#' * `conc`: the given vector of concentraitons
#' * `act`: the given vector of responses
#' * `weights`: the vector of measurement weights used in minimizing the sum
#' of squared errors
#' * `coefficients`: the full four-parameter Hill parameter vector (accessible
#' by the function `coef()`)
#' * `par`: the vector of paramters that were actually fit
#' * `fitted.values`: the predicted responses of the best fit model (accessible
#' by the functoin `fitted()`)
#' * `residuals`: the difference between the actual responses and the predicted
#' responses (accessible by the function `residuals()`)
#' * `model`: the vector of values between 1 and 4 specifying the precise model
#' that was fit
#' * `mname`: a character string naming the precise model that was fit. One of
#' "m2p", "m3plc", "m3puc", or "m4p"
#' * `start`: a four-value parameter vector used as the starting value for the
#' model fit
#' * `direction`: the direction constraint used in the fit
#' * `pbounds`: a two-by-four matrix of values specifying the lower and upper
#' bounds used in the fit
#'
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#' data <- data.frame(conc=conc,response=response,weight=c(0.5,rep(1,10),0.1))
#'
#' hfit <- fitHillModel(conc,response,c(1,2,3,4),start=c(0.5,1,0,100))
#' hfit2 <- fitHillModel(response~conc,data,c(1,2,4),weight,start=c(0.5,1,0,100),
#'                       direction=0,lower=c(NA,NA,0,0))
fitHillModel <- function(formula,data,model,weights=NULL,start=NULL,
						 direction=0,lower=NULL,upper=NULL) UseMethod("fitHillModel")

#' @export
fitHillModel.formula <- function(formula,data,model,weights=NULL,start=NULL,
								 direction=0,lower=NULL,upper=NULL) {
	mf <- stats::model.frame(formula=formula, data=data)
	conc <- stats::model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(conc,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { conc <- conc[,-i] }
	}
	act <- stats::model.response(mf)
	weights <- eval(substitute(weights),data)
	hfit <- fitHillModel.default(conc,act,model,weights,start,
									  direction,lower,upper)
	hfit$call <- match.call()
	return(hfit)
}

#' @export
fitHillModel.default <- function(formula,data,model,weights=NULL,start=NULL,
								 direction=0,lower=NULL,upper=NULL) {
	conc <- formula
	act <- data

	prel <- conc>0 & is.finite(conc)
	model <- sort(unique(model))

	if (is.null(weights)) { weights <- rep(1,length(conc)) }
	else if (length(weights)==1) { weights <- rep(weights,length(conc)) }

	# Pick appropriate parameter bounds
	tlower <- c(exp(1.5*log(min(conc[prel]))-0.5*log(max(conc[prel]))), 0.1, -Inf, -Inf)
	tupper <- c(exp(1.5*log(max(conc[prel]))-0.5*log(min(conc[prel]))),  10,  Inf,  Inf)
	if (!is.null(lower)) {
		if (length(lower)==4) { tlower[!is.na(lower)] <- lower[!is.na(lower)] }
		else if (length(lower)==length(model)) { tlower[model[!is.na(lower)]] <- lower[!is.na(lower)] }
		else { stop("Length of parameter 'lower' must equal 4 or the number of free parameters specified by 'model'.") }
	}
	if (!is.null(upper)) {
		if (length(upper)==4) { tupper[!is.na(upper)] <- upper[!is.na(upper)] }
		else if (length(upper)==length(model)) { tupper[model[!is.na(upper)]] <- upper[!is.na(upper)] }
		else { stop("Length of parameter 'upper' must equal 4 or the number of free parameters specified by 'model'.") }
	}
	pbounds <- rbind(tlower,tupper)

	# Pick appropriate starting values
	if (is.null(start)) {
		lmf <- mean(log(conc[prel])*act[prel])-mean(log(conc[prel]))*mean(act[prel])
		if (lmf>0) { start <- c(exp(mean(log(conc[prel]))), 1, min(act), max(act)) }
		else       { start <- c(exp(mean(log(conc[prel]))), 1, max(act), min(act)) }

		start <- pmin(pmax(start,pbounds[,1]),pbounds[,2])
	} else {
		if (is.null(lower)) { pbounds[1,] <- pmin(pbounds[1,],start) }
		if (is.null(upper)) { pbounds[2,] <- pmax(pbounds[2,],start) }
	}

	if (!(1%in%model)||!(2%in%model)) {
		stop("For this function, valid models must allow dose of median effect ",
			 "(parameter 1) and Hill slope (parameter 2) to vary.")
	}
	if (3%in%model) {
		if (4%in%model) { fit <- fitHill4Par(conc,act,weights,start,direction,pbounds) }
		else { fit <- fitHill3ParU(conc,act,weights,start,direction,pbounds) }
	} else {
		if (4%in%model) { fit <- fitHill3ParL(conc,act,weights,start,direction,pbounds) }
		else { fit <- fitHill2Par(conc,act,weights,start,direction,pbounds) }
	}

	hcoef <- fit$fullpar
	names(hcoef) <- c("IDM","n","E0","Ef")
	hpar <- fit$par
	names(hpar) <- names(hcoef)[fit$model]
	fitv <- hcoef[["E0"]]+(hcoef[["Ef"]]-hcoef[["E0"]])*evalHillModel_sf(conc,hcoef[c("IDM","n")])
	hfit <- structure(
		list(conc=conc, act=act, weights=weights,
			 par=hpar, coefficients=hcoef,
			 fitted.values=fitv, residuals=(act-fitv),
			 mname=fit$mname, model=fit$model,
			 start=start, direction=direction, pbounds=pbounds),
		class="hillrm"
	)
	return(hfit)
}

fitHill2Par <- function(conc,act,weights,start,direction,pbounds) {
	if ((direction>0 && start[[3]]>start[[4]]) ||
			(direction<0 && start[[3]]<start[[4]])) {
		stop("Initial values do not satisfy constraints")
	}

	nbounds <- log(pbounds[,1:2])
	nstart <- log(start[1:2])

	valderivfunc <- function(parv) {
		sfpar <- exp(parv)
		sfres <- evalHillModel_sf(conc,sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		sfact <- start[[3]]+(start[[4]]-start[[3]])*sfact - act

		ovalue <- sum((weights*sfact)^2)
		derivs <- (2*(start[[4]]-start[[3]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)*sfpar
		return(list(value=ovalue,derivatives=derivs))

	}
	parv2fullpar <- function(parv) {
		return(c(exp(parv),start[3:4]))
	}

	nls <- runBoundedOptim(valderivfunc,parv2fullpar,nstart,nbounds)
	nls$model <- c(1,2)
	nls$par <- nls$fullpar[nls$model]
	nls$mname <- "m2p"
	return(nls)
}
fitHill3ParU <- function(conc,act,weights,start,direction,pbounds) {
	if (direction>0) { pbounds[2,3] <- min(pbounds[2,3],start[[4]]) }
	else if (direction<0) { pbounds[1,3] <- min(pbounds[1,3],start[[4]]) }

	nbounds <- log(pbounds[,1:2])
	nstart <- log(start[1:2])

	wt2 <- (weights^2)/mean(weights^2)

	valderivfunc <- function(parv) {
		sfpar <- exp(parv)
		sfres <- evalHillModel_sf(conc,sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*(1-sfact)),mean(wt2*((1-sfact)^2)),
				 mean(wt2*act),mean(wt2*(1-sfact)*act))
		ebnd <- boundedOpt1d(mnv,start[[4]],pbounds[,3])[[1]]
		sfact <- ebnd+(start[[4]]-ebnd)*sfact - act

		ovalue <- sum((weights*sfact)^2)
		derivs <- (2*(start[[4]]-ebnd))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)*sfpar
		return(list(value=ovalue,derivatives=derivs))

	}
	parv2fullpar <- function(parv) {
		sfpar <- exp(parv)
		sfact <- evalHillModel_sf(conc,sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*(1-sfact)),mean(wt2*((1-sfact)^2)),
				 mean(wt2*act),mean(wt2*(1-sfact)*act))
		ebnd <- boundedOpt1d(mnv,start[[4]],pbounds[,3])[[1]]
		hpar <- c(sfpar,ebnd,start[[4]])
		return(hpar)

	}

	nls <- runBoundedOptim(valderivfunc,parv2fullpar,nstart,nbounds)
	nls$model <- c(1,2,3)
	nls$par <- nls$fullpar[nls$model]
	nls$mname <- "m3puc"
	return(nls)
}
fitHill3ParL <- function(conc,act,weights,start,direction,pbounds) {
	if (direction>0) { pbounds[1,4] <- max(pbounds[1,4],start[[3]]) }
	else if (direction<0) { pbounds[2,4] <- min(pbounds[2,4],start[[3]]) }

	nbounds <- log(pbounds[,1:2])
	nstart <- log(start[1:2])

	wt2 <- (weights^2)/mean(weights^2)

	valderivfunc <- function(parv) {
		sfpar <- exp(parv)
		sfres <- evalHillModel_sf(conc,sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*sfact),mean(wt2*(sfact^2)),
				 mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[3]],pbounds[,4])[[1]]
		sfact <- start[[3]]+(ebnd-start[[3]])*sfact - act

		ovalue <- sum((weights*sfact)^2)
		derivs <- (2*(ebnd-start[[3]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)*sfpar
		return(list(value=ovalue,derivatives=derivs))

	}
	parv2fullpar <- function(parv) {
		sfpar <- exp(parv)
		sfact <- evalHillModel_sf(conc,sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*sfact),mean(wt2*(sfact^2)),
				 mean(wt2*act),mean(wt2*sfact*act))
		ebnd <- boundedOpt1d(mnv,start[[3]],pbounds[,4])[[1]]
		hpar <- c(sfpar,start[[3]],ebnd)
		return(hpar)

	}

	nls <- runBoundedOptim(valderivfunc,parv2fullpar,nstart,nbounds)
	nls$model <- c(1,2,4)
	nls$par <- nls$fullpar[nls$model]
	nls$mname <- "m3plc"
	return(nls)
}
fitHill4Par <- function(conc,act,weights,start,direction,pbounds) {
	ebounds <- pbounds[,3:4]
	obounds <- getHillOuterBounds(direction,ebounds)

	nbounds <- log(pbounds[,1:2])
	nstart <- log(start[1:2])

	wt2 <- (weights^2)/mean(weights^2)

	valderivfunc <- function(parv) {
		sfpar <- exp(parv)
		sfres <- evalHillModel_sf(conc,sfpar,calcderivs=TRUE)
		sfact <- sfres$value

		mnv <- c(mean(wt2*sfact),mean(wt2*(sfact^2)),
				 mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,obounds)
		sfact <- ebnds[[1]]+(ebnds[[2]]-ebnds[[1]])*sfact - act

		ovalue <- sum((weights*sfact)^2)
		derivs <- (2*(ebnds[[2]]-ebnds[[1]]))*as.vector(rbind(weights*weights*sfact)%*%sfres$derivatives)*sfpar
		return(list(value=ovalue,derivatives=derivs))
	}
	parv2fullpar <- function(parv) {
		sfpar <- exp(parv)
		sfact <- evalHillModel_sf(conc,sfpar,calcderivs=FALSE)

		mnv <- c(mean(wt2*sfact),mean(wt2*(sfact^2)),
				 mean(wt2*act),mean(wt2*sfact*act))
		ebnds <- boundedOpt2d(mnv,obounds)
		hpar <- c(sfpar,ebnds[1:2])
		return(hpar)
	}

	nls <- runBoundedOptim(valderivfunc,parv2fullpar,nstart,nbounds)
	nls$model <- c(1,2,3,4)
	nls$par <- nls$fullpar[nls$model]
	nls$mname <- "m4p"
	return(nls)
}

getHillOuterBounds <- function(direction,spbounds) {
	base_obounds <- emptyBound2d()
	if (is.finite(spbounds[1,1])) { base_obounds <- addBound2d(base_obounds,c(1,0,spbounds[1,1]),"E0min") }
	if (is.finite(spbounds[2,1])) { base_obounds <- addBound2d(base_obounds,c(-1,0,-spbounds[2,1]),"E0max") }
	if (is.finite(spbounds[1,2])) { base_obounds <- addBound2d(base_obounds,c(0,1,spbounds[1,2]),"Efmin") }
	if (is.finite(spbounds[2,2])) { base_obounds <- addBound2d(base_obounds,c(0,-1,-spbounds[2,2]),"Efmax") }
	if (direction!=0) { base_obounds <- addBound2d(base_obounds,c(-direction,direction,0),"Direction") }
	if (is.null(base_obounds)) { stop("Bounds specified by 'bounds' and 'direction' parameters cannot be satisified.") }
	return(base_obounds)
}

runBoundedOptim <- function(vdfunc,fpfunc,nstart,nbounds,parscale=NULL) {
	derivs <- NULL
	parfunc <- function(p) {
		p <- pmin(pmax(p,nbounds[1,]),nbounds[2,])
		vd <- vdfunc(p)
		derivs <<- vd$derivatives
		return(vd$value)
	}
	derivfunc <- function(p) { return(derivs) }

	if (is.null(parscale)) { control <- list(maxit=1000) }
	else { control <- list(maxit=1000,parscale=parscale) }

	nls <- stats::optim(nstart,parfunc,derivfunc,method="L-BFGS-B",lower=nbounds[1,],
				 upper=nbounds[2,],control=control)
	nls$par <- pmin(pmax(nls$par,nbounds[1,]),nbounds[2,])
	nls$fullpar <- fpfunc(nls$par)
	return(nls)
}
