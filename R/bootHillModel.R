

#' Estimate Bootstrapped Confidence Intervals on Hill Model Parameters
#'
#' By bootsttrapping a large number of vectors, this function estimates
#' confidence intervals on the paramters of the given Hill model.  If the model
#' already has confidence intervals estimated, they will be replaced with a
#' warning
#'
#' @param hfit An object of class `hillrm`
#' @param ciLevs The lower and upper p-values for the estimated confidence
#' interval.  The default values, 0.025 and 0.975, produce a 95% confidence
#' interval.
#' @param numBoot The number of bootstrapped coefficients to be sampled.  If
#' `NULL` (the default) will be selected to have at least 10 samples lie
#' the selected interval, with a minimum of 100 and a maximum of 1000.
#'
#' @return An object of class `hillrm`, containing all the values found in any
#' `hillrm` object (see [fitHillModel()]) as well as the following three
#' values:
#' * `ciLevs`: the values used to set the bounds of the confidence intervals
#' * `ciCoefs`: a width-4 array of all bootstrapped Hill model coefficents
#' sampled by the function
#' * `ciMat`: a 2-by-4 array of values representing the estimated confidence
#' intervals on the four Hill model parameters
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#'
#' hfit <- fitHillModel(conc,response,c(1,2,3,4),start=c(0.5,1,0,100))
#' cihfit <- calcHillBootstrap(hfit)
calcHillBootstrap <- function(hfit,ciLevs=c(0.025,0.975),numBoot=NULL) {
	if (!inherits(hfit,"hillrm")) {
		stop("Object 'hfit' must be of class 'hillrm'.")
	}
	if (!is.null(hfit$ciLevs)) {
		warning("This Hill fit already has parameter confidence intervals estimated. These will be deleted and replaced.")
		hfit$ciLevs <- NULL
		hfit$ciCoefs <- NULL
		hfit$ciMat <- NULL
	}
	if (is.null(numBoot)) { numBoot <- round(max(min(10/(1-ciLevs[2]+ciLevs[1]),1000),100)) }

	bcoefs <- array(NA,dim=c(numBoot,4))
	for (i in 1:numBoot) {
		bact <- hfit$fitted.values+sample(hfit$residuals,length(hfit$residuals),replace=TRUE)
		tfit <- try(fitHillModel(hfit$conc,bact,hfit$model,hfit$coefficients,hfit$direction,hfit$pbounds[1,],hfit$pbounds[2,]),silent=TRUE)
		if (!inherits(tfit,"try-error")) { bcoefs[i,] <- tfit$coefficients }
	}
	bcoefs <- bcoefs[!is.na(bcoefs[,1]),]
	qmat <- t(apply(bcoefs,2,stats::quantile,probs=ciLevs))

	structure(
		c(unclass(hfit),list(ciLevs=ciLevs,ciCoefs=bcoefs,ciMat=qmat)),
		class="hillrm"
	)
}


#' Estimate an confidence interval on a Hill model property
#'
#' Given a function from Hill model parameters to one or more model properties,
#' this function produces a confidence interval on that value or those values
#' using the bootstrapped model coefficents produced by [calcHillBootstrap()].
#' This is useful for estimating confidence intervals on other values like
#' IC50, or generating confidence intervals on fitted values for plots.
#'
#' @param hfit An object of class `hillrm`, with the `ciCoefs` property
#' produced by [calcHillBootstrap()]
#' @param parfunc A function from a four parameter Hill model vector (see
#' [evalHillModel()]) to a single value or a vector of values with a consisten
#' length
#' @param civals An optional set of upper and lower bounds on the confidence
#' interval to be estimated.  If `NULL`, the default, the `ciLevs` property
#' from [calcHillBootstrap()] will be used.
#'
#' @return An n-by-3 array, where n is the length of the vector produced by
#' `parfunc`.  The first row is the lower bound of the confidence interval,
#' the second row is the function evaluated at the best-fit Hill model, and the
#' third row is the upper bound of the confidence interval.
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,75)
#' response <- evalHillModel(conc, hpar) + rnorm(length(conc),sd=7.5)
#'
#' hfit <- fitHillModel(conc,response,c(1,2,3,4),start=c(0.5,1,0,100))
#' cihfit <- calcHillBootstrap(hfit)
#'
#' ic50_ci <- calcHillConfInt(cihfit,function(h) invertHillModel(50,h))
calcHillConfInt <- function(hfit,parfunc,civals=NULL) {
	if (!inherits(hfit,"hillrm")) {
		stop("Object 'hfit' must be of class 'hillrm'.")
	}
	if (!inherits(parfunc,"function")) { stop("Input 'parfunc' must be a function of one variable.") }
	if (is.null(hfit$ciLevs)) { stop("Input 'hfit' must have bootstrapped coefficients.") }
	if (is.null(civals)) { civals <- hfit$ciLevs }

	outval <- parfunc(hfit$coefficients)
	outmat <- array(outval,dim=c(length(outval),1))
	outmat <- cbind(array(NA,dim=c(length(outval),nrow(hfit$ciCoefs))))
	for (b in 1:nrow(hfit$ciCoefs)) { outmat[,b] <- parfunc(hfit$ciCoefs[b,]) }

	outci <- apply(outmat,1,stats::quantile,probs=civals)
	fullout <- cbind(as.vector(outci[1,]),outval,as.vector(outci[2,]))
	return(fullout)
}
