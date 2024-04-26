
#' Evaluate a Hill dose response model
#'
#' @param conc A vector of concentrations (including 0 or Inf)
#' @param hpar A four parameter vector specifying the Hill model. The values of
#' the parameter vector are, in order, the dose of median effect (also often
#' referred to as the EC50), the Hill slope, the minimal effect (observed when
#' no drug or dose is present), and the maximal effect (theoretically observed
#' when the drug or dose is infinite).
#'
#' @return A vector of response values the same length as `conc`
#' @export
#'
#' @examples
#' conc <- c(0,2^(-6:3),Inf)
#' hpar <- c(1,3,0,100)
#'
#' response <- evalHillModel(conc, hpar)
evalHillModel <- function(conc,hpar) {
	IDM <- hpar[[1]]
	n <- hpar[[2]]
	E0 <- hpar[[3]]
	Ef <- hpar[[4]]

	return(E0 + (Ef-E0)*evalHillModel_sf(conc,c(IDM,n),calcderivs=FALSE))
}


#' Calculates the required doses of a Hill dose response model
#'
#' This funciton takes one or more desired response values and determines the
#' doses that will produce the desired effects given a particular Hill dose
#' response model. This is useful for estimating things like IC50.
#'
#' @param effect A vector of desired response values
#' @param hpar A four parameter vector specifying the Hill model. Parameter
#' details are found in the documentation for [evalHillModel()]
#' @param invalidNA Specifies what to do with values that are outside the range
#' of the given Hill model.  If `FALSE` (the default), values "below" the given
#' range will be set to zero, and values "above" the given range will be set to
#' Inf.  If `TRUE`, invalid values will be set to `NA`.
#'
#' @return A vector of concentrations the same length as `effect`.
#' @export
#'
#' @examples
#' invertHillModel(0.5, c(1,2,0,0.7))
#'
#' invertHillModel(seq(0.1,0.9,by=0.1), c(0.1,4,0,0.65), invalidNA=TRUE)
invertHillModel <- function(effect,hpar,invalidNA=FALSE) {
	IDM <- hpar[[1]]
	n <- hpar[[2]]
	E0 <- hpar[[3]]
	Ef <- hpar[[4]]

	if (invalidNA) {
		lowValue <- NA
		highValue <- NA
	} else {
		lowValue <- 0
		highValue <- Inf
	}

	sc_effect <- (effect-E0)/(Ef-E0)
	fr_effect <- sc_effect/(1-sc_effect)
	conc <- ifelse(
		sc_effect < 0,
		lowValue,
		ifelse(
			sc_effect > 1,
			highValue,
			IDM*(fr_effect^(1/n))
		)
	)

	return(conc)
}


#' Calculate the area under the curve for a Hill dose response model
#'
#' The area under the curve, or AUC, is a commonly used and robust metric
#' for evaluating and comparing dose response models. The area is calculated in
#' a log-concentration space, and so is dependent not only on the concentration
#' bounds, but also on the base of the logarithm used.  By default, this
#' function follows the common convention of using base 10.
#'
#' @param hpar A four parameter vector specifying the Hill model. Parameter
#' details are found in the documentation for [evalHillModel()]
#' @param range A two element vector specifying the lower and upper bound of
#' area being calculated
#' @param baseline The reference baseline response around which the area is
#' being calculated. The default value of 0 is generally the most intuitive
#' choice, but a value of 1 (or 100 in percent) could be used if the
#' dose-response model is fitting relative survival.
#' @param logbase The base of the logarithm applied to concentrations
#'
#' @return A single value specifying the area under the curve in the given range
#' @export
#'
#' @examples
#' auc <- calculateHillAUC(c(1,3,0,75), c(0.05,10))
#'
#' # Area *over* the curve in survival studies
#' aoc <- -calculateHillAUC(c(0.1,2,1,0.1), c(0.01,1), baseline=1)
calculateHillAUC <- function(hpar,range,baseline=0,logbase=10) {
	IDM <- hpar[[1]]
	n <- hpar[[2]]
	E0 <- hpar[[3]]
	Ef <- hpar[[4]]

	((E0-baseline)*diff(log(range),base=logbase) +
	 	(Ef-E0)*log(
	 		(range[[2]]^n + IDM^n)/(range[[1]]^n + IDM^n),
	 		base=logbase)/n
	)
}

evalHillModel_sf <- function(conc,sfpar,calcderivs=FALSE) {
	IDM <- clip_positive(sfpar[[1]])
	n <- clip_positive(sfpar[[2]])

	pow <- (conc/IDM)^n   # [0, Inf]
	R0 <- 1/(1+pow)       # [0, 1]
	Rf <- 1-R0            # [0, 1]

	if(!calcderivs) { return(Rf) }
	padj <- (Rf^2)/pmax(pow,.Machine$double.xmin)

	dRfdIDM <- -padj*clip_positive(n/IDM)
	dRfdn <- padj*clip_finite(log(conc/IDM))

	return(list(value=Rf,derivatives=cbind(dRfdIDM,dRfdn)))
}

clip_lo <- function(n) pmax(n,-.Machine$double.xmax)
clip_hi <- function(n) pmin(n,.Machine$double.xmax)
clip_finite <- function(n) pmin(pmax(n,-.Machine$double.xmax),.Machine$double.xmax)
clip_positive <- function(n) pmin(pmax(n,.Machine$double.xmin),.Machine$double.xmax)
