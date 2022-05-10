#--------------------------------------------------------------------------------------------- Libraries

library(purrr) 
library(plyr) 
library(dplyr)
library(stringr)
library(nlme)
library(matrixStats)
library(ggplot2)
library(RConics)
library(zoo)
#--------------------------------------------------------------------------------------------- Functions

#' Negation of %in% function
#' 
`%notin%` <- Negate(`%in%`)

#--------------------------------------------------------------------------------------------- 

#' Creates data.frame of Midpoints vs Rates for
#'  individual subjects curves
#'  
#' @param BuildDictionary (list) A list of arguments from  \emph{FitDiseaseProgressionCurve} function
#' @return (data.frame) A data frame with columns ID, Midpoints, Rates
#' 

EstimateMeanSlope <- function(BuildDictionary) {
  method.logical  <- BuildDictionary[["individual.lm"]]
  data            <- BuildDictionary[["data"]]
  formula.fixed   <- BuildDictionary[["formula.fixed"]]
  formula.random  <- BuildDictionary[["formula.random"]]
  model.control   <- BuildDictionary[["lmeControl"]]
  rate.vec        <- c()
  midpoint.vec.av <- midpoint.vec.rate <- c()
  
  if(method.logical) { # this part of the function makes a linear regression
    split.data      <- split(data, data$ID)
    ids             <- names(split.data)
    for(i in 1:length(split.data)) {
      subj          <- split.data[[i]]
      subj.lm       <- lm(formula = as.formula(formula.fixed), data = subj)
      int           <- coef(subj.lm)[["(Intercept)"]]
      rate          <- coef(subj.lm)[["Time_Since_Baseline"]]
      pred.model    <- predict(subj.lm)
      rate.vec      <- append(rate.vec, rate)
      midpoint.rate <- int + (rate * (.5 * max(subj$Time_Since_Baseline)))
      midpoint.vec.rate  <- append(midpoint.vec.rate, midpoint.rate)
    }
    mean.slope.data <- data.frame("ID"        = ids,
                                  "Rates"     = rate.vec,
                                  "Midpoints" = midpoint.vec.rate)
    
  } else { # this part of the function creates a mixed effects model 
    model      <- lme(as.formula(formula.fixed), 
                      random  = as.formula(formula.random), 
                      data    = data,
                      control = model.control)
    re         <- nlme :: ranef(model)
    fe         <- nlme :: fixef(model)
    splitdata  <- split(data, data$ID)
    ids        <- names(splitdata)
    for(i in ids) {
      subj       <- splitdata[[i]]
      int        <- re["(Intercept)"][i,]
      int        <- int + fe[["(Intercept)"]]
      rate       <- re["Time_Since_Baseline"][i, ]
      rate       <- rate + fe[["Time_Since_Baseline"]]
      rate.vec   <- append(rate.vec, rate)
      pred.val   <- stats::predict(object = model)
      keeps      <- which(names(pred.val) == i)
      pred.val   <- pred.val[keeps]
      
      midpoint.rate <- int + (rate * (.5 * max(subj$Time_Since_Baseline)))
      midpoint.vec.rate  <- append(midpoint.vec.rate, midpoint.rate)
      
    }
    mean.slope.data <- data.frame("ID"          = ids,
                                  "Rates"       = rate.vec,
                                  "Midpoints"   = midpoint.vec.rate)
    
  }
  
  return(mean.slope.data)
}
#---------------------------------------------------------------------------------------------

#' Fits 3rd degree poly over Midpoints/Rates data
#' 
#' @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#' @return (list) Coefficients of fitted 3rd degree polynomial
#' 
FitPolynomial <- function(EstimateMeanSlopeOutput) {
  data <- EstimateMeanSlopeOutput
  third.degree.poly <- lm(Rates ~ poly(Midpoints, 3, raw = TRUE), 
                          data = data)
  poly.coefs <- coef(third.degree.poly)
  poly.coefs <- list("coef.x3"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)3"]], 
                     "coef.x2"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)2"]], 
                     "coef.x1"  = poly.coefs[["poly(Midpoints, 3, raw = TRUE)1"]], 
                     "int"      = poly.coefs[["(Intercept)"]])
  
  return(poly.coefs)
}

#---------------------------------------------------------------------------------------------
#' Calculates the roots (real) of the fitted polynomial
#' 
#' @param FitPolynomialOutput (list) Output from \emph{FitPolynomial} function
#' @return (vector) Vector of real roots of fitted polynomial
#' 
FindRealRoots <- function(FitPolynomialOutput) {
  poly.coefs  <- FitPolynomialOutput
  poly.coefs  <- c(poly.coefs[["coef.x3"]], poly.coefs[["coef.x2"]], 
                   poly.coefs[["coef.x1"]], poly.coefs[["int"]])
  roots        <- RConics::cubic(poly.coefs)
  roots        <- as.complex(roots)
  roots.keep   <- which(Im(roots) == 0)
  roots        <- roots[roots.keep]
  roots        <- Re(roots)
  return(roots)
}

#---------------------------------------------------------------------------------------------
#' Checks to see whether the roots satisfy criterion for integration
#'  
#' @param FindRealRootsOutput (vector) Output from \emph{FindRealRoots} function
#' @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#' @return (list) List of all real roots with logical value indicating
#'  whether they satisfy integration conditions
#'  
CheckRealRoots <- function(FindRealRootsOutput, EstimateMeanSlopeOutput) {
  all.roots     <- c()
  roots.satisfy <- list()
  roots         <- FindRealRootsOutput
  midpoints     <- EstimateMeanSlopeOutput$Midpoints
  min.midpoint  <- min(midpoints)
  max.midpoint  <- max(midpoints)
  if(length(roots) == 0) {
    roots.satisfy[[1]] <- list("root" = NA,
                               "satisfies_condition" = TRUE)
    all.roots <- TRUE
  } else {
    for(i in 1:length(roots)) {
      if(roots[i] >= min.midpoint & roots[i] <= max.midpoint) {
        logi               <- FALSE
        all.roots          <- append(all.roots, logi)
        roots.satisfy[[i]] <- list("root" = roots[i],
                                   "meets_condition" = logi)
      } else {
        logi               <- TRUE
        all.roots          <- append(all.roots, logi)
        roots.satisfy[[i]] <- list("root" = roots[i],
                                   "meets_condition" = logi)
      }
    }
  }
  all.roots              <- all(all.roots)
  roots.satisfy[["all_roots_satisfy?"]] <- all.roots
  return(roots.satisfy)
}

#---------------------------------------------------------------------------------------------
#' Defines functions for 3rd degree polynomial and reciprocal of
#'  3rd degree polynomial from fitted polynomial coefficients
#'  
#'  @param FitPolynomialOutput (list) Output from \emph{FitPolynomial} function
#'  @return (list) List of functions for checking polynomial fit visually 
#'   and integrating polynomial reciprocal
#'   
DefinePolynomialCurveAndReciprocal <- function(FitPolynomialOutput) {
  poly.coefs      <- FitPolynomialOutput
  PolynomialCurve <- function(x) {
    (poly.coefs[["coef.x3"]] * (x^3)) +
      (poly.coefs[["coef.x2"]] * (x^2)) +
      (poly.coefs[["coef.x1"]] * (x))   +
      poly.coefs[["int"]]
  }
  ReciprocalCurve <- function(x) {
    1 / ((poly.coefs[["coef.x3"]] * (x^3)) +
           (poly.coefs[["coef.x2"]] * (x^2)) +
           (poly.coefs[["coef.x1"]] * (x))   +
           poly.coefs[["int"]])
  }
  
  return(list("Polynomial_Function" = PolynomialCurve,
              "Reciprocal_Function" = ReciprocalCurve))
}

#---------------------------------------------------------------------------------------------
#' Calculates the window of integration based on min and max midpoints,
#'  roots of the fitted polynomial, and the direction of the curve.
#'  
#'  @param CheckRealRootsOutput (list) Output from \emph{CheckRealRoots} function
#'  @param EstimateMeanSlopeOutput (data.frame) Output from \emph{EstimateMeanSlope} function
#'  @param DefinePolynomialCurveAndReciprocalOutput (list) Output from \emph{DefinePolynomialCurveAndReciprocal} function
#'  @param seq.by (numeric) # of times to integrate along curve (default is 1000)
#'  @return (list) a list of values containing the integration domain,
#'   whether the domain was subsetted due to roots of the polynomial, and
#'   the direction of the curve (positive/negative)
#'   
CalculateBoundsofIntegration <- function(CheckRealRootsOutput, 
                                         EstimateMeanSlopeOutput, 
                                         DefinePolynomialCurveAndReciprocalOutput) {
  
  polynomial.curve <- DefinePolynomialCurveAndReciprocalOutput[["Polynomial_Function"]]
  min.midpoint.row <- which.min(EstimateMeanSlopeOutput$Midpoints)
  max.midpoint.row <- which.max(EstimateMeanSlopeOutput$Midpoints)
  min.midpoint     <- EstimateMeanSlopeOutput["Midpoints"][min.midpoint.row, ]
  max.midpoint     <- EstimateMeanSlopeOutput["Midpoints"][max.midpoint.row, ]
  min.midpoint.y   <- polynomial.curve(min.midpoint)
  max.midpoint.y   <- polynomial.curve(max.midpoint)
  number.positive.rates <- length(which(EstimateMeanSlopeOutput$Rates > 0))
  number.negative.rates <- length(which(EstimateMeanSlopeOutput$Rates < 0))
  a <- min.midpoint
  b <- max.midpoint
  head_subset <- FALSE
  tail_subset <- FALSE
  if(number.positive.rates < number.negative.rates) {
    direction <- "decreasing"
  } else {
    direction <- "increasing"
  }
  roots <- CheckRealRootsOutput
  roots[["all_roots_satisfy?"]] <- NULL
  roots.vals         <-  unlist(map(roots, pluck, 1))
  integration.points <- c("min"=min.midpoint, "roots"=roots.vals, "max"= max.midpoint)
  if(direction=="increasing") {
    integration.points <- integration.points[order(integration.points, decreasing = FALSE)]
  } else {
    integration.points <- integration.points[order(integration.points, decreasing =TRUE)]
  }
  means.points <- rollapply(integration.points, 
                            width = 2, by = 1, 
                            FUN = mean, align = "left")
  direc <- polynomial.curve(means.points)
  direc <- unlist(Map(function(x) {if(x > 0) {"increasing"} else {"decreasing"}}, direc))
  direc <- c(NA, direc)
  a <- which(names(integration.points) == "min")
  b <- which(names(integration.points) == "max")
  keep.points <- data.frame("integration_points" = integration.points,
                            "direction"          =  direc)
  keep.points <- keep.points[a:b,]
  keep.points["direction"][1,] <- NA
  if(direction=="increasing") {
    if(nrow(keep.points) == 2) {
      keep.points$integration_start <- keep.points["integration_points"][1,]
      keep.points$integration_end   <- keep.points["integration_points"][2,]
    } else if(nrow(keep.points)==3) {
      integration.end.row <- which(keep.points$direction==direction)
      keep.points$integration_start <- keep.points["integration_points"][integration.end.row - 1,]
      keep.points$integration_end <- keep.points["integration_points"][integration.end.row,]
    } else if(nrow(keep.points) > 3) {
      integration.end.row <- which(keep.points$direction==direction)
      integration.end.row <- min(integration.end.row)
      keep.points$integration_start <- keep.points["integration_points"][integration.end.row - 1,]
      keep.points$integration_end <- keep.points["integration_points"][integration.end.row,]
    }
  } else {
    if(nrow(keep.points) == 2) {
      keep.points$integration_start <- keep.points["integration_points"][1,]
      keep.points$integration_end   <- keep.points["integration_points"][2,]
    } else if(nrow(keep.points)==3) {
      integration.end.row <- which(keep.points$direction==direction)
      keep.points$integration_start <- keep.points["integration_points"][integration.end.row,]
      keep.points$integration_end <- keep.points["integration_points"][integration.end.row + 1,]
    } else if(nrow(keep.points) > 3) {
      integration.end.row <- which(keep.points$direction==direction)
      integration.end.row <- min(integration.end.row)
      keep.points$integration_start <- keep.points["integration_points"][integration.end.row,]
      keep.points$integration_end <- keep.points["integration_points"][integration.end.row + 1,]
    }
  }
  return(list("roots.frame" = keep.points, "direction" = direction))
}


#---------------------------------------------------------------------------------------------
#' Calculates integral of polynomial along polynomial domain
#' 
#' @param DefinePolynomialCurveAndReciprocalOutput (list) Output of \emph{DefinePolynomialCurveAndReciprocal} function
#' @param CalculateBoundsofIntegrationOutput (list) Output of \emph{CalculateBoundsofIntegration} function
#' @return (vector) Vector of integrated values 
#' 
IntegratePolynomial <- function(integration.domain, polyfunction) {
  integrated.values.vector  <- c()
  for(k in 1:length(integration.domain)) {
    integration.val <- integrate(polyfunction,
                                 lower = integration.domain[1],
                                 upper = integration.domain[k])$val
    integration.val <- suppressWarnings(as.numeric(integration.val))
    integrated.values.vector <- suppressWarnings(append(integrated.values.vector, integration.val))
    if(is.na(integration.val)) {
      warning("NA value generated during integration")
    }
  }
  return(integrated.values.vector)
}

#---------------------------------------------------------------------------------------------
#' Calculates mean and standard error at each point along the domain of the
#'  curve from boostrap iterations
#'  
#' @param CalculateBoundsofIntegrationOutput (list) Output of \emph{CalculateBoundsofIntegration} function
#' @param BootStrappedDF (data.frame) Data frame of bootstrapped values
#' @return (data.frame) Disease progression model data
#' 

CalculateSE <- function(integration.domain, BootStrappedDF, n_iter) {
  BootStrappedDF <- as.matrix(BootStrappedDF)
  response       <- integration.domain
  Conf.Low       <- rep(NA, nrow(BootStrappedDF))
  Conf.Hi        <- rep(NA, nrow(BootStrappedDF))
  domain         <- rep(NA, nrow(BootStrappedDF))
  n              <- n_iter
  for(i in 1 : nrow(BootStrappedDF)) {
    row.mean      <- mean(BootStrappedDF[i,], na.rm = TRUE)
    variance      <- sd(BootStrappedDF[i,], na.rm = TRUE)
    Conf.Low[i]   <- quantile(BootStrappedDF[i,], .05, na.rm = TRUE)
    Conf.Hi[i]    <- quantile(BootStrappedDF[i,], .95, na.rm = TRUE)
    domain[i]     <- row.mean
    
  }
  disease.progression.data <- cbind(response, domain, Conf.Low, Conf.Hi)
  colnames(disease.progression.data) <- c("Response", "Domain", "CI_Low", "CI_Hi")
  return(disease.progression.data)
}

#---------------------------------------------------------------------------------------------
#' Determines whether to reorder the disease progression model data depending on if the curve
#'  is decreasing
#'  
#' @param CalculateSEOutput (matrix) Output from \emph{CalculateSE} function
#' @param CalculateBoundsofIntegrationOutput (list) Output from \emph{CalculateBoundsofIntegration} function
#' @return (data.frame) Disease progression data reordered if necessary
#' 
ReorderIfDecreasing <- function(CalculateSEOutput, direction) {
  disease.progression.data <- as.data.frame(CalculateSEOutput)
  if(direction == "decreasing") {
    disease.progression.data["Domain"]              <- disease.progression.data["Domain"] * -1
    disease.progression.data["CI_Low"]              <- disease.progression.data["CI_Low"] * -1
    disease.progression.data["CI_Hi"]               <- disease.progression.data["CI_Hi"]  * -1
    disease.progression.data["Response"]            <- disease.progression.data["Response"][nrow(disease.progression.data):1,]
  }
  return(disease.progression.data)  
}

#--------------------------------------------------------------------------------------------- Generate plots
#' Plots disease progression curve
#' 
#' @param data (data.frame) Disease progression data frame
#' @return (ggplot) Plot of disease progression curve
#' 
PlotCurve <- function(data) {
  na.remove       <- which(is.na(data$Domain))
  if(length(na.remove) >= 1) {
    warning("NA values were removed generating plot")
    data            <- data[-na.remove,]
  }
  init.plot       <- ggplot(data = data, aes(x = Domain,  y = Response)) + geom_point()
  plot.curve      <- init.plot  + geom_line(aes(x=CI_Low, y = Response), linetype="dashed")
  plot.curve      <- plot.curve + geom_line(aes(x=CI_Hi, y = Response), linetype="dashed")
  return(plot.curve)
}

#---------------------------------------------------------------------------------------------
#' Plots Midpoints vs Rates
#' 
#' @param data (data.frame) Midpoints vs Rates data frame (EstimateMeanSlopeOutput)
#' @return (ggplot) Plot of Midpoints vs Rates
#' 
PlotMeanSlope <- function(data) {
  plot        <- ggplot(data = data, aes(x = Midpoints, y = Rates)) + geom_point()
  return(plot)
}

#---------------------------------------------------------------------------------------------
