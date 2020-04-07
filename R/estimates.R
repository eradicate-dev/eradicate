
#' calcN
#'
#' @description
#' \code{calcN} estimates abundance for a defined region from a fitted model.  The default is for the population defined by the sampling units.  The user can optionally supply a data.frame
#' of covariate values for any spatial unit of interest.
#'
#' @param obj A fitted model object.
#' @param newdata An (optional) \code{data.frame} of covariates for spatial units of interest.
#'    There must by covariate values for every parameter in \code{obj}.
#' @param off.set Either a scalar offset value to apply to each spatial unit
#' for prediction (e.g. cell area) or a vector of the same length as \code{nrow(newdata)}.
#'
#' @return a \code{data.frame} giving the predictions for each spatial unit
#'  as well as the overall abundance estimate for the region with associated
#'  SE and confidence intervals.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  mod <- nmix(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
calcN <- function(obj, ...){
  # method generic
  UseMethod("calcN", obj)
}

linearComb <- function(obj, ...){
  # method generic
  UseMethod("linearComb", obj)
}

backTransform <- function(obj, ...){
  # method generic
  UseMethod("backTransform", obj)
}

#' SE
#'
#' @description extracts standard errors for efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
SE <- function(obj, ...){
  # method generic
  UseMethod("SE", obj)
}

#' summary.efit
#'
#' \code{summary} summarises an efit object giving the estimated parameters for the
#'  model in \code{obj} (on the link scale) with associated SE and confidence
#'  intervals.
#'
#' @param obj A fitted model object.
#'
#' @return a \code{data.frame}
#'
#' @examples
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  mod <- nmix(~1, ~1, data=emf)
#'  summary(emf)
#'
#' @export
#'
summary.efit<- function(object, ...)
{
  type<- object$types
  out.list<- list()
  for(i in 1:length(type)) {
    ests <- object[[type[i]]]$estimates
    SEs <- SE(object, type[i])
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    invlink <- object[[type[i]]]$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit",
                   cloglog = "cloglog")
    cat(object[[type[i]]]$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, digits = 3)
    out.list[[type[i]]]<- outDF
    invisible(outDF)
    cat(" ","\n")
  }
  invisible(out.list)
}

#' summary.efitGP
#'
#' \code{summary} summarises an efitGP object giving the estimated parameters for N
#'  and catchability in \code{obj} (on the link scale) with associated SE and confidence
#'  intervals.
#'
#' @param obj A fitted model object.
#'
#' @return a \code{data.frame}
#'
#' @examples
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  mod <- nmix(~1, ~1, data=emf)
#'  summary(emf)
#'
#' @export
#'
summary.efitGP<- function(object, ...)
{
  type<- object$types
  out.list<- list()
  for(i in 1:length(type)) {
    ests <- object[[type[i]]]$estimates
    SEs <- SE(object, type[i])
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    invlink <- object[[type[i]]]$invlink
    link <- switch(invlink,
                   exp = "log",
                   cloglog = "cloglog",
                   logistic = "logit")
    cat(object[[type[i]]]$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, digits = 3)
    out.list[[type[i]]]<- outDF
    invisible(outDF)
    cat(" ","\n")
  }
  invisible(out.list)
}

# Compute linear combinations of estimates using coefficients

linearComb.efit<- function(obj, coefficients, off.set = NULL, ...)
{
  estimates<- obj$state$estimates
  covMat<- obj$state$covMat
  if(!is(coefficients, "matrix"))
    coefficients <- t(as.matrix(coefficients))
  if(ncol(coefficients) != length(estimates)) stop("error - wrong number of covariates")
  if (is.null(off.set))
    off.set <- rep(0, nrow(coefficients))
  e <- as.vector(coefficients %*% estimates) + off.set
  v <- coefficients %*% covMat %*% t(coefficients)
  umlc <- list(estimates = e, covMat = v, coefficients = coefficients,
               invlink = obj$state$invlink, invlinkGrad = obj$state$invlinkGrad)
  class(umlc)<- c("efit", class(umlc))
  umlc
}

backTransform.efit<- function(obj, ...) {
  # Delta method VAR
  invlink<- obj$invlink
  invlinkGrad<- obj$invlinkGrad
  estimate<- obj$estimates
  covMat<- obj$covMat
  e <- do.call(invlink,list(estimate))
  grad <- do.call(invlinkGrad,list(estimate))

  if(length(estimate) > 1) {
    v <- diag(grad) %*% covMat %*% diag(grad)
  } else {
    v <- grad^2 * covMat
  }
  umbt <- list(estimates = e, covMat = v)
  umbt
}

backTransform.efitGP<- function(obj, type="state", ...) {
  # Delta method VAR
  ests<- obj[[type]]
  invlink<- ests$invlink
  invlinkGrad<- ests$invlinkGrad
  estimate<- ests$estimates
  covMat<- ests$covMat
  e <- do.call(invlink,list(estimate))
  grad <- do.call(invlinkGrad,list(estimate))

  if(length(estimate) > 1) {
    v <- diag(grad) %*% covMat %*% diag(grad)
  } else {
    v <- grad^2 * covMat
  }
  umbt <- list(estimates = e, covMat = v)
  umbt
}

#' @rdname calcN
#' @export
calcN.efit<- function(obj, newdata, off.set=NULL, CI.level=0.95, ...) {
  if(missing(newdata) || is.null(newdata)) {
    origdata <- obj$data
    M <- numSites(origdata)
    if(is.null(siteCovs(origdata))) {
      newdata <- data.frame(Intercept = rep(1, M))
    } else {
      newdata <- siteCovs(origdata)
    }
  }
  design <- getDesign(obj, newdata)
  X<- design$X
  M<- nrow(X)
  if(!is.null(off.set) & length(off.set) == 1) off.set<- rep(off.set, M)
  lc<- linearComb(obj, coefficients=X, off.set=off.set)
  est<- backTransform(lc)
  V<- est$covMat
  Nhat<- sum(est$estimates)
  varN<- sum(est$covMat)
  seN<- sqrt(varN)
  cv<- sqrt(varN)/Nhat
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr<- Nhat*z
  upr<- Nhat/z

  bigN<- data.frame(N=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1))
  list(cellpreds=est$estimates, Nhat=bigN)
}

#' @rdname calcN
#' @export
calcN.efitM<- function(obj, newdata, off.set=NULL, CI.level=0.95, ...) {
  if(missing(newdata) || is.null(newdata)) {
    origdata <- obj$data
    M <- numSites(origdata)
    if(is.null(siteCovs(origdata))) {
      newdata <- data.frame(Intercept = rep(1, M))
    } else {
      newdata <- siteCovs(origdata)
    }
  }
  design <- getDesign(obj, newdata)
  X<- design$X
  M<- nrow(X)
  if(!is.null(off.set) & length(off.set) == 1) off.set<- rep(off.set, M)
  lc<- linearComb(obj, coefficients=X, off.set=off.set)
  est<- backTransform(lc)
  V<- est$covMat
  Pocc<- sum(est$estimates)/ M # mean occupancy
  varN<- sum(est$covMat) * (1/M^2) # delta method VAR
  seN<- sqrt(varN)
  cv<- seN/Pocc
  z <- qnorm((1-CI.level)/2, lower.tail = FALSE)
  lwr<- Pocc - seN*z
  upr<- Pocc + seN*z
  bigN<- data.frame(Pocc=round(Pocc,2),se=round(seN,3), lcl=round(lwr,2), ucl=round(upr,2))
  list(cellpreds=est$estimates, Occ=bigN)
}

#' @rdname calcN
#' @export
calcN.efitR<- function(obj, off.set=NULL, CI.level=0.95, ...) {
  # Currently only makes sense to predict to original data
  origdata <- obj$data
  M <- numSites(origdata)
  if(is.null(siteCovs(origdata))) {
    newdata <- data.frame(Intercept = rep(1, M))
  } else {
    newdata <- siteCovs(origdata)
  }
  design<- getDesign(obj, newdata)
  tot.rem<- sum(origdata$y, na.rm=TRUE)
  X<- design$X
  M<- nrow(X)
  if(!is.null(off.set) & length(off.set) == 1) off.set<- rep(off.set, M)
  lc<- linearComb(obj, coefficients=X, off.set=off.set)
  est<- backTransform(lc)
  V<- est$covMat
  Nhat<- sum(est$estimates)
  Nresid<- Nhat - tot.rem
  varN<- sum(est$covMat)
  seN<- sqrt(varN)
  cv<- sqrt(varN)/Nhat
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr<- Nhat*z
  upr<- Nhat/z
  lwr1<- Nresid*z
  upr1<- Nresid/z

  bigN<- data.frame(N=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1))
  littleN<- data.frame(N = round(Nresid,1),se=round(seN,1), lcl=round(lwr1,1), ucl=round(upr1,1))
  list(cellpreds=est$estimates, Nhat=bigN, Nresid=littleN)
}

#' @rdname calcN
#' @export
calcN.efitGP<- function(obj, CI.level=0.95, ...) {
  # Only need to calc sensible CI using the methods in
  # Chao (1989) Biometrics 45(2), 427-438
  x <- obj$data
  R<- sum(x$catch)
  ests<- backTransform(obj, type="state")
  N<- ests$estimates
  se.N<- sqrt(diag(ests$covMat))
  cv.N<- se.N/N
  z <- qt((1-CI.level)/2, length(x$catch) - 2, lower.tail = FALSE)
  asymp.N <- exp(z * sqrt(log(1 + ((se.N^2)/((N - R)^2)))))
  lcl.N <- R + (N - R)/asymp.N
  ucl.N <- R + (N - R)*asymp.N
  Nr<- N - R
  lcl.Nr<- (N - R)/asymp.N
  ucl.Nr<- (N - R)*asymp.N
  cv.Nr<- se.N/Nr
  bigN<- data.frame(N=round(N), se=round(se.N,1), lcl=round(lcl.N,1), ucl=round(ucl.N,1))
  littleN<- data.frame(N = round(Nr), se=round(se.N,1),lcl=round(lcl.Nr,1), ucl=round(ucl.Nr,1))
  list(Nhat=bigN, Nresid=littleN)
}

#' @rdname SE
#' @export
SE.efit<- function(obj, type, ...){
  v<- obj[[type]]$covMat
  sqrt(diag(v))
}


#' profileCI
#'
#' @description extracts profile likelihood confidence intervals for parameters
#' on the link scale.
#' @param obj A fitted model object.
#' @param type parameter type (i.e. "state","det")
#' @param level confidence level
#'
#' @export
#'
profileCI<- function(object, type, level = 0.95) {
  parm <- 1:length(object[[type]]$estimates)
  nllFun <- nllFun(object)
  ests <- mle(object)
  nP <- length(parm)
  ci <- matrix(NA, nP, 2)
  types <- object$types
  numbertable <- list()
  for(i in 1:length(types)) {
    length.est <- length(object[[types[i]]]$estimates)
    numbertable[[i]] <- data.frame(type = rep(types[i], length.est), num = seq(length.est))
  }
  numbertable<- do.call('rbind', numbertable)
  allparms<- which(numbertable$type == type & numbertable$num %in% parm)
  multiple<- c(2,4,8,12)
  for(i in 1:nP) {
    cat("Profiling parameter",i,"of",nP,"...")
    se <- SE(object, type)
    whichPar<- allparms[i]
    for(mult in multiple) {
      ci[i,] <- calc.profileCI(nllFun, whichPar=whichPar, MLE=ests,
                               interval=ests[whichPar] + mult*se[i]*c(-1,1), level=level)
      if(all(is.finite(ci[i,]))) break
    }
    cat(" done.\n")
  }
  colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
  if(any(!is.finite(ci)))
    warning("At least one endpoint of profile confidence interval is on the boundary.",
            call. = FALSE)
  return(ci)
}
