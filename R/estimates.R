
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


#' coef
#'
#' @description extracts coefficients for efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
coef <- function(obj, ...){
  # method generic
  UseMethod("coef", obj)
}

#' vcmat
#'
#' @description extracts variance covaraiance matrix for efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
vcmat <- function(obj, ...){
  # method generic
  UseMethod("vcmat", obj)
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

#' calcP
#'
#' @description calculates detection probability from a fitted model.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
calcP <- function(obj, ...){
  # method generic
  UseMethod("calcP", obj)
}

#' occTraject
#'
#' @description extracts occupancy trajectories from occuMS objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
occTraject <- function(obj, ...){
  # method generic
  UseMethod("occTraject", obj)
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
  type<- names(object$estimates)
  out.list<- list()
  for(i in 1:length(type)) {
    ests <- object$estimates[[type[i]]]$estimates
    SEs <- SE(object, type[i])
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    invlink <- object$estimates[[type[i]]]$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit",
                   cloglog = "cloglog")
    cat(object$estimates[[type[i]]]$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, digits = 3)
    out.list[[type[i]]]<- outDF
    invisible(outDF)
    cat(" ","\n")
  }
  invisible(out.list)
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
  sites<- design$retained.sites
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj, "state")
  covMat<- obj$estimates$state$covMat
  if(ncol(X) != length(estimates)) stop("error - wrong number of covariates")
  if (is.null(off.set)) off.set <- rep(1, M)
  else if(length(off.set) == 1) off.set<- rep(off.set, M)
  # cellwise estimates
  eta <- as.vector(X %*% estimates)
  vc <- rowSums((X %*% covMat) * X) # equivalent to diag(X %*% covMat %*% t(X))
  ests<- do.call(invlink,list(eta))
  grad <- do.call(invlinkGrad,list(eta))
  v<- (grad * off.set)^2 * vc
  cellpreds<- data.frame(N = ests * off.set, se = sqrt(v), site = sites)
  # overall estimate
  Nhat<- off.set %*% ests
  gradN<- off.set %*% (grad * X)
  varN<- gradN %*% covMat %*% t(gradN)
  seN<- sqrt(varN)
  cv<- sqrt(varN)/Nhat
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr<- Nhat*z
  upr<- Nhat/z

  bigN<- data.frame(N=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1))
  list(cellpreds=cellpreds, Nhat=bigN)
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
  sites<- design$retained.sites
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj,"state")
  covMat<- obj$estimates$state$covMat
  if(ncol(X) != length(estimates)) stop("error - wrong number of covariates")
  if (is.null(off.set)) off.set <- rep(1, M)
  else if(length(off.set) == 1) off.set<- rep(off.set, M)
  # cellwise estimates
  eta <- as.vector(X %*% estimates)
  vc <- rowSums((X %*% covMat) * X) # equivalent to diag(X %*% covMat %*% t(X))
  ests<- do.call(invlink,list(eta))
  grad <- do.call(invlinkGrad,list(eta))
  v<- (grad * off.set)^2 * vc
  cellpreds<- data.frame(Pocc = ests * off.set, se = sqrt(v), site = sites)
  # Overall estimate
  Pocc<- off.set %*% ests / M # mean occupancy
  gradN<- off.set %*% (grad * X) / M
  varN<- gradN %*% covMat %*% t(gradN)
  seN<- sqrt(varN)
  cv<- seN/Pocc
  z <- qnorm((1-CI.level)/2, lower.tail = FALSE)
  lwr<- Pocc - seN*z
  upr<- Pocc + seN*z
  bigN<- data.frame(Pocc=round(Pocc,2),se=round(seN,3), lcl=round(lwr,2), ucl=round(upr,2))
  list(cellpreds=cellpreds, Occ=bigN)
}

#' @rdname calcN
#' @export
calcN.efitR<- function(obj, newdata, off.set=NULL, CI.level=0.95, ...) {
  if(missing(newdata) || is.null(newdata)) {
    origdata <- obj$data
    M <- numSites(origdata)
    if(is.null(siteCovs(origdata))) {
      newdata <- data.frame(Intercept = rep(1, M))
    } else {
      newdata <- siteCovs(origdata)
    }
  }
  tot.rem<- sum(obj$data$y, na.rm=TRUE)
  design <- getDesign(obj, newdata)
  X<- design$X
  M<- nrow(X)
  sites<- design$retained.sites
  mixture<- obj$mixture
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj, "state")
  covMat<- vcmat(obj, "state")
  if(mixture == "ZIP") {
    logit.psi<- coef(obj, "zeroinfl")
    logit.psi.var<- vcmat(obj, "zeroinfl")
    psi<- do.call(obj$estimates$zeroinfl$invlink, list(logit.psi))
    psi.grad<- 1/(exp(-logit.psi) + 1) # derivative of log(psi)
    psi.var<- psi.grad^2 * logit.psi.var
  }
  if(ncol(X) != length(estimates)) stop("error - wrong number of covariates")
  if (is.null(off.set)) off.set <- rep(1, M)
  else if(length(off.set) == 1) off.set<- rep(off.set, M)
  # cellwise estimates
  if(mixture == "ZIP") {
    eta <- as.vector(X %*% estimates) + log(1-psi)
    vc <- rowSums((X %*% covMat) * X) + as.vector(rep(psi.var, M))
    ests<- do.call(invlink,list(eta))
    grad <- do.call(invlinkGrad,list(eta))
    v<- (grad * off.set)^2 * vc
    cellpreds<- data.frame(N = ests * off.set, se = sqrt(v), site = sites)
  }
  else{
    eta <- as.vector(X %*% estimates)
    vc <- rowSums((X %*% covMat) * X) # equivalent to diag(X %*% covMat %*% t(X))
    ests<- do.call(invlink,list(eta))
    grad <- do.call(invlinkGrad,list(eta))
    v<- (grad * off.set)^2 * vc
    cellpreds<- data.frame(N = ests * off.set, se = sqrt(v), site = sites)
  }
  # overall estimate
  Nhat<- off.set %*% ests
  gradN<- off.set %*% (grad * X)
  varN<- gradN %*% covMat %*% t(gradN)
  Nresid<- Nhat - tot.rem
  seN<- sqrt(varN)
  cv<- sqrt(varN)/Nhat
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr<- Nhat*z
  upr<- Nhat/z
  lwr1<- Nresid*z
  upr1<- Nresid/z

  bigN<- data.frame(N=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1))
  littleN<- data.frame(N = round(Nresid,1),se=round(seN,1), lcl=round(lwr1,1), ucl=round(upr1,1))
  list(cellpreds=cellpreds, Nhat=bigN, Nresid=littleN)
}

#' @rdname calcN
#' @export
calcN.efitGP<- function(obj, CI.level=0.95, ...) {
  # Only need to calc sensible CI using the methods in
  # Chao (1989) Biometrics 45(2), 427-438
  x <- obj$data
  R<- sum(x$catch)
  invlink = obj$state$invlink
  eta<- obj$estimates$state$estimates
  covMat<- obj$estimates$state$covMat
  N<- do.call(invlink, list(eta))
  se.N<- sqrt(diag(covMat))
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


#' @rdname coef
#' @export
coef.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  obj$estimates[[type]]$estimates
}

#' @rdname vcmat
#' @export
vcmat.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  obj$estimates[[type]]$covMat
}

#' @rdname SE
#' @export
SE.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  v<- obj$estimates[[type]]$covMat
  sqrt(diag(v))
}

#' @rdname calcP
#' @export
calcP.efitR<- function(obj, na.rm = TRUE) {
  detformula <- as.formula(obj$detformula)
  lamformula <- as.formula(obj$lamformula)
  emf <- obj$data
  designMats <- getDesign.eFrame(emf, lamformula, detformula, na.rm = na.rm)
  y <- designMats$y
  V <- designMats$V
  V.offset <- designMats$V.offset
  if (is.null(V.offset))
    V.offset <- rep(0, nrow(V))
  M <- nrow(y)
  J <- ncol(y)
  pars <- coef(obj, type = "det")
  p <- plogis(V %*% pars + V.offset)
  p <- matrix(p, M, J, byrow = TRUE)
  pi <- do.call(removalPiFun, list(p = p))
  return(pi)
}

#' @rdname occTraject
#' @export
occTraject.efitMS<- function(obj, type = c("projected","smoothed"), mean=TRUE, ...){
  type <- match.arg(type, c("projected","smoothed"))
  if(identical(type,"projected")) {
    if(mean) obj$projected.mean
    else obj$projected
  }
  else if(identical(type,"smoothed")) {
    if(mean) obj$smoothed.mean
    else obj$smoothed
  }
  else stop("type must be one of 'smoothed' or 'projected' ")
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
  parm <- 1:length(object$estimates[[type]]$estimates)
  nllFun <- nllFun(object)
  ests <- mle(object)
  nP <- length(parm)
  ci <- matrix(NA, nP, 2)
  types <- object$types
  numbertable <- list()
  for(i in 1:length(types)) {
    length.est <- length(object$estimates[[types[i]]]$estimates)
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
