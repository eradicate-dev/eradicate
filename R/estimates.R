
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
coef.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  obj$estimates[[type]]$estimates
}

#' AIC
#'
#' @description extracts the Akaike Information Criterion from efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
AIC.efit <- function(obj){
  # method generic
  return(obj$AIC)
}

#' vcov
#'
#' @description extracts variance covariance matrix from efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
vcov.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  obj$estimates[[type]]$covMat
}

#' se
#'
#' @description extracts standard errors for efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
se <- function(obj, ...){
  # method generic
  UseMethod("se", obj)
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

#' fitted
#'
#' @description calculates the fitted values for a efit model.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
fitted.efitR<- function(obj, na.rm = TRUE) {
  detformula <- as.formula(obj$detformula)
  lamformula <- as.formula(obj$lamformula)
  emf <- obj$data
  designMats <- getDesign.eFrame(emf, lamformula, detformula, na.rm = na.rm)
  X <- designMats$X
  X.offset <- designMats$X.offset
  if (is.null(X.offset)) {
    X.offset <- rep(0, nrow(X))
  }
  mixture<- obj$mixture
  invlink = obj$estimates$state$invlink
  estimates<- coef(obj, "state")
  eta <- as.vector(X %*% estimates + X.offset)
  ests<- do.call(invlink,list(eta))
  p <- calcP(obj, na.rm = na.rm) # P(detection | presence)
  fitted <- ests * p  # true for models with E[Y] = p * E[X]
  fitted
}


#' residuals
#'
#' @description calculates working residuals from an efit model object.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
residuals.efitR<- function(obj) {
  y <- getY(obj$data)
  e <-  fitted(obj, na.rm=FALSE)
  r <- y - e
  return(r)
}

#' traject
#'
#' @description extracts trajectories for open population models
#' \code{occuMS} and \code{remMNO} objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
traject <- function(obj, ...){
  # method generic
  UseMethod("traject", obj)
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
    SEs <- se(object, type[i])
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

#' summary.efitGPlist
#'
#' \code{summary} summarises an efitGPlist object giving the estimated parameters for the
#'  model in \code{obj} (on the link scale) with associated SE and confidence
#'  intervals.
#'
#' @param object A list of efitGP objects.
#'
#' @return a \code{list}
#'
#' @examples
#'  emf <- eFrameGP(removals, effort, session=session)
#'  mod <- remGP(emf)
#'  summary(emf)
#'
#' @export
#'
summary.efitGPlist<- function(object, ...)
{
  out.list<- lapply(object, summary.efit)
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
  row.names(bigN)<- "Total"
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
  row.names(bigN)<- "Total"
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
  M<- nrow(X) # for stacked data
  sites<- design$retained.sites
  mixture<- obj$mixture
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj, "state")
  covMat<- vcov(obj, "state")
  if(mixture == "ZIP") {
    logit.psi<- coef(obj, "zeroinfl")
    logit.psi.var<- vcov(obj, "zeroinfl")
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
  row.names(bigN)<- "Total"
  row.names(littleN)<- "Residual"
  list(cellpreds=cellpreds, Nhat=bigN, Nresid=littleN)
}

#' @rdname calcN
#' @export
calcN.efitGP<- function(obj, CI.level=0.95, CI.calc = c("norm","lnorm","boot"), nboot=500, ...) {
  # CI.calc = "LN" calculated using the log normal method in
  # Chao (1989) Biometrics 45(2), 427-438
  CI.calc <- match.arg(CI.calc)
  x <- obj$data
  R<- sum(x$catch)
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  eta<- obj$estimates$state$estimates
  grad <- do.call(invlinkGrad,list(eta))
  covMat<- obj$estimates$state$covMat
  N<- do.call(invlink, list(eta))
  Nr<- N - R
  V<- grad^2 * diag(covMat)
  se.N<- sqrt(V)
  cv.N<- se.N/N
  z <- qt((1-CI.level)/2, length(x$catch) - 2, lower.tail = FALSE)
  if(CI.calc == "norm") {
    C <- z*se.N
    lcl.N <- max(R, N - C)
    ucl.N <- N + C
    lcl.Nr<- max(0, (N-R) - C)
    ucl.Nr<- (N - R) + C
  }
  else if(CI.calc == "lnorm") {
    C <- exp(z * sqrt(log(1 + ((se.N^2)/((N - R)^2)))))
    lcl.N <- R + (N - R)/C
    ucl.N <- R + (N - R)*C
    lcl.Nr<- (N - R)/C
    ucl.Nr<- (N - R)*C
  }
  else if(CI.calc == "boot") {
    Bout<- rep(NA, nboot)
    k<- obj$estimates$catch$estimates
    p0<- plogis(k)
    pb<- 1-(1 - p0)^x$effort
    pi<- removalPiFun(pb)
    nllFun<- nllFun(obj)
    for(j in 1:nboot) {
      newc <- rmultinom(1, N, c(pi, 1 - sum(pi)))
      Bdata <- data.frame(catch = newc[-length(newc)],effort = x$effort)
      Bdata$cpue <- Bdata$catch/Bdata$effort
      Bdata$cumcatch<- cumsum(Bdata$catch) - Bdata$catch
      cf <- coef(lm(cpue ~ cumcatch, data = Bdata))
      cstart<- log(-cf[2])
      nstart<- log(-cf[1]/cf[2])
      starts<- c(nstart,cstart)
      if(sum(Bdata$catch) > 0) {
        m<- optim(starts, nllFun, x=Bdata)
        Bout[j] <- do.call(invlink, list(m$par[1]))
      }
      else Bout[j]<- 0
    }
    se.N <- sd(Bout)
    C <- z*se.N
    lcl.N <- max(R, N - C)
    ucl.N <- N + C
    lcl.Nr<- max(0, (N-R) - C)
    ucl.Nr<- (N - R) + C
  }
  else stop("Unknown CI method")


  bigN<- data.frame(N=round(N), se=round(se.N,1), lcl=round(lcl.N,1), ucl=round(ucl.N,1))
  littleN<- data.frame(N = round(Nr), se=round(se.N,1),lcl=round(lcl.Nr,1), ucl=round(ucl.Nr,1))
  row.names(bigN)<- "Total"
  row.names(littleN)<- "Residual"
  list(Nhat=bigN, Nresid=littleN)
}

#' @rdname calcN
#' @export
calcN.efitGPlist<- function(obj, CI.level=0.95, CI.calc = c("norm","lnorm","boot"), nboot=500, ...){
  out<- lapply(obj, calcN.efitGP, CI.level=CI.level, CI.calc=CI.calc, nboot=nboot, ...)
  out
}


#' @rdname calcN
#' @export
calcN.efitMS<- function(obj, newdata, off.set=NULL, CI.level=0.95, npost=500, ...) {
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
  T <- origdata$numPrimary
  J <- ncol(origdata$y)/T
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
  # Overall Occupancy estimate (random effects)
  re<- raneffects(obj)
  pp<- postSamples(re, npost)
  pp.sum<- apply(pp, c(2,3), mean)
  Nhat<- apply(pp.sum, 1, mean)
  seN<- apply(pp.sum, 1, sd)
  lwr<- apply(pp.sum, 1, quantile, (1-CI.level)/2)
  upr<- apply(pp.sum, 1, quantile, 1-((1-CI.level)/2))
  df<- data.frame(N = round(Nhat,2),se=round(seN,2),
                       lcl=round(lwr,2), ucl=round(upr,2),.season = 1:T)
  list(cellpreds=cellpreds, Nhat=df)
}


#' @rdname calcN
#' @export
calcN.efitMNO<- function(obj, newdata, off.set=NULL, CI.level=0.95, npost=500, ...) {
  if(missing(newdata) || is.null(newdata)) {
    origdata <- obj$data
    M <- numSites(origdata)
    if(is.null(siteCovs(origdata))) {
      newdata <- data.frame(Intercept = rep(1, M))
    } else {
      newdata <- siteCovs(origdata)
    }
  }
  T<- obj$data$numPrimary
  num.removed<- obj$data$num.removed
  design <- getDesign(obj, newdata)
  X<- design$X
  M<- nrow(X)
  sites<- design$retained.sites
  mixture<- obj$mixture
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj, "state")
  covMat<- vcov(obj, "state")
  if(mixture == "ZIP") {
    logit.psi<- coef(obj, "zeroinfl")
    logit.psi.var<- vcov(obj, "zeroinfl")
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
  seN<- sqrt(varN)
  cv<- sqrt(varN)/Nhat
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr<- Nhat*z
  upr<- Nhat/z
  # Residual N estimate (random effects)
  re<- raneffects(obj)
  pp<- postSamples(re, npost)
  pp.sum<- apply(pp, c(2,3), sum)
  Nr<- apply(pp.sum, 1, mean)
  seR<- apply(pp.sum, 1, sd)
  lwr1<- apply(pp.sum, 1, quantile, (1-CI.level)/2)
  upr1<- apply(pp.sum, 1, quantile, 1-((1-CI.level)/2))
  bigN<- data.frame(N=round(Nhat),se=round(seN,1), lcl=round(lwr), ucl=round(upr))
  littleN<- data.frame(N = round(Nr),se=round(seR,1), .season = 1:T, lcl=round(lwr1), ucl=round(upr1))

  Nresid<- Nr[T] - num.removed[T]
  cv<- seR[T]/Nresid
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr2<- Nresid*z
  upr2<- Nresid/z
  residN <- data.frame(N=round(Nresid), se=round(seR[T],1), lcl=round(lwr2), ucl=round(upr2))
  row.names(bigN)<- "Initial"
  row.names(residN)<- "Residual"
  list(cellpreds=cellpreds, Nhat=bigN, Nseason=littleN, Nresid=residN)
}

#' @rdname calcN
#' @export
calcN.efitMNS<- function(obj, newdata, off.set=NULL, CI.level=0.95, ...) {
  # get original data
  origdata <- obj$data
  T <- origdata$numPrimary
  delta <- origdata$delta
  if(missing(newdata) || is.null(newdata)) {
    M <- numSites(origdata)
    season <- data.frame(.season = as.factor(rep(1:T, each = M)))
    trend<- data.frame(.trend = rep(delta, each = M))
    if(is.null(siteCovs(origdata))) {
      scovs <- data.frame(Intercept = rep(1, M))
    } else {
      scovs <- siteCovs(origdata)
    }
    newdata <- cbind(season, trend, scovs[rep(1:M, T),,drop=FALSE])
  } else {
    M<- nrow(newdata)
    season <- data.frame(.season = as.factor(rep(1:T, each = M)))
    trend<- data.frame(.trend = rep(delta, each = M))
    newdata <- cbind(season, trend, newdata[rep(1:M, T),,drop=FALSE])
  }

  design <- getDesign(obj, newdata)
  X<- design$X
  M<- nrow(X)
  num.removed <- obj$data$num.removed
  sites<- design$retained.sites
  mixture<- obj$mixture
  invlink = obj$estimates$state$invlink
  invlinkGrad = obj$estimates$state$invlinkGrad
  estimates<- coef(obj, "state")
  covMat<- vcov(obj, "state")
  if(mixture == "ZIP") {
    logit.psi<- coef(obj, "zeroinfl")
    logit.psi.var<- vcov(obj, "zeroinfl")
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
    cellpreds<- data.frame(N = ests * off.set, se = sqrt(v), site = sites, .season = season)
  }
  else{
    eta <- as.vector(X %*% estimates)
    vc <- rowSums((X %*% covMat) * X) # equivalent to diag(X %*% covMat %*% t(X))
    ests<- do.call(invlink,list(eta))
    grad <- do.call(invlinkGrad,list(eta))
    v<- (grad * off.set)^2 * vc
    cellpreds<- data.frame(N = ests * off.set, se = sqrt(v), site = sites, .season = season)
  }
  # overall estimate
  bigN<- matrix(NA, nrow=T, ncol=5)
  for(i in 1:T) {
    off.seas<- off.set
    off.seas[which(season != i)]<- 0
    Nhat<- off.seas %*% ests
    gradN<- off.seas %*% (grad * X)
    varN<- gradN %*% covMat %*% t(gradN)
    seN<- sqrt(varN)
    cv<- sqrt(varN)/Nhat
    z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
    lwr<- Nhat*z
    upr<- Nhat/z
    bigN[i, ]<- c(round(Nhat), i, round(seN,1), round(lwr), round(upr))
  }
  bigN<- as.data.frame(bigN)
  Nresid<- Nhat - num.removed[T]
  cv<- seN/Nresid
  z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2)))
  lwr2<- Nresid*z
  upr2<- Nresid/z
  residN <- data.frame(N=round(Nresid), se=round(seN,1), lcl=round(lwr2), ucl=round(upr2))
  row.names(residN)<- "Residual"
  names(bigN)<- c("Nhat",".season","se","lwr","upr")
  list(cellpreds=cellpreds, Nhat=bigN, Nresid=residN)
}

#-------------------------------------------
#' @rdname se
#' @export
se.efit<- function(obj, type, ...){
  if(is.null(type)) stop("estimate type required")
  v<- obj$estimates[[type]]$covMat
  sqrt(diag(v))
}
#--------------------------------------------------------------------
# calcP methods
#--------------------------------------------------------------------
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

#' @rdname calcP
#' @export
calcP.efitGP<- function(obj, na.rm = TRUE) {
  emf<- obj$data
  eff<- emf$effort
  pars <- coef(obj, type = "catch")
  p0<- plogis(pars)
  p <- 1 - (1-p0)^eff
  pi <- removalPiFun(p)
  return(pi)
}

#' @rdname calcP
#' @export
calcP.efitMNO<- function(obj, na.rm = TRUE) {
# multinomial open
  data <- obj$data
  D <- getDesign(data, obj$lamformula, obj$gamformula, obj$omformula, obj$detformula,
                 obj$iotaformula, na.rm=na.rm)
  detparms <- coef(obj, 'det')
  Xp.offset <- D$Xp.offset
  if(is.null(Xp.offset)) Xp.offset <- rep(0, nrow(D$Xp))
  plong <- plogis(D$Xp %*% detparms + Xp.offset)

  M <- nrow(D$y)
  T <- data$numPrimary
  J <- ncol(D$y) / T

  pmat <- aperm(array(plong, c(J,T,M)), c(3,1,2))

  pout <- array(NA, c(M,J,T))
  for (t in 1:T){
    pout[,,t] <- do.call(data$piFun, list(p=pmat[,,t]))
  }
  matrix(aperm(pout,c(2,3,1)), M, J*T, byrow=TRUE)

}

#' @rdname calcP
#' @export
calcP.efitMS<- function(obj, na.rm = TRUE) {
  data <- obj$data
  detParms <- coef(obj, 'det')
  D <- getDesign(data, obj$lamformula, obj$gamformula, obj$epsformula, obj$detformula, na.rm=na.rm)
  y <- D$y
  V <- D$V

  M <- nrow(y)	# M <- nrow(X.it)
  nY <- data$numPrimary
  J <- numY(data)/nY

  p <- plogis(V %*% detParms)
  p <- array(p, c(J, nY, M))
  p <- aperm(p, c(3, 1, 2))
  p <- matrix(p, nrow=M)
  return(p)
}
#-------------------------------------------------------
# update methods
#-------------------------------------------------------
#' update
#'
#' @description updates efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
update.efitMS<- function(obj, lamformula., gamformula., epsformula.,
                          detformula., ..., evaluate = TRUE) {

  if(is.null(call <- obj$call)) stop("need an object with call component")
  lamformula <- as.formula(call[['lamformula']])
  gamformula <- as.formula(call[['gamformula']])
  epsformula <- as.formula(call[['epsformula']])
  detformula <- as.formula(call[['detformula']])

  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(lamformula.)) {
    uplamformula <- update.formula(lamformula, lamformula.)
    call[['lamformula']] <- uplamformula
  }
  if (!missing(gamformula.)) {
    upgamformula <- update.formula(gamformula, gamformula.)
    call[['gamformula']] <- upgamformula
  }
  if (!missing(epsformula.)) {
    upepsformula <- update.formula(epsformula, epsformula.)
    call[['epsformula']] <- upepsformula
  }
  if (!missing(detformula.)) {
    updetformula <- update.formula(detformula, detformula.)
    call[['detformula']] <- updetformula
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
}



#' update
#'
#' @description updates efit model objects.
#'
#' @param obj A fitted model object.
#'
#' @export
#'
update.efitMNO<- function(obj, lamformula., detformula., gamformula.,
                          omformula., ..., evaluate = TRUE) {

  if(is.null(call <- obj$call)) stop("need an object with call component")
  lamformula <- as.formula(call[['lamformula']])
  detformula <- as.formula(call[['detformula']])
  gamformula <- as.formula(call[['gamformula']])
  omformula <- as.formula(call[['omformula']])

  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(lamformula.)) {
    upLamformula <- update.formula(lamformula,lamformula.)
    call[['lamformula']] <- upLamformula
  }
  if (!missing(gamformula.)) {
    upGamformula <- update.formula(gamformula, gamformula.)
    call[['gamformula']] <- upGamformula
  }
  if (!missing(omformula.)) {
    upOmformula <- update.formula(omformula, omformula.)
    call[['omformula']] <- upOmformula
  }
  if (!missing(detformula.)) {
    upDetformula <- update.formula(detformula, detformula.)
    call[['detformula']] <- upDetformula
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
}


#-------------------------------------------------------
#' @rdname traject
#' @export
traject.efitMS<- function(obj, mean=TRUE, ...){

  backward <- function(detParams, phis) {
    beta <- array(NA, c(K + 1, nY, M))
    for (i in 1:M) {
      backP <- rep(1, K + 1)
      for (t in nY:1) {

        beta[, t, i] <- backP

        detVec <- rep(1, K + 1)
        for (j in 1:J) {
          if(!is.na(y.arr[i,t,j])) {
            mp <- V.arr[,,i,t,j] %*% detParams
            detVecObs <- gSingleDetVec(y.arr[i,t,j], mp)

            detVec <- detVec * detVecObs
          }
        }
        if (t > 1)
          backP <- t(phis[,,t-1,i]) %*% (detVec * backP)
      }
    }
    return(beta)
  }

  data <- obj$data
  nY <- data$numPrimary
  ests <- obj$estimates

  lamformula <- obj$lamformula
  gamformula <- obj$gamformula
  epsformula <- obj$epsformula
  detformula <- obj$detformula
  K <- 1
  designMats <- getDesign(data, lamformula, gamformula, epsformula, detformula)
  V.itjk <- designMats$V
  X.it.gam <- designMats$X.gam
  X.it.eps <- designMats$X.eps
  W.i <- designMats$W

  detParms <- colnames(V.itjk)
  gamParms <- colnames(X.it.gam)
  epsParms <- colnames(X.it.eps)
  psiParms <- colnames(W.i)

  psiParams <- ests$state$estimates
  colParams <- ests$col$estimates
  extParams <- ests$ext$estimates
  detParams <- ests$det$estimates

  y <- designMats$y
  M <- nrow(y)
  J <- ncol(y)/nY

  ## remove final year from X.it
  X.it.gam <- as.matrix(X.it.gam[-seq(nY,M*nY,by=nY),])
  X.it.eps <- as.matrix(X.it.eps[-seq(nY,M*nY,by=nY),])

  X.gam <- X.it.gam %x% c(-1,1)
  X.eps <- X.it.eps %x% c(-1,1)
  phis <- array(NA,c(2,2,nY-1,M))

  psis <- plogis(W.i %*% psiParams)

  ## computed projected estimates
  phis[,1,,] <- plogis(X.gam %*% colParams)
  phis[,2,,] <- plogis(X.eps %*% -extParams)

  projected <- array(NA, c(2, nY, M))
  projected[1,1,] <- 1 - psis
  projected[2,1,] <- psis
  for(i in 1:M) {
    for(t in 2:nY) {
      projected[,t,i] <- phis[,,t-1,i] %*% projected[,t-1,i]
    }
  }
  return(t(projected[2,,]))
}



#' @rdname traject
#' @export
#'
traject.efitMNO <- function(obj, na.rm=FALSE) {
  dynamics <- obj$dynamics
  mixture <- obj$mixture
  fix <- obj$fix
  immigration <- obj$immigration
  data <- obj$data
  D <- getDesign(data, obj$lamformula, obj$gamformula, obj$omformula, obj$detformula,
                 obj$iotaformula, na.rm=na.rm)
  Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xiota <- D$Xiota
  Xlam.offset <- D$Xlam.offset; Xgam.offset <- D$Xgam.offset
  Xom.offset <- D$Xom.offset
  Xiota.offset <- D$Xiota.offset
  delta <- D$delta

  y <- D$y
  M <- nrow(y)
  T <- data$numPrimary
  J <- ncol(y) / T

  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
  if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
  if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
  if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

  lambda <- exp(Xlam %*% coef(obj, 'state') + Xlam.offset)
  if(identical(mixture, "ZIP")) {
    psi <- plogis(coef(obj, type="psi"))
    lambda <- (1-psi)*lambda
  }
  if (fix == 'omega'){
    omega <- matrix(1, M, T-1)
  } else if(!identical(dynamics, "trend")) {
      omega <- matrix(plogis(Xom %*% coef(obj, 'omega') + Xom.offset),
                      M, T-1, byrow=TRUE)
  }
  if(fix == "gamma"){
    gamma <- matrix(0, M, T-1)
  } else {
    gamma <- matrix(exp(Xgam %*% coef(obj, 'gamma') + Xgam.offset),
                    M, T-1, byrow=TRUE)
  }
  if(immigration)
    iota <- matrix(exp(Xiota %*% coef(obj, 'iota') + Xiota.offset),
                   M, T-1, byrow=TRUE)
  else
    iota <- matrix(0, M, T-1)

  N <- matrix(NA, M, T)
  for(i in 1:M) {
    N[i, 1] <- lambda[i]
    if(delta[i, 1] > 1) {
      for(d in 2:delta[i ,1]) {
        if(identical(dynamics, "autoreg"))
          N[i, 1] <- N[i, 1] * (omega[i,1] + gamma[i, 1]) + iota[i, 1]
        else if(identical(dynamics, "trend"))
          N[i,1] <- N[i,1] * gamma[i,1] + iota[i, 1]
        else
          N[i,1] <- N[i,1] * omega[i,1] + gamma[i,1]
      }
    }
    for(t in 2:T) {
      if(identical(dynamics, "autoreg"))
        N[i, t] <- N[i, t-1] * (omega[i, t-1] + gamma[i, t-1]) +
          iota[i, t-1]
      else if(identical(dynamics, "trend"))
        N[i,t] <- N[i,t-1] * gamma[i,t-1] + iota[i, t-1]
      else
        N[i,t] <- N[i,t-1] * omega[i,t-1] + gamma[i,t-1]
      if(delta[i, t] > 1) {
        for(d in 2:delta[i, t]) {
          if(identical(dynamics, "autoreg"))
            N[i, t] <- N[i, t] * (omega[i, t-1] + gamma[i, t-1]) +
              iota[i, t-1]
          else if(identical(dynamics, "trend"))
            N[i, t] <- N[i, t] * gamma[i, t-1] + iota[i, t-1]
          else
            N[i,t] <- N[i,t] * omega[i, t-1] + gamma[i, t-1]
        }
      }
    }
  }
  N
}

