
#' calcN
#'
#' @description
#' \code{calcN} estimates abundance for a defined region from a fitted model.
#' The user must supply the number of spatial units (cells) or a data.frame
#' of covariate values for each spatial unit.
#'
#' @param obj A fitted model object.
#' @param covs A \code{data.frame} of covariates for each spatial unit.
#'    There must by covariates for every parameter in \code{obj}.
#'    if no covariates are used then \code{covs} should be a data.frame
#'    with one column indicating cell number.
#' @param off.set Either a scalar offset value to apply to each spatial unit
#' for prediction (e.g. cell area) or a vector of the same length as \code{nrow(covs)}.
#'
#' @return a \code{data.frame} giving the predictions for each spatial unit
#'  as well as the overall abundance estimate for the region with associated
#'  SE and confidence intervals.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  site.df<- san_nic_pre$traps
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  mod <- nmix(~1, ~1, data=emf)
#'  Nhat<- calcN(mod, site.df)
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
  # State
    ests <- object$state$estimates
    SEs <- SE(object, "state")
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
        !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object$state$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit")
    cat(object$state$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
    cat(" ","\n")
  # detection
    ests <- object$det$estimates
    SEs <- SE(object, "det")
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
      !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object$det$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit")
    cat(object$det$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
}
#' @rdname summary.efit
#' @export
summary.efitREST<- function(object, ...)
{
  # State
    ests <- object$state$estimates
    SEs <- SE(object, "state")
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
        !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object$state$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit")
    cat(object$state$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
    cat(" ","\n")
  # staying time
    ests <- object$stay$estimates
    SEs <- SE(object, "stay")
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
      !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object$stay$invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit")
    cat(object$stay$name, " (", link, "-scale)", ":\n", sep="")
    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
}

# Compute linear combinations of estimates using coefficients

linearComb.efit<- function(obj, coefficients, off.set = NULL, ...)
{
    estimates<- obj$state$estimates
    covMat<- obj$state$covMat
    if(!is(coefficients, "matrix"))
        coefficients <- t(as.matrix(coefficients))
    if(ncol(coefficients) != length(estimates)) stop("error - wrong number oc covariates")
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

#' @rdname calcN
#' @export
calcN.efit<- function(obj, covs, off.set=NULL, CI.level=0.95, ...) {
  design <- getDesign(obj, covs)
  X<- design$X
  M<- nrow(X)
  if(!is.null(off.set) & length(off.set) == 1) off.set<- rep(off.set, M)
  lc<- linearComb(obj, coefficients=X, off.set=off.set)
  est<- backTransform(lc)
  V<- est$covMat
  Nhat<- sum(est$estimates)
  varN<- sum(est$covMat)
  seN<- sqrt(varN)
  z <- qnorm((1-CI.level)/2, lower.tail = FALSE)
  lwr<- Nhat - seN*z
  upr<- Nhat + seN*z
  cv<- sqrt(varN)/Nhat
  za <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2))) # asymptotic CI
  bigN<- data.frame(Nhat=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1),
                    lcl.ln=round(Nhat*za,1), ucl.ln=round(Nhat/za,1))
  list(cellpreds=est$estimates, Nhat=bigN)
}

#' @rdname calcN
#' @export
calcN.efitM<- function(obj, covs, off.set=NULL, CI.level=0.95, ...) {
  design <- getDesign(obj, covs)
  X<- design$X
  M<- nrow(X)
  if(!is.null(off.set) & length(off.set) == 1) off.set<- rep(off.set, M)
  lc<- linearComb(obj, coefficients=X, off.set=off.set)
  est<- backTransform(lc)
  V<- est$covMat
  Pocc<- sum(est$estimates)/ M # proportion of sites occupied
  varN<- sum(est$covMat) * (1/M^2) # delta method VAR
  seN<- sqrt(varN)
  z <- qnorm((1-CI.level)/2, lower.tail = FALSE)
  lwr<- Pocc - seN*z
  upr<- Pocc + seN*z
  bigN<- data.frame(Pocc=round(Pocc,2),se=round(seN,3), lcl=round(lwr,2), ucl=round(upr,2))
  list(cellpreds=est$estimates, Occ=bigN)
}

#' @rdname SE
#' @export
SE.efit<- function(obj, type=c("state","det"), ...){
  type <- match.arg(type, c("state", "det"))
  if(identical(type,"state"))
    v<- obj$state$covMat
  else if(identical(type,"det"))
    v<- obj$det$covMat
    sqrt(diag(v))
}

#' @rdname SE
#' @export
SE.efitREST<- function(obj, type=c("state","stay"), ...){
  type <- match.arg(type, c("state", "stay"))
  if(identical(type,"state"))
    v<- obj$state$covMat
  else if(identical(type,"stay"))
    v<- obj$stay$covMat
    sqrt(diag(v))
}
