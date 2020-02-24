
#' calcN
#'
#' @description
#' \code{calcN} estimates abundance for a defined region from a fitted model.
#' The user must supply the number of spatial units (cells) or a data.frame
#' of covariate values for each spatial unit.
#'
#' @param obj A fitted model object.
#' @param cells The number of spatial units used for prediction
#' @param covs A \code{data.frame} of covariates for each spatial unit.
#'    There must by covariates for every parameter in \code{obj}.
#'    Either one of \code{cells} or \code{covs} must be supplied.
#'
#' @return a \code{data.frame} giving the overall abundance estimate for
#'  the region with associated SE and confidence intervals.
#'
#' @examples
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  mod <- nmix(~1, ~1, data=emf)
#'  Nhat<- calcN(mod, ncells=55)
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
#' \code{summary} summarises an efit object.
#'
#' @param obj A fitted model object.
#'
#' @return a \code{data.frame} giving the estimated parameters for the
#'  model in \code{obj} (on the link scale) and associated SE and 95% CI.
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
#' @describeIn summary.efit
#'
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

linearComb.efit<- function(obj, coefficients = NULL, offset = NULL, ...)
{
    estimates<- obj$state$estimates
    covMat<- obj$state$covMat
    if(is.null(coefficients) & length(estimates == 1))
      coefficients<- t(as.matrix(1))
    if(!is(coefficients, "matrix"))
        coefficients <- t(as.matrix(coefficients))
    if(ncol(coefficients) != length(estimates)) stop("error - model contains covariates")
    if (is.null(offset))
        offset <- rep(0, nrow(coefficients))
    e <- as.vector(coefficients %*% estimates) + offset
    v <- coefficients %*% covMat %*% t(coefficients)

    umlc <- list(estimates = e, covMat = v, coefficients = coefficients,
                 invlink = obj$state$invlink, invlinkGrad = obj$state$invlinkGrad)
    class(umlc)<- c("efit", class(umlc))
    umlc
}


backTransform.efit<- function(obj, ...) {

            ## In general, MV delta method is Var=J*Sigma*J^T where J is Jacobian
            ## In this case, J is diagonal with elements = gradient
            ## Reduces to scaling the rows then columns of Sigma by the gradient
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
calcN.default<- function(obj, ncells=NULL, covs=NULL, CI.level=0.95, ...) {

  if(is.null(ncells) & is.null(covs)) stop("one of ncells or covs must be supplied")
  if(!is.null(ncells) & is.null(covs)) {
    lc<- linearComb(obj)
  } else if(is.null(ncells) & !is.null(covs)) {
    lc<- linearComb(obj, coefficients=covs)
  }
  else stop("problem")
  est<- backTransform(lc)
  if(length(est$estimates == 1) & (!is.null(ncells))) {
    V<- as.vector(est$covMat)
    Nhat<- est$estimates * ncells
    varN<- ncells^2 * V
    seN<- sqrt(varN)
    z <- qnorm((1-CI.level)/2, lower.tail = FALSE)
    lwr<- Nhat - seN*z
    upr<- Nhat + seN*z
    data.frame(Nhat=round(Nhat,1),se=round(seN,1), lcl=round(lwr,1), ucl=round(upr,1))
  }
  else if(length(est$estimates) > 1 & !is.null(covs)) {
    V<- est$covMat
    Nhat<- sum(est$estimates)
    varN<- sum(est$covMat)
    seN<- sqrt(varN)
    cv<- sqrt(varN)/Nhat
    z <- exp(qnorm((1-CI.level)/2) * sqrt(log(1+cv^2))) # asymptotic CI
    bigN<- data.frame(Nhat=round(Nhat,1),se=round(seN,1), lcl=round(Nhat*z,1), ucl=round(Nhat/z,1))
    list(cellpreds=est$estimates, Nhat=bigN)
  }
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
