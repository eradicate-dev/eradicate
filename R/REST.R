#' REST
#'
#' @name REST
#'
#' @description
#' \code{REST} fits the random encounter and staying time model of
#' Nakashima et al (2018), which models the relationship between camera
#' encounter rate, staying time and density. The model assume that cameras
#' sample a viewshed of area \code{A} with perfect detection.
#'
#' @usage REST(formula, data, starts, method="BFGS", se=TRUE, ...)
#'
#' @param formula formula for animal density.
#' @param data A \code{eFrameREST} object containing the number of encounters
#' for each camera and the staying times. see \code{\link{eFrameREST}} for how to format
#'  the required data.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  y<- san_nic_rest$y
#'  stay<- san_nic_rest$stay
#'  cens<- san_nic_rest$cens
#'  area<- san_nic_rest$area
#'  active_hours<- san_nic_rest$active_hours
#'
#'  emf <- eFrameREST(y, stay, cens, area, active_hours)
#'  mod <- REST(~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
REST <- function(lamformula, data, starts, method = "BFGS",
    se = TRUE, ...)
{
    if(!is(data, "eFrameREST"))
		    stop("Data is not a data frame or eFrameREST.")
    designMats <- getDesign(data, lamformula)
    X <- designMats$X
    y <- designMats$y
    stay <- designMats$stay
    cens <- designMats$cens
    eff<- designMats$effort
    A<- designMats$area
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    R <- ncol(y)
    M <- nrow(y)

    lamParms <- colnames(X)
    stayParm <- "(Intercept)"

    nAP <- ncol(X)
    nDP<- 1
    nP <- nDP + nAP

    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(t(y))
    navec <- is.na(yvec)
    evec <- as.numeric(t(eff))

    nll <- function(parms) {
        lambda <- exp(X %*% parms[1 : nAP] + X.offset)
        lambda.ij<- rep(lambda, each = R)
        logrho<- parms[nAP+nDP]
        logL.uc <- dexp(stay[cens==1], exp(logrho), log=TRUE)
        logL.c<- pexp(stay[cens==0], exp(logrho), lower.tail=FALSE, log=TRUE)
        logL1<- sum(logL.uc) + sum(logL.c)

        logmu<- log(A) + log(evec) + logrho + log(lambda.ij)
        logL2<- dpois(yvec, exp(logmu), log=TRUE)
        logL2[navec]<- 0
        ll<- logL1 + sum(logL2)
        (-1)*ll
      }


    if(missing(starts))
        starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    opt <- fm
    if(se) {
        tryCatch(covMat <- solve(fm$hessian),
                 error=function(x) stop(simpleError("Hessian is singular.
                                        Try providing starting values or using fewer covariates.")))
    } else {
        covMat <- matrix(NA, nP, nP)
      }
    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP
    names(ests)<- c(lamParms, stayParm)

    stateEstimates <- list(name = "Abundance",
                                   short.name = "lambda",
                                   estimates = ests[1:nAP],
                                   covMat = as.matrix(covMat[1:nAP,1:nAP]),
                                   invlink = "exp",
                                   invlinkGrad = "exp")

    rhoEstimates <- list(name = "Staying time", short.name = "rho",
                                 estimates = ests[(nAP + nDP)],
                                 covMat = as.matrix(covMat[(nAP + 1):nP, (nAP + 1):nP]),
                                 invlink = "exp",
                                 invlinkGrad = "exp")

    estimates<- list(state=stateEstimates, stay=rhoEstimates)

    efit <- list(fitType = "REST model",
        call = match.call(), lamformula = lamformula, estimates=estimates,
        sitesRemoved = designMats$removed.sites,AIC = fmAIC, opt = opt,
        negLogLike = fm$value, nllFun = nll, data = data)
    class(efit) <- c('efitREST','efit','list')
    return(efit)
}
