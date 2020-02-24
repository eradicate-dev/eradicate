
#' remPois
#'
#' @name remPois
#'
#' @description
#' \code{remPois} fits the multinomial-Poisson removal model to data from
#' a number of primary periods where individuals removed are recorded for each site.
#'
#' @usage remPois(lamformula, detformula, data, starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the removal detection component.  Only
#'  site-level covariates are allowed for the removal detection component.
#'  This differs from the similar model in \code{unmarked}. Currently
#'  an intercept-only model is assumed for the occupancy component.
#' @param data A \code{eFrameR} object containing the response (counts)
#'  and site-level covariates. see \code{\link{eFrameR}} for how to format
#'  the required data.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  data(snc)
#'  emf <- eFrameR(y=removed, siteCovs=site.df)
#'  mod <- remPois(~1, ~1, data=emf)
#'  Nhat<- calcN(mod, ncells=55)
#'
#' @export
#'
remPois <- function(lamformula, detformula, data, starts, method = "BFGS", se = TRUE, ...) {

    if(!is(data, "eFrameR"))
		    stop("Data is not a eFrameR.")

    designMats <- getDesign(data, lamformula, detformula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
        }
    J <- ncol(y)
    R <- ncol(y)
    M <- nrow(y)
    piFun <- data$piFun

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    nP <- nDP + nAP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(y)
    navec <- is.na(yvec)

    nll <- function(parms) {
        lambda <- exp(X %*% parms[1 : nAP] + X.offset)
        p <- plogis(V %*% parms[(nAP + 1) : nP] + V.offset)
        p.matrix <- matrix(p, M, R, byrow = TRUE)
        pi <- do.call(piFun, list(p = p.matrix))
        logLikeSite <- dpois(y, matrix(lambda, M, J) * pi, log = TRUE)
        logLikeSite[navec] <- 0
        -sum(logLikeSite)
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
    names(ests) <- c(lamParms, detParms)

    stateName <- "Abundance"

    stateEstimates <- list(name = stateName,
                           short.name = "lambda",
                           estimates = ests[1:nAP],
                           covMat = as.matrix(
                           covMat[1:nAP,1:nAP]),
                           invlink = "exp",
                           invlinkGrad = "exp")

    detEstimates <- list(name = "Detection", short.name = "p",
                         estimates = ests[(nAP + 1) : nP],
                         covMat = as.matrix(covMat[(nAP + 1):nP, (nAP + 1):nP]),
                         invlink = "logistic",
                         invlinkGrad = "logistic.grad")

    efit <- list(fitType = "removal Poisson",
        call = match.call(), formula = formula,state=stateEstimates,det=detEstimates,
        sitesRemoved = designMats$removed.sites, AIC = fmAIC, opt = opt,
        negLogLike = fm$value, nllFun = nll)
    class(efit) <- c('efitN','efit','list')
    return(efit)
}
