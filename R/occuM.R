
#' occuM
#'
#' @name occuM
#'
#' @description
#' \code{occuM} fits the occupancy model of McKenzie et al. (2002).
#'
#' @usage occuM(lamformula, detformula, data, knownocc, starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the detection component.  Only
#'  site-level covariates are allowed for the detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrame} object containing the response (counts)
#'  and site-level covariates. see \code{\link{eFrame}} for how to format
#'  the required data.
#' @param knownocc Vector of row numbers of sites that are known to be occupied.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  mod <- occuM(~1, ~1, data=emf)
#'  Nhat<- calcN(mod, ncells=55)
#'
#' @export
#'
occuM<- function(lamformula, detformula, data, knownOcc = numeric(0), starts,
                 method = "BFGS", se = TRUE, ...) {

    if(!is(data, "eFrame")) stop("Data is not an eFrame")

    designMats <- getDesign(data, lamformula, detformula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    removed <- designMats$removed.sites
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if(is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    y <- truncateToBinary(y)
    J <- ncol(y)
    M <- nrow(y)

    ## convert knownOcc to logical so we can correctly to handle NAs.
    knownOccLog <- rep(FALSE, numSites(data))
    knownOccLog[knownOcc] <- TRUE
    if(length(removed)>0)
        knownOccLog <- knownOccLog[-removed]

    occParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nOP <- ncol(X)

    nP <- nDP + nOP
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(t(y))
    navec <- is.na(yvec)
    nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

    nll <- function(params) {
        psi <- plogis(X %*% params[1 : nOP] + X.offset)
        psi[knownOccLog] <- 1
        pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
        cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
        cp[navec] <- 1 # so that NA's don't modify likelihood
        cpmat <- matrix(cp, M, J, byrow = TRUE) #
        loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
        -sum(loglik)
        }


    if(missing(starts)) starts <- rep(0, nP)
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
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)

    state <- list(name = "Occupancy", short.name = "psi",
                              estimates = ests[1:nOP],
                              covMat = as.matrix(covMat[1:nOP,1:nOP]),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad")

    det <- list(name = "Detection", short.name = "p",
                            estimates = ests[(nOP + 1) : nP],
                            covMat = as.matrix(covMat[(nOP + 1) : nP,
                                                      (nOP + 1) : nP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")


    efit <- list( fitType = "occuM", call = match.call(),
                 formula = formula, state=state, det=det,
                 sitesRemoved = designMats$removed.sites,
                 AIC = fmAIC, opt = opt, negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog)
    class(efit) <- c('efitN','efit','list')
    return(efit)
}
