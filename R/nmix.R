
#' nmix
#'
#' @name nmix
#'
#' @description
#' \code{nmix} fits the N-mixture model of Royle et al (2004)
#'
#' @usage nmix(lamformula, detformula, data, K, mixture=c("P", "NB"),
#' starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the detection component.  Only
#'  site-level covariates are allowed for the detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrame} object containing the response (counts)
#'  and site-level covariates. see \code{\link{eFrame}} for how to format
#'  the required data.
#' @param K Integer upper index of integration for abundance. This should be
#'  set high enough so that it does not affect the parameter estimates. Note
#'  that computation time will increase with K
#' @param mixture Distribution modelfor the latent abundance, either Poisson (P)
#'  or Negative-binomial (NB).
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  mod <- nmix(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
nmix <- function(lamformula, detformula, data, K, mixture = c("P", "NB"), starts,
                   method = "BFGS", se = TRUE, ...)
{
    mixture <- match.arg(mixture, c("P", "NB"))
    if(!is(data, "eFrame")) stop("Data is not an eFrame")

    designMats <- getDesign(data, lamformula, detformula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    NAmat <- is.na(y)

    J <- ncol(y)
    M <- nrow(y)

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)

    if(missing(K)) {
        K <- max(y, na.rm = TRUE) + 100
        warning("K was not specified and was set to ", K, ".")
    }
    if(K <= max(y, na.rm = TRUE))
        stop("specified K is too small. Try a value larger than any observation")
    k <- 0:K
    lk <- K+1
    M <- nrow(y)
    J <- ncol(y)
    k.ik <- rep(k, M)
    k.ijk <- rep(k, M*J)

    nP <- nAP + nDP + (mixture != "P")
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    y.ij <- as.numeric(t(y))
    y.ijk <- rep(y.ij, each = K + 1)
    navec <- is.na(y.ijk)
    ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
    ijk.to.ikj <- with(ijk, order(i, k, j))
    nll <- function(parms) {
        theta.i <- exp(X %*% parms[1 : nAP] + X.offset)
        p.ij <- plogis(V %*% parms[(nAP + 1) : (nAP + nDP)] + V.offset)
        theta.ik <- rep(theta.i, each = K + 1)
        p.ijk <- rep(p.ij, each = K + 1)

        bin.ijk <- dbinom(y.ijk,k.ijk,p.ijk)
        bin.ijk[which(is.na(bin.ijk))] <- 1
        bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (K + 1), J,
                             byrow = TRUE)
        g.ik <- rowProds(bin.ik.mat)

        if(identical(mixture,"P")) {
            f.ik <- dpois(k.ik,theta.ik)
        }
        else if (identical(mixture,"NB")){
            f.ik <- dnbinom(k.ik, mu = theta.ik, size = exp(parms[nP]))
        }
        dens.i.mat <- matrix(f.ik * g.ik, M, K + 1, byrow = TRUE)
        dens.i <- rowSums(dens.i.mat)  # sum over the K

        -sum(log(dens.i))
        }


    if(missing(starts)) starts <- rep(0, nP)
    fm <- optim(starts, nll, method=method, hessian=se, ...)
    opt <- fm
    ests <- fm$par
    nbParm <- switch(mixture,
                     NB = "alpha",
                     P = character(0))
    names(ests) <- c(lamParms, detParms, nbParm)
    if(se) {
        tryCatch(covMat <- solve(fm$hessian), error=function(x)
                 stop(simpleError("Hessian is singular.  Try using fewer covariates.")))
    } else {
        covMat <- matrix(NA, nP, nP)
    }
    fmAIC <- 2 * fm$value + 2 * nP

    typeNames<- c("state","det")
    stateEstimates <- list(name="Abundance", short.name="lam",
        estimates = ests[1:nAP],
        covMat = as.matrix(covMat[1:nAP,1:nAP]),
	invlink = "exp", invlinkGrad = "exp")

    detEstimates <- list(name = "Detection", short.name = "p",
        estimates = ests[(nAP + 1) : (nAP + nDP)],
        covMat = as.matrix(covMat[(nAP + 1):(nAP + nDP),
                                  (nAP + 1):(nAP + nDP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")


    if(identical(mixture,"NB")) {
        stateEstimates$alpha <- list(name="Dispersion",
            short.name = "alpha", estimates = ests[nP],
            covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
            invlinkGrad = "exp")
    }

    efit <- list(fitType="nmix", call=match.call(), types=typeNames,
                 lamformula = lamformula, detformula=detformula,
                 sitesRemoved = designMats$removed.sites,
                 state=stateEstimates, det=detEstimates, AIC = fmAIC, opt = opt,
                 negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
                 data = data)
    class(efit) <- c('efit','list')
    return(efit)
}
