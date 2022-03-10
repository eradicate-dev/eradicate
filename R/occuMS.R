
#' occMS
#'
#' @name occuMS
#'
#' @description
#' \code{occMS} fits the dynamic multi-season occupancy model of MacKenzie et al. (2003).
#' This is a port of the \code{colext} function in \code{unmarked}.
#'
#' @usage occMS(lamformula, gamformula, epsformula, detformula,
#' data, starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent occupancy component.
#' @param gamformula formula for the latent colonisation component. If season
#' specific colonisation/extinction parameters are desired, then use the
#' reserved keyword \code{.season}, (i.e. \code{~.season}).
#' @param epsformula formula for the latent survival component.
#' @param detformula formula for the detection component.  Only
#'  site- or \code{.season}-level covariates are allowed for
#'  the detection component. This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrameMS} object containing the detection/non-detection
#' data (0/1) and site-level covariates. see \code{\link{eFrameMS}} for how to format
#'  the required data. count data will get trunctated to (0/1) by \code{occMS()}.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  mod <- occMS(~1, ~1, ~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
occuMS <- function(lamformula = ~ 1, gamformula = ~ 1, epsformula = ~ 1, detformula = ~ 1,
                   data, starts, method = "BFGS", se = TRUE, ...) {

    if(!is(data, "eFrameMS"))
        stop("Data is not a eFrameMS.")

    designMats <- getDesign(data, lamformula, gamformula, epsformula, detformula)
    y <- designMats$y
    V.itjk <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W

    detParms <- colnames(V.itjk)
    gamParms <- colnames(X.it.gam)
    epsParms <- colnames(X.it.eps)
    psiParms <- colnames(W.i)
    parm.names <- c(psiParms, gamParms, epsParms, detParms)
    M <- nrow(y)
    nY <- data$numPrimary
    J <- ncol(y)/nY
    n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))

    ## remove final year from X.it
    X.it.gam <- as.matrix(X.it.gam[-seq(nY,M*nY,by=nY),])
    X.it.eps <- as.matrix(X.it.eps[-seq(nY,M*nY,by=nY),])

    nDP <- length(detParms)
    nGP <- length(gamParms)
    nEP <- length(epsParms)
    nSP <- length(psiParms)
    nDMP <-  1
    K <- 1

    nP <- nDP + nSP + nGP + nEP  # total number of parameters
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    y.itj <- as.numeric(t(y))

    # get ragged array indices
    y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
    J.it <- rowSums(!is.na(y.it))

    V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
    V.arr <- aperm(V.arr, c(2,1,5,4,3))

    y.arr <- array(y.itj, c(J, nY, M))
    y.arr <- aperm(y.arr, c(3:1))

    alpha <- array(NA, c(K + 1, nY, M))

    forward <- function(detParms, phis, psis, storeAlpha = FALSE) {
        negloglike <- 0
        psiSite <- matrix(c(1-psis,psis), K + 1, M, byrow = TRUE)
        mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
        for(t in 1:nY) {
            detVecs <- gDetVecs(y.arr, mp, J.it[seq(from = t, to = length(J.it)-nY+t,
                                                    by=nY)], t)
            psiSite <- psiSite * detVecs
            if(storeAlpha) alpha[,t,] <<- psiSite[,]
            if(t < nY) {
                for(i in 1:M) {
                    psiSite[,i] <- phis[,,t,i] %*% psiSite[,i]
                }
            } else {
                negloglike <- negloglike - sum(log(colSums(psiSite)))
            }
        }
        negloglike
    }

    X.gam <- X.it.gam %x% c(-1,1)
    X.eps <- X.it.eps %x% c(-1,1)
    phis <- array(NA,c(2,2,nY-1,M))

    nll <- function(params) {
        psis <- plogis(W.i %*% params[1:nSP])
        colParams <- params[(nSP + 1) : (nSP + nGP)]
        extParams <- params[(nSP + nGP + 1) : (nSP + nGP + nEP)]
        detParams <- params[(nSP + nGP + nEP + 1) : nP]
        # these are in site-major, year-minor order
        phis[,1,,] <- plogis(X.gam %*% colParams)
        phis[,2,,] <- plogis(X.eps %*% -extParams)

        forward(detParams, phis, psis) + 0.001*sqrt(sum(params^2))
    }

    if(missing(starts)) starts <- rep(0,nP)
    fm <- optim(starts, nll, method=method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    names(ests) <- parm.names
    fmAIC <- 2 * fm$value + 2 * nP # + 2*nP*(nP + 1)/(M - nP - 1)

    psiParams <- ests[1:nSP]
    colParams <- ests[(nSP + 1) : (nSP + nGP)]
    extParams <- ests[(nSP + nGP + 1) : (nSP + nGP + nEP)]
    detParams <- ests[(nSP + nGP + nEP + 1) : nP]

    psi <- list(name = "Initial", short.name = "psi",
                            estimates = psiParams,
                            covMat = as.matrix(covMat[1:nSP,1:nSP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    col <- list(name = "Colonization", short.name = "col",
                            estimates = colParams,
                            covMat = as.matrix(covMat[(nSP + 1) :
                                                          (nSP + nGP), (nSP + 1) : (nSP + nGP)]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    ext <- list(name = "Extinction", short.name = "ext",
                            estimates = extParams,
                            covMat = as.matrix(covMat[(nSP + nGP + 1) :
                                                          (nSP + nGP + nEP),
                                                      (nSP + nGP + 1) :
                                                          (nSP + nGP + nEP)]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    det <- list(name = "Detection", short.name = "p",
                            estimates = detParams,
                            covMat = as.matrix(covMat[(nSP+nGP+nEP+1):nP,
                                                      (nSP+nGP+nEP+1):nP]),
                            invlink = "logistic",
                            invlinkGrad = "logistic.grad")

    estimates <- list(state = psi, col = col, ext = ext, det=det)

    efit <- list(fitType = "occMS",
                 call = match.call(),
                 lamformula = lamformula,
                 gamformula = gamformula,
                 epsformula = epsformula,
                 detformula = detformula,
                 data = data, sitesRemoved = designMats$removed.sites,
                 estimates = estimates,
                 AIC = fmAIC, opt = fm, negLogLike = fm$value,
                 nllFun = fm$nll)

    class(efit) <- c('efitMS','efit','list')
    return(efit)
}


## Helper functions

gDetVec2<- function(y, detVec, mp) {
    if(y==0) {
        detVec[2]<- detVec[2] * (1/(1 + exp(mp)))
    } else {
        detVec[1]<- 0
        detVec[2]<- detVec[2] * (exp(mp)/(1 + exp(mp)))
    }
    return(detVec)
}

gSingleDetVec<- function(y, mp) {
    K<- 2
    detVec<- rep(1, K)
    detVec<- gDetVec2(y, detVec, mp)
    return(detVec)
}


gDetVecs<- function(y, mp, Ji, tin) {
    ndim<- dim(mp)
    nDMP<- ndim[1]
    J<- ndim[2]
    nY<- ndim[3]
    M<- ndim[4]
    K<- 2
    t<- tin - 1
    detVec<- rep(0, K*M)
    dind<- 1
    for(i in 0:(M-1)) {
        detVec[dind:(dind+1)]<- 1

        for(j in 0:(Ji[i+1]-1)) {
            yind<- i+t*M + j*M*nY + 1
            mpind<- j*nDMP + t*nDMP*J + i*nDMP*J*nY + 1
            if((j >= 0) && !is.na(y[yind])) {
                detVec[dind:(dind+1)]<- gDetVec2(y[yind], detVec[dind:(dind+1)], mp[mpind])
            }
        }
        dind<- dind + K
    }
    return(detVec)
}

