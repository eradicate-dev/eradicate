
#' occuMS
#'
#' @name occuMS
#'
#' @description
#' \code{occuMS} fits the dynamic multi-season occupancy model of MacKenzie et al. (2003).
#' This is a port of the \code{colext} function in \code{unmarked}.
#'
#' @usage occuMS(psiformula, gamformula, epsformula, detformula,
#' data, starts, method="BFGS", se=TRUE, ...)
#'
#' @param psiformula formula for the latent occupancy component.
#' @param gamformula formula for the latent colonisation component. If season
#' specific colonisation/extinction parameters are desired, then use the
#' reserved keyword \code{.season}, (i.e. \code{~.season}).
#' @param epsformula formula for the latent survival component.
#' @param detformula formula for the detection component.  Only
#'  site- or \code{.season}-level covariates are allowed for
#'  the detection component. This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrameMS} object containing the detection/non-detection
#' data (0/1) and site-level covariates. see \code{\link{eFrameMS}} for how to format
#'  the required data. count data will get trunctated to (0/1) by \code{occuMS()}.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  mod <- occuMS(~1, ~1, ~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @useDynLib eradicate, .registration=TRUE
#' @export
#'
#'
occuMS <- function(psiformula = ~ 1, gamformula = ~ 1, epsformula = ~ 1, detformula = ~ 1,
                   data, starts, method = "BFGS", se = TRUE, ...)
{

    if(!is(data, "eFrameMS"))
        stop("Data is not a eFrameMS.")

    K <- 1
    ## truncate at K
    data$y[data$y > K] <- K
    y <- getY(data)
    J <- data$obsPerSeason

    M <- nrow(y)
    nY <- ncol(y)/J
    n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))

    fc <- match.call()
    fc[[1]] <- as.name("occuMS.fit")
    fc$psiformula = psiformula
    fc$gamformula = gamformula
    fc$epsformula = epsformula
    fc$detformula = detformula
    fc$data <- as.name("data")
    fc$J <- as.name("J")
    fc$method <- as.name("method")
    fc$getHessian <- as.name("se")
    fc$se <- NULL
    if(missing(starts)) {
        fc$starts <- NULL
    } else {
        fc$starts <- eval(starts)
    }

    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(fc)))
        for (a in names(extras)[existing]) fc[[a]] <- extras[[a]]
        if (any(!existing)) {
            fc <- as.call(c(as.list(fc), extras[!existing]))
        }
    }

    fm <- eval(fc)

    fm$n.det <- n.det
    opt <- fm$opt
    nP <- fm$nP; M <- fm$M; nDP <- fm$nDP; nGP <- fm$nGP
    nEP <- fm$nEP; nSP <- fm$nSP

    covMat <- invertHessian(opt, nP, se)
    ests <- opt$par
    names(ests) <- fm$mle$names
    fmAIC <- 2 * opt$value + 2 * nP # + 2*nP*(nP + 1)/(M - nP - 1)

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

    estimates <- list(psi = psi, col = col, ext = ext, det=det)

    efit <- list(fitType = "occuMS",
                 call = match.call(),
                 psiformula = psiformula,
                 gamformula = gamformula,
                 epsformula = epsformula,
                 detformula = detformula,
                 data = data, sitesRemoved = fm$designMats$removed.sites,
                 estimates = estimates,
                 AIC = fmAIC, opt = opt, negLogLike = opt$value,
                 nllFun = fm$nll,
                 projected = fm$projected,
                 projected.mean = fm$projected.mean,
                 smoothed = fm$smoothed, smoothed.mean = fm$smoothed.mean)
    class(efit) <- c('efitMS','efit','list')
    return(efit)
}


occuMS.fit <- function(psiformula, gamformula, epsformula, detformula, data, J,
                       starts=NULL, method, getHessian = TRUE, wts, ...) {
    K <- 1
    designMats <- getDesign(data, psiformula, gamformula, epsformula, detformula)
    V.itjk <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W

    detParms <- colnames(V.itjk)
    gamParms <- colnames(X.it.gam)
    epsParms <- colnames(X.it.eps)
    psiParms <- colnames(W.i)

    y <- designMats$y
    M <- nrow(y)
    nY <- ncol(y)/J
    if(missing(wts)) wts <- rep(1, M)

    ## remove final year from X.it
    X.it.gam <- as.matrix(X.it.gam[-seq(nY,M*nY,by=nY),])
    X.it.eps <- as.matrix(X.it.eps[-seq(nY,M*nY,by=nY),])

    nDP <- length(detParms)
    nGP <- length(gamParms)
    nEP <- length(epsParms)
    nSP <- length(psiParms)
    nDMP <-  1

    nP <- nDP + nSP + nGP + nEP  # total number of parameters

    y.itj <- as.numeric(t(y))

    ## replace NA's with 99 before passing to C++
    ## TODO: need better missing data passing mechanism (maybe NaN of Inf?)
    y.itj[is.na(y.itj)] <- 99
    V.itjk[is.na(V.itjk)] <- 9999
    # get ragged array indices
    y.it <- matrix(t(y), nY*M, J, byrow = TRUE)
    J.it <- rowSums(!is.na(y.it))

    V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
    V.arr <- aperm(V.arr, c(2,1,5,4,3))

    y.arr <- array(y.itj, c(J, nY, M))
    y.arr <- aperm(y.arr, c(3:1))
    storage.mode(J.it) <- storage.mode(y.arr) <- storage.mode(K) <- "integer"

    alpha <- array(NA, c(K + 1, nY, M))

    forward <- function(detParms, phis, psis, storeAlpha = FALSE) {

        negloglike <- 0
        psiSite <- matrix(c(1-psis,psis), K + 1, M, byrow = TRUE)

        mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
        for(t in 1:nY) {
            storage.mode(t) <- "integer"
            detVecs <- .Call("getDetVecs", y.arr, mp,
                             J.it[seq(from = t, to = length(J.it)-nY+t,
                                      by=nY)], t, K,
                             PACKAGE = "eradicate")
            psiSite <- psiSite * detVecs
            if(storeAlpha) alpha[,t,] <<- psiSite[,]
            if(t < nY) {
                for(i in 1:M) {
                    psiSite[,i] <- phis[,,t,i] %*% psiSite[,i]
                }
            } else {
                negloglike <- negloglike - sum(wts*log(colSums(psiSite)))
            }
        }
        negloglike
    }

    backward <- function(detParams, phis) {
        beta <- array(NA, c(K + 1, nY, M))
        for (i in 1:M) {
            backP <- rep(1, K + 1)
            for (t in nY:1) {

                beta[, t, i] <- backP

                detVec <- rep(1, K + 1)
                for (j in 1:J) {
                    if(y.arr[i,t,j] != 99) {
                        mp <- V.arr[,,i,t,j] %*% detParams
                        detVecObs <- .Call("getSingleDetVec",
                                           y.arr[i,t,j], mp, K,
                                           PACKAGE = "eradicate")
                        detVec <- detVec * detVecObs
                    }
                }
                if (t > 1)
                    backP <- t(phis[,,t-1,i]) %*% (detVec * backP)
            }
        }
        return(beta)
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

    if(is.null(starts)) starts <- rep(0,nP)
    fm <- optim(starts, nll, method=method, hessian = getHessian, ...)
    mle <- fm$par

    psis <- plogis(W.i %*% mle[1:nSP])
    colParams <- mle[(nSP + 1) : (nSP + nGP)]
    extParams <- mle[(nSP + nGP + 1) : (nSP + nGP + nEP)]
    detParams <- mle[(nSP + nGP + nEP + 1) : nP]

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
    projected.mean <- apply(projected, 1:2, mean)
    rownames(projected.mean) <- c("unoccupied","occupied")
    colnames(projected.mean) <- 1:nY

    ## smoothing
    forward(detParams, phis, psis, storeAlpha = TRUE)
    beta <- backward(detParams, phis)
    gamma <- array(NA, c(K + 1, nY, M))
    for(i in 1:M) {
    for(t in 1:nY) {
        gamma[,t,i] <- alpha[,t,i]*beta[,t,i] / sum(alpha[,t,i]*beta[,t,i])
    }}
    smoothed.mean <- apply(gamma, 1:2, mean)
    rownames(smoothed.mean) <- c("unoccupied","occupied")
    colnames(smoothed.mean) <- 1:nY

    parm.names <- c(psiParms, gamParms, epsParms, detParms)
    mle.df <- data.frame(names = parm.names, value = mle)
    rownames(mle.df) <- paste(c(rep("psi", nSP), rep("col", nGP),
                                rep("ext", nEP), rep("det", nDP)),
                              c(1:nSP,1:nGP,1:nEP, 1:nDP))

    list(mle = mle.df, opt=fm, nP = nP, M = M, nDP = nDP, nGP = nGP,
         nEP = nEP, nSP = nSP,
         nllFun = nll, designMats = designMats,
         projected = projected, projected.mean = projected.mean,
         smoothed = gamma,
         smoothed.mean = smoothed.mean)
}
