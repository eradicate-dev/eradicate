
#' remPoisM
#'
#' @name remPoisM
#'
#' @description
#' \code{remPoisM} fits the combined removal model + royle-nichols
#' occupancy/abundance model to two sources of monitoring data. The first
#' (main) source of data is removal type data consisting of a number of primary
#' periods where individuals removed are recorded for each site. The second
#' source of data consists of standard occupancy data which is collected at the same
#' set of sites (or a proportion thereof) using a different method (e.g. cameras).
#' The model then estimates abundance in light of detections of individuals arising from both
#' sources of data.
#'
#' @usage remPoisM(lamformula, detformula, data, K, starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the removal detection component.  Only
#'  site-level covariates are allowed for the removal detection component.
#'  This differs from the similar model in \code{unmarked}. Currently
#'  an intercept-only model is assumed for the detection component of the occupancy data.
#' @param data A \code{eFrameRM} object containing the response (counts)
#'  and site-level covariates. see \code{\link{eFrameRM}} for how to format
#'  the required data.
#' @param K Integer upper index of integration for abundance. This should be
#'  set high enough so that it does not affect the parameter estimates. Note
#'  that computation time will increase with K.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  ## uses san_nic_rem
#'  emf <- eFrameRM(y=rem, y1=y1, cells=cells, Z=nights)
#'  mod <- remPoisM(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remPoisM <- function(lamformula, detformula, data, K=25, starts, method = "BFGS",
    se = TRUE, ...)
{
    if(!is(data, "eFrameRM"))
		    stop("Data is not a data frame or eFrameRM.")
    D <- getDesign(data, lamformula, detformula)
    y <- D$y
    ym <- D$ym
    X <- D$X
    Xm<- D$Xm
    V <- D$V;
    Vm<- D$Vm
    Z<- D$Z
    X.offset <- D$X.offset
    Xm.offset <- D$Xm.offset
    V.offset <- D$V.offset
    y.sites<- D$y.sites
    ym.sites<- D$ym.sites

    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if (is.null(Xm.offset)) {
      Xm.offset <- rep(0, nrow(Xm))
    }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    M <- nrow(y)
    R <- ncol(y)
    P <- nrow(ym)
    Q <- ncol(ym)

    rsites<- ym.sites %in% y.sites # sites with both removal and monitoring

    piFun <- data$piFun
    lamParms <- colnames(X)
    detParms <- c(colnames(V), colnames(Vm))

    nDP <- ncol(V)
    nAP <- ncol(X)
    nDP1 <- ncol(Vm)

    nP <- nDP + nAP +  nDP1

    n <- 0:K

    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    yvec <- as.numeric(y)
    navec <- is.na(yvec)
    ymvec <- as.numeric(ym)
    navecm <- is.na(ymvec)

    nll <- function(parms) {
        lambda <- exp(X %*% parms[1 : nAP] + X.offset)
        p <- plogis(V %*% parms[(nAP + 1) : (nDP + nAP)] + V.offset)
        p.matrix <- matrix(p, M, R, byrow = TRUE)
        pi <- do.call(piFun, list(p = p.matrix))
        logLikeR <- dpois(y, matrix(lambda, M, R) * pi, log = TRUE)
        logLikeR[navec] <- 0
        logLikeR<- sum(logLikeR)
        # add in extra monitoring in ym
        r.ij <- matrix(plogis(Vm %*% parms[(nDP + nAP + 1) : nP]), P, Q,
                       byrow = TRUE)
        r.pi <- do.call(piFun, list(p = r.ij))
        p.ij.list <- lapply(n, function(k) 1 - (1 - r.ij)^k)
        cp.ij.list <- lapply(p.ij.list, function(pmat) pmat^ym * (1-pmat)^(Z-ym))
        cp.ij.list <- lapply(cp.ij.list, function(cpmat) {
          cpmat[navecm] <- 1
          cpmat
        })
        ## compute P(N = n | lambda * pi) along i
        lambda.m <- exp(Xm %*% parms[1 : nAP] + Xm.offset)
        lambda.i <- matrix(lambda.m, P, Q)
        lambda.i[rsites,] <- lambda.i[rsites,] * r.pi[rsites,] # marginal lambda for removal sites
        lambda.in <- sapply(n, function(x) dpois(x, lambda.i))
        cp.mat <- sapply(cp.ij.list, as.vector)
        cp.in <- rowSums(cp.mat * lambda.in)
        loglikeM<- sum(log(cp.in))
        -(logLikeR + loglikeM)
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

    typeNames<- c("state","det")

    stateEstimates <- list(name = "Abundance",
                                   short.name = "lambda",
                                   estimates = ests[1:nAP],
                                   covMat = as.matrix(covMat[1:nAP,1:nAP]),
                                   invlink = "exp",
                                   invlinkGrad = "exp")

    detEstimates <- list(name = "Detection", short.name = "p",
                                 estimates = ests[(nAP + 1) : nP],
                                 covMat = as.matrix(covMat[(nAP + 1):nP, (nAP + 1):nP]),
                                 invlink = "logistic",
                                 invlinkGrad = "logistic.grad")

    efit <- list(fitType = "Multinomial Removal + Monitoring",
        call = match.call(), types=typeNames, lamformula = lamformula, detformula=detformula,
        state=stateEstimates,det=detEstimates, sitesRemoved = D$removed.sites,
        AIC = fmAIC, opt = opt, negLogLike = fm$value, nllFun = nll, data = data)
    class(efit) <- c('efitR','efit','list')
    return(efit)
}
