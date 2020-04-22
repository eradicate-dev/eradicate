
#' remMN
#'
#' @name remMN
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
#'  This differs from the similar model in \code{unmarked}.
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
#'  rem<- san_nic_rem$rem
#'  emf <- eFrameR(y=rem)
#'  mod <- remPois(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remMN <- function(lamformula, detformula, data, mixture=c("P","NB","ZIP"), starts,
                    method = "BFGS", se = TRUE, ...) {

    if(!is(data, "eFrameR"))
		    stop("Data is not a eFrameR.")
    mixture <- match.arg(mixture)
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
    M <- nrow(y)
    piFun <- data$piFun

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    nP <- nDP + nAP + (mixture != "P")
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    namat <- is.na(y)
    ylogfact<- lgamma(y + 1)

    if(identical(mixture, "P"))  {
        nll <- function(parms) {
            lambda <- exp(X %*% parms[1 : nAP] + X.offset)
            lambda.mat <- matrix(lambda, M, J)
            p <- plogis(V %*% parms[(nAP + 1) : nP] + V.offset)
            p.matrix <- matrix(p, M, J, byrow = TRUE)
            pi <- do.call(piFun, list(p = p.matrix))
            log.lam.pi<- sweep(log(pi),1,log(lambda),'+')
            lam.pi<- sweep(pi, 1, lambda, '*')
            logLike <- y * (log(lambda.mat) + log(pi)) - (lambda.mat * pi) - ylogfact
            logLike[namat] <- 0
            -sum(logLike)
         }
    }
    else if(identical(mixture, "NB")){
        nll <- function(parms) {
            lambda <- exp(X %*% parms[1 : nAP] + X.offset)
            lambda.mat <- matrix(lambda, M, J)
            p <- plogis(V %*% parms[(nAP + 1) : nP] + V.offset)
            alpha <- exp(parms[nP])
            p.matrix <- matrix(p, M, J, byrow = TRUE)
            pi <- do.call(piFun, list(p = p.matrix))
            yrow<- rowSums(y)
            yrlogfact<- rowSums(ylogfact)
            ptot<- rowSums(pi)
            logLike<- rep(NA, M)
            for(i in 1:M){
                ll1<- lgamma(alpha+yrow[i]) - (lgamma(alpha) + yrlogfact[i])
                notna<- !namat[i,]
                ll2<- y[i,notna] * (log(lambda[i]) + log(pi[i,notna]) - log(alpha + lambda[i]*ptot[i]))
                ll3<- alpha * (log(alpha) - log(alpha + lambda[i] * ptot[i]))
                logLike[i]<- ll1 + sum(ll2) + ll3
            }
            -sum(logLike)
        }
    }
    if(missing(starts))
        starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)

    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    fmAIC <- 2 * fm$value + 2 * nP

    if(identical(mixture,"NB"))
        names(ests)<- c(lamParms,detParms,"alpha")
    else if(identical(mixture,"ZIP"))
        names(ests)<- c(lamParms,detParms,"phi")
    else
        names(ests)<- c(lamParms, detParms)

    typeNames<- c("state","det")

    stateEstimates <- list(name = "Abundance",
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

    if(identical(mixture,"NB")) {
        typeNames<- c(typeNames,"disp")
        dispEstimates <- list(name="Dispersion", short.name="disp",
                              estimates = ests[nP],
                              covMat = as.matrix(covMat[nP,nP]),
                              invlink = "exp", invlinkGrad = "exp")
    }
    else dispEstimates<- NULL

    if(identical(mixture,"ZIP")) {
        typeNames<- c(typeNames,"zeroinf")
        zeroinfEstimates <- list(name="Zero inflation", short.name="phi",
                              estimates = ests[nP],
                              covMat = as.matrix(covMat[nP,nP]),
                              invlink = "logistic", invlinkGrad = "logistic.grad")
    }
    else zeroinfEstimates<- NULL

    efit <- list(fitType = "Multinomial removal",
        call = match.call(), types=typeNames,lamformula = lamformula,detformula=detformula,
        state=stateEstimates,det=detEstimates, disp=dispEstimates,zeroinf=zeroinfEstimates,
        sitesRemoved = designMats$removed.sites,
        AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll, data = data)
    class(efit) <- c('efitR','efit','list')
    return(efit)
}
