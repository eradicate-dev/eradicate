
#' remGR
#'
#' @name remGR
#'
#' @description
#' \code{remGR} fits the generalized removal model to data collected from
#' repeated removal episodes from M sites each consisting of J
#' secondary periods.
#'
#' @usage remGR(lamformula, detformula, data, mixture = c("P", "NB"), K,
#'                   starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the removal detection component.  Only
#'  site-level covariates are allowed for the removal detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrameR} object containing the response (counts)
#'  and site-level covariates. see \code{\link{eFrameR}} for how to format
#'  the required data.
#' @param model for abundance, either Poisson 'P' or negative binomial 'NB'
#' @param K upper bound for superpopulation abundance
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  emf <- eFrameGR(y=rem)
#'  mod <- remGR(~1, ~1, K=100, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remGR <- function(lamformula, detformula, data, mixture=c('P', 'NB'),
                  K, starts, method = "BFGS", se = TRUE, ...)
{
  if(!is(data, "eFrameGR"))
    stop("Data is not a eFrameGR.")

  mixture <- match.arg(mixture)

  D <- getDesign(data, lamformula, detformula)

  Xlam <- D$Xlam
  Xdet <- D$Xdet
  y <- D$y

  Xlam.offset <- D$Xlam.offset
  Xdet.offset <- D$Xdet.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

  if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 100
  k <- 0:K
  lk <- length(k)
  M <- nrow(y)
  R <- ncol(y)

  yt <- apply(y, 1, function(x) {
    if(all(is.na(x)))
      return(NA)
    else return(sum(x, na.rm=TRUE))
  })

  piFun <- data$piFun

  lamParms <- colnames(Xlam)
  detParms <- colnames(Xdet)

  nLP <- ncol(Xlam)
  nDP <- ncol(Xdet)
  nP <- nLP + nDP + (mixture=='NB')
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))


  lfac.k <- lgamma(k+1)
  kmyt <- matrix(NA, M, lk)
  lfac.kmyt <- matrix(0, M, lk)
  fin <- matrix(NA, M, lk)
  naflag <- matrix(NA, M, R)
  for(i in 1:M) {
    fin[i, ] <- k - max(yt[i], na.rm=TRUE) >= 0
    naflag[i,] <- is.na(y[i,])
    if(!all(naflag[i,])) {
        kmyt[i,] <- k - yt[i]
        lfac.kmyt[i, fin[i,]] <- lgamma(kmyt[i, fin[i,]] + 1)
      }
    }

  nll <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    p <- plogis(Xdet %*% pars[(nLP+1):(nLP+nDP)] + Xdet.offset)
    p <- matrix(p, nrow=M, byrow=TRUE)
    cp <- matrix(as.numeric(NA), M, R+1)

    pi <- do.call(piFun, list(p = p))
    cp[,1:R]<- pi
    cp[,1:R][is.na(y)]<- NA
    cp[,R+1] <- 1 - apply(cp[,1:R,drop=FALSE], 1, sum, na.rm=TRUE)

    switch(mixture,
           P = f <- sapply(k, function(x) dpois(x, lambda)),
           NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[nP]))))
    g <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
        na <- naflag[i,]
        if(!all(na))
          A <- lfac.k - lfac.kmyt[i, ] +
            sum(y[i, !na] * log(cp[i, which(!na)])) +
            kmyt[i,] * log(cp[i, R+1])

      g[i,] <- exp(A)
    }
    f[!fin] <- g[!fin] <- 0
    ll <- rowSums(f*g)
    -sum(log(ll))
  }

  if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll, method = method, hessian = se, ...)
  covMat <- invertHessian(fm, nP, se)
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP

  if(identical(mixture,"NB"))
    names(ests)<- c(lamParms,detParms,"alpha")
  else
    names(ests)<- c(lamParms,detParms)

  stateEstimates <- list(name = "Abundance", short.name = "lambda",
                         estimates = ests[1:nLP],
                         covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
                         invlinkGrad = "exp")

  detEstimates <- list(name = "Detection", short.name = "p",
                       estimates = ests[(nLP+1):(nLP+nDP)],
                       covMat = as.matrix(covMat[(nLP+1):(nLP+nDP),(nLP+1):(nLP+nDP)]),
                       invlink = "logistic",
                       invlinkGrad = "logistic.grad")

  estimates<- list(state=stateEstimates,det=detEstimates)

  if(identical(mixture,"NB")) {
    dispEstimates <- list(name = "Dispersion", short.name = "disp",
                           estimates = ests[nP],
                           covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
                           invlinkGrad = "exp")
    estimates$disp<- dispEstimates
  }

  efit <- list(fitType = "generalised removal",
               call = match.call(), lamformula = lamformula, detformula=detformula,
               estimates=estimates, sitesRemoved = D$removed.sites,
               AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll,
               mixture=mixture, K=K, data = data)
  class(efit) <- c('efitGR','efitR','efit','list')

  return(efit)
}




