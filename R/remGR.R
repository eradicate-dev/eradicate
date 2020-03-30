
#' remGR
#'
#' @name remGR
#'
#' @description
#' \code{remGR} fits the generalized removal model to data collected from
#' repeated removal episodes from M sites over T primary periods with each primary consisting of J
#' secondary periods.
#'
#' @usage remGR(lamformula, phiformula, detformula, data, mixture = c("P", "NB"), K,
#'                   starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param phiformula formula for availability
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
#'  emf <- eFrameR(y=rem)
#'  mod <- remPois(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remGR <- function(lamformula, phiformula, detformula, data, mixture=c('P', 'NB'),
                  K, starts, method = "BFGS", se = TRUE, ...)
{
  if(!is(data, "eFrameGR"))
    stop("Data is not a eFrameGR.")

  mixture <- match.arg(mixture)

  D <- getDesign(data, lamformula, phiformula, detformula)

  Xlam <- D$Xlam
  Xphi <- D$Xphi
  Xdet <- D$Xdet
  y <- D$y  # MxJT

  Xlam.offset <- D$Xlam.offset
  Xphi.offset <- D$Xphi.offset
  Xdet.offset <- D$Xdet.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
  if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

  if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 100
  k <- 0:K
  lk <- length(k)
  M <- nrow(y)
  T <- data$numPrimary
  R <- ncol(y) / T

  y <- array(y, c(M, R, T))
  y <- aperm(y, c(1,3,2))
  yt <- apply(y, 1:2, function(x) {
    if(all(is.na(x)))
      return(NA)
    else return(sum(x, na.rm=TRUE))
  })

  piFun <- data$piFun

  lamPars <- colnames(Xlam)
  detPars <- colnames(Xdet)
  nLP <- ncol(Xlam)
  if(T==1) {
    nPP <- 0
    phiPars <- character(0)
  } else if(T>1) {
    nPP <- ncol(Xphi)
    phiPars <- colnames(Xphi)
  }
  nDP <- ncol(Xdet)
  nP <- nLP + nPP + nDP + (mixture=='NB')
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))


  lfac.k <- lgamma(k+1)
  kmyt <- array(NA, c(M, T, lk))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- matrix(NA, M, lk)
  naflag <- array(NA, c(M, T, R))
  for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
    for(t in 1:T) {
      naflag[i,t,] <- is.na(y[i,t,])
      if(!all(naflag[i,t,])) {
        kmyt[i,t,] <- k - yt[i,t]
        lfac.kmyt[i, t, fin[i,]] <- lgamma(kmyt[i, t, fin[i,]] + 1)
      }
    }
  }

  nll <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    if(T==1)
      phi <- 1
    else if(T>1)
      phi <- drop(plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset))
      p <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)

    phi.mat <- matrix(phi, M, T, byrow=TRUE)
    phi <- as.numeric(phi.mat)

    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, R, T))
    p <- aperm(p, c(1,3,2))
    cp <- array(as.numeric(NA), c(M, T, R+1))

    for(t in 1:T) cp[,t,1:R] <- do.call(piFun, list(p[,t,]))
    cp[,,1:R] <- cp[,,1:R] * phi
    cp[,, 1:R][is.na(y)]<- NA
    cp[,,R+1] <- 1 - apply(cp[,,1:R,drop=FALSE], 1:2, sum, na.rm=TRUE)

    switch(mixture,
           P = f <- sapply(k, function(x) dpois(x, lambda)),
           NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[nP]))))
    g <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
      A <- matrix(0, lk, T)
      for(t in 1:T) {
        na <- naflag[i,t,]
        if(!all(na))
          A[, t] <- lfac.k - lfac.kmyt[i, t,] +
            sum(y[i, t, !na] * log(cp[i, t, which(!na)])) +
            kmyt[i, t,] * log(cp[i, t, R+1])
      }
      g[i,] <- exp(rowSums(A))
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

  typeNames<- c("state","det")

    stateEstimates <- list(name = "Abundance", short.name = "lambda",
                         estimates = ests[1:nLP],
                         covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
                         invlinkGrad = "exp")
  if(T>1) {
    typeNames<- c(typeNames,"avail")
    availEstimates <- list(name = "Availability",
                         short.name = "phi",
                         estimates = ests[(nLP+1):(nLP+nPP)],
                         covMat = as.matrix(covMat[(nLP+1):(nLP+nPP),(nLP+1):(nLP+nPP)]),
                         invlink = "logistic",
                         invlinkGrad = "logistic.grad")
  }
    else availEstimates<- NULL

    detEstimates <- list(name = "Detection", short.name = "p",
                       estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
                       covMat = as.matrix(covMat[(nLP+nPP+1):(nLP+nPP+nDP),(nLP+nPP+1):(nLP+nPP+nDP)]),
                       invlink = "logistic",
                       invlinkGrad = "logistic.grad")

  if(identical(mixture,"NB")) {
    stateEstimates$alpha <- list(name="Dispersion",
                                 short.name = "alpha", estimates = ests[nP],
                                 covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
                                 invlinkGrad = "exp")
  }

  efit <- list(fitType = "generalised removal",
               call = match.call(), types=typeNames,lamformula = lamformula, detformula=detformula,
               phiformula=phiformula, state=stateEstimates,det=detEstimates, avail=availEstimates,
               sitesRemoved = D$removed.sites,AIC = fmAIC, opt = fm,
               negLogLike = fm$value, nllFun = nll, mixture=mixture, K=K, data = data)
  class(efit) <- c('efitR','efit','list')

  return(efit)
}




