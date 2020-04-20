
#' remGRM
#'
#' @name remGRM
#'
#' @description
#' \code{remGRM} fits the generalized removal model to data collected from
#' repeated removal episodes from M sites over T primary periods with each primary consisting of J
#' secondary periods. The model also facilitates the analysis of index (count) data collected
#' in conjunction with the removal data to make joint inference on abundance.
#'
#' @usage remGRM(lamformula, phiformula, detformula, mdetformula, data, mixture = c("P", "NB"), K,
#'                   starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param phiformula formula for availability
#' @param detformula formula for the removal detection component.  Only
#'  site-level covariates are allowed for the removal detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param mdetformula formula for the index detection component.  Only
#'  site-level covariates are allowed for the index detection component.
#' @param data A \code{eFrameR} object containing the response (counts, index)
#'  and site-level covariates. see \code{\link{eFrameGRM}} for how to format
#'  the required data.
#' @param model for abundance, either Poisson 'P' or negative binomial 'NB'
#' @param K upper bound for superpopulation abundance
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efitGRM} model object.
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  ym<- san_nic_rem$ym
#'  emf <- eFrameGRM(rem, ym, numPrimary=1)
#'  mod <- remGRM(~1, ~1, ~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remGRM <- function(lamformula, phiformula, detformula, mdetformula, data, mixture=c('P', 'NB'),
                  K, starts, method = "BFGS", se = TRUE, ...)
{
  if(!is(data, "eFrameGRM"))
    stop("Data is not a eFrameGRM.")

  mixture <- match.arg(mixture)

  D <- getDesign(data, lamformula, phiformula, detformula, mdetformula)

  Xlam <- D$Xlam
  Xphi <- D$Xphi
  Xdet <- D$Xdet
  Xdetm <- D$Xdetm
  y <- D$y  # MxJT
  ym <- D$ym

  Xlam.offset <- D$Xlam.offset
  Xphi.offset <- D$Xphi.offset
  Xdet.offset <- D$Xdet.offset
  Xdetm.offset <- D$Xdetm.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
  if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))
  if(is.null(Xdetm.offset)) Xdet.offset <- rep(0, nrow(Xdet))

  if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 100
  k <- 0:K
  lk <- length(k)
  M <- nrow(y)
  T <- data$numPrimary
  R <- ncol(y) / T

  y <- array(y, c(M, R, T))
  y <- aperm(y, c(1,3,2))

  ym <- array(ym, c(M, R, T))
  ym <- aperm(ym, c(1,3,2))

  cumy<- y
  for(i in 1:M){
    for(t in 1:T) {
      if(all(is.na(y[i,t,])))
        cumy[i,t,]<- NA
      else {
        for(r in 1:R){
          cumy[i,t,r]<- sum(y[i,t,1:r], na.rm=TRUE) - y[i,t,r]
        }
      }
    }
  }

  yt<- as.matrix(cumy[,,R])

   piFun <- data$piFun

  lamParms <- colnames(Xlam)
  detParms <- colnames(Xdet)
  detmParms<- colnames(Xdetm)
  detmParms[1]<- "(mIntercept)"
  detParms<- c(detParms,detmParms)

  nLP <- ncol(Xlam)
  if(T==1) {
    nPP <- 0
    phiParms <- character(0)
  } else if(T>1) {
    nPP <- ncol(Xphi)
    phiParms <- colnames(Xphi)
  }
  nDP <- ncol(Xdet)
  nDPM<- ncol(Xdetm)

  nP <- nLP + nPP + nDP + nDPM + (mixture=='NB')
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))


  lfac.k <- lgamma(k+1)
  kmyt <- array(NA, c(M, T, lk))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- matrix(NA, M, lk)
  naflag <- array(NA, c(M, T, R))
  # remaining in each i,t,r, for every k
  krem<- array(NA, c(M ,T, R, lk))
  namflag<- array(NA, c(M, T, R))

  for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
    for(t in 1:T) {
      naflag[i,t,] <- is.na(y[i,t,])
      namflag[i,t,]<- is.na(ym[i,t,])
      if(!all(naflag[i,t,])) {
        kmyt[i,t,] <- k - yt[i,t]
        lfac.kmyt[i, t, fin[i,]] <- lgamma(kmyt[i, t, fin[i,]] + 1)
        for(r in 1:R) {
          krem[i,t,r,]<- ifelse(k - cumy[i,t,r] < 0, NA, k - cumy[i,t,r])
        }
      }
    }
  }

  nll <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    if(T==1)
      phi <- 1
    else if(T>1)
      phi <- drop(plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset))
      p <- cloglog(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)
      pm <- exp(Xdetm %*% pars[(nLP+nPP+nDP+1):(nLP+nPP+nDP+nDPM)])

    phi.mat <- matrix(phi, M, T, byrow=TRUE)
    phi <- as.numeric(phi.mat)

    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, R, T))
    p <- aperm(p, c(1,3,2))

    pm<- matrix(pm, nrow=M, byrow=TRUE)
    pm <- array(pm, c(M, R, T))
    pm <- aperm(pm, c(1,3,2))

    cp <- array(as.numeric(NA), c(M, T, R+1))

    for(t in 1:T) cp[,t,1:R] <- do.call(piFun, list(p[,t,]))
    cp[,,1:R] <- cp[,,1:R] * phi
    cp[,, 1:R][is.na(y)]<- NA
    cp[,,R+1] <- 1 - apply(cp[,,1:R,drop=FALSE], 1:2, sum, na.rm=TRUE)

    switch(mixture,
           P = ff <- sapply(k, function(x) dpois(x, lambda)),
           NB = ff <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[nP]))))
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
    ff[!fin] <- g[!fin] <- 0
    ll <- rowSums(ff*g)
    # Add index-removal LL for ym
    gm <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
      AM <- array(0, c(lk, T, R))
      for(t in 1:T) {
        na <- naflag[i,t,]
        nam <- namflag[i,t,]
        isr<- which(!as.logical(na + nam))
        if(!all(na))
          for(r in isr) {
            AM[, t, r] <- dpois(ym[i,t,r], krem[i,t,r,]*pm[i,t,r], log=TRUE)
          }
        }
      gm[i,] <- exp(apply(AM, 1:2, sum))
      gm[i,is.na(gm[i,])]<- 0  # account for -ve krem[i,t,r, ]
      if(all(nam)) ff[i,]<- 1
    }
    gm[!fin] <- 0
    llm <- rowSums(ff*gm)

     (-1)*(sum(log(ll))+sum(log(llm)))
  }

  if(missing(starts)) starts <- rep(0, nP)
  fm <- optim(starts, nll, method = method, hessian = se, ...)
  covMat <- invertHessian(fm, nP, se)
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP

  if(identical(mixture,"NB"))
    names(ests)<- c(lamParms,phiParms,detParms,"alpha")
  else
    names(ests)<- c(lamParms,phiParms,detParms)

  typeNames<- c("state","det")

  stateEstimates <- list(name = "Abundance", short.name = "lambda",
                         estimates = ests[1:nLP],
                         covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
                         invlinkGrad = "exp")

  if(identical(mixture,"NB")) {
    typeNames<- c(typeNames,"disp")

    dispEstimates <- list(name = "Dispersion", short.name = "disp",
                          estimates = ests[nP],
                          covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
                          invlinkGrad = "exp")
  }
  else dispEstimates<- NULL

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
                       estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP+nDPM)],
                       covMat = as.matrix(covMat[(nLP+nPP+1):(nLP+nPP+nDP+nDPM),(nLP+nPP+1):(nLP+nPP+nDP+nDPM)]),
                       invlink = "cloglog",
                       invlinkGrad = "cloglog.grad")

  efit <- list(fitType = "generalised removal",
               call = match.call(), types=typeNames,lamformula = lamformula, detformula=detformula,
               phiformula=phiformula, state=stateEstimates,det=detEstimates, avail=availEstimates,
               disp=dispEstimates, sitesRemoved = D$removed.sites,AIC = fmAIC, opt = fm,
               negLogLike = fm$value, nllFun = nll, mixture=mixture, K=K, data = data)
  class(efit) <- c('efitR','efit','list')

  return(efit)
}




