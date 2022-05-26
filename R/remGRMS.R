
#' remGRMS
#'
#' @name remGRMS
#'
#' @description
#' \code{remGRMS} fits the generalized removal model to 'stacked' data collected from
#' repeated removal episodes from M sites over T primary periods with each primary consisting of J
#' secondary periods. The model also facilitates the analysis of index (count) data collected
#' in conjunction with the removal data to make joint inference on abundance.  Currently
#' supported models include the Poisson and Negative Binomial
#'
#' @usage remGRMS(lamformula, detformula, mdetformula, data, mixture = c("P", "NB"), K,
#'                   starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the removal detection component.
#' @param mdetformula formula for the index detection component.
#' @param data A \code{eFrameGRM} object containing the response (counts, index)
#'  and site-level covariates. see \code{\link{eFrameGRM}} for how to format
#'  the required data.
#' @param model for abundance, either Poisson 'P' or negative binomial 'NB'
#' @param K upper bound for superpopulation abundance
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efitGRMS} model object.
#'
#' @examples
#'  rem<- san_nic_open$removal
#'  ym<- san_nic_open$index
#'  emf <- eFrameGRMS(rem, ym)
#'  mod <- remGRMS(~1, ~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remGRMS <- function(lamformula, detformula, mdetformula, data, mixture=c('P', 'NB'),
                  K, starts, method = "BFGS", se = TRUE, ...)
{
  if(!is(data, "eFrameGRMS"))
    stop("Data is not a eFrameGRMS.")

  mixture <- match.arg(mixture)

  D <- getDesign(data, lamformula, detformula, mdetformula)

  Xlam <- D$Xlam
  Xdet <- D$Xdet
  Xdetm <- D$Xdetm
  y <- D$y  # MxJT
  ym <- D$ym

  Xlam.offset <- D$Xlam.offset
  Xdet.offset <- D$Xdet.offset
  Xdetm.offset <- D$Xdetm.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))
  if(is.null(Xdetm.offset)) Xdet.offset <- rep(0, nrow(Xdet))

  if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 100
  k <- 0:K
  lk <- length(k)
  M <- nrow(y)
  R <- ncol(y)

  cumy<- y
  for(i in 1:M){
      if(all(is.na(y[i,])))
        cumy[i,]<- NA
      else {
        for(r in 1:R){
          cumy[i,r]<- sum(y[i,1:r], na.rm=TRUE) - y[i,r]
        }
    }
  }

  yt<- as.matrix(cumy[,R])

  piFun <- data$piFun

  lamParms <- colnames(Xlam)
  detParms <- colnames(Xdet)
  detmParms<- colnames(Xdetm)

  nLP <- ncol(Xlam)
  nDP <- ncol(Xdet)
  nDPM<- ncol(Xdetm)

  nP <- nLP + nDP + nDPM + (mixture=='NB')
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))


  lfac.k <- lgamma(k+1)
  kmyt <- matrix(NA, M, lk)
  lfac.kmyt <- matrix(0, M, lk)
  fin <- matrix(NA, M, lk)
  naflag <- matrix(NA, M, R)
  # remaining in each i,r, for every k
  krem<- array(NA, c(M ,R, lk))
  namflag<- matrix(NA, M, R)

  for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
      naflag[i,] <- is.na(y[i,])
      namflag[i,]<- is.na(ym[i,])
      if(!all(naflag[i,])) {
        kmyt[i,] <- k - yt[i]
        lfac.kmyt[i, fin[i,]] <- lgamma(kmyt[i, fin[i,]] + 1)
        for(r in 1:R) {
          krem[i,r,]<- ifelse(k - cumy[i,r] < 0, NA, k - cumy[i,r])
        }
      }
    }

  nll <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    p <- plogis(Xdet %*% pars[(nLP+1):(nLP+nDP)] + Xdet.offset)
    pm <- exp(Xdetm %*% pars[(nLP+nDP+1):(nLP+nDP+nDPM)])

    p <- matrix(p, M, R, byrow=TRUE)
    pm<- matrix(pm, M, R, byrow=TRUE)

    cp <- matrix(as.numeric(NA), M, R+1)

    pi <- do.call(piFun, list(p = p))
    cp[,1:R]<- pi
    cp[,1:R][is.na(y)]<- NA
    cp[,R+1] <- 1 - apply(cp[,1:R,drop=FALSE], 1, sum, na.rm=TRUE)

    switch(mixture,
           P = ff <- sapply(k, function(x) dpois(x, lambda)),
           NB = ff <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[nP]))))
    g <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
      na <- naflag[i,]
      if(!all(na))
      A<- lfac.k - lfac.kmyt[i, ] +
            sum(y[i, !na] * log(cp[i, which(!na)])) +
            kmyt[i, ] * log(cp[i, R+1])
      g[i,] <- exp(A)
    }
    ff[!fin] <- g[!fin] <- 0
    ll <- rowSums(ff*g)
    # Add index-removal LL for ym
    gm <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
      AM <- matrix(0, lk,  R)
      na <- naflag[i,]
      nam <- namflag[i,]
      isr<- which(!as.logical(na + nam))
      if(!all(na))
         for(r in isr) {
            AM[, r] <- dpois(ym[i,r], krem[i,r,]*pm[i,r], log=TRUE)
          }

      gm[i,] <- exp(apply(AM, 1, sum))
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
    names(ests)<- c(lamParms,detParms,detmParms,"alpha")
  else
    names(ests)<- c(lamParms,detParms,detmParms)


  stateEstimates <- list(name = "Abundance", short.name = "lambda",
                         estimates = ests[1:nLP],
                         covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
                         invlinkGrad = "exp")

  detEstimates <- list(name = "Detection", short.name = "p",
                       estimates = ests[(nLP+1):(nLP+nDP)],
                       covMat = as.matrix(covMat[(nLP+1):(nLP+nDP),(nLP+1):(nLP+nDP)]),
                       invlink = "logistic",
                       invlinkGrad = "logictic.grad")

  mdetEstimates <- list(name = "mDetection", short.name = "pm",
                       estimates = ests[(nLP+nDP+1):(nLP+nDP+nDPM)],
                       covMat = as.matrix(covMat[(nLP+nDP+1):(nLP+nDP+nDPM),
                                                 (nLP+nDP+1):(nLP+nDP+nDPM)]),
                       invlink = "exp",
                       invlinkGrad = "exp")

  estimates<- list(state=stateEstimates, det=detEstimates, detm=mdetEstimates)

  if(identical(mixture,"NB")) {
    dispEstimates <- list(name = "Dispersion", short.name = "disp",
                          estimates = ests[nP],
                          covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
                          invlinkGrad = "exp")
    estimates$disp<- dispEstimates
  }

  efit <- list(fitType = "generalised removal",
               call = match.call(), lamformula = lamformula,
               detformula=detformula, mdetformula=mdetformula,
                estimates=estimates, sitesRemoved = D$removed.sites,
               AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll,
               mixture=mixture, K=K, data = data)
  class(efit) <- c('efitGRMS','efitR','efit','list')

  return(efit)
}




