#' remMNO
#'
#' @name remMNO
#'
#' @description
#' \code{remMNO} fits the multinomial, open population removal model to data collected from
#' repeated removal episodes from M sites over T primary periods with each primary consisting of J
#' secondary periods. This is a port of the similar function in \code{unmarked}.
#'
#' @usage remMNO(lamformula, gamformula, omformula, detformula, data, mixture = c("P", "NB","ZIP"),
#'                K, dynamics = c("constant", "autoreg", "notrend", "trend", "ricker", "gompertz"),
#'                fix = c("none", "gamma", "omega"), immigration=FALSE, iotaformula = ~1,
#'                starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param gamformula formula for availability
#' @param omformula formula for availability
#' @param detformula formula for the removal detection component.  Only
#'  site-level covariates are allowed for the removal detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrameMNO} object containing the removal data for each primary
#' period and site-level covariates. see \code{\link{eFrameMNO}} for how to format
#'  the required data.
#' @param mixture for abundance, either Poisson 'P', negative binomial 'NB' or zero-inflated
#' poisson 'ZIP'.
#' @param K upper bound for superpopulation abundance
#' @param dynamics Character string describing the type of population dynamics. "constant" indicates that there
#' is no relationship between omega and gamma. "autoreg" is an auto-regressive model in which recruitment is
#' modeled as gamma*N[i,t-1]. "notrend" models gamma as lambda*(1-omega) such that there is no temporal trend.
#' "trend" is a model for exponential growth, N[i,t] = N[i,t-1]*gamma, where gamma in this case is finite rate of
#' increase (normally referred to as lambda). "ricker" and "gompertz" are models for density-dependent population
#' growth. "ricker" is the Ricker-logistic model, N[i,t] = N[i,t-1]*exp(gamma*(1-N[i,t-1]/omega)), where gamma is
#' the maximum instantaneous population growth rate (normally referred to as r) and omega is the equilibrium
#' abundance (normally referred to as K). "gompertz" is a modified version of the Gompertz-logistic model, N[i,t]
#' = N[i,t-1]*exp(gamma*(1-log(N[i,t-1]+1)/log(omega+1))), where the interpretations of gamma and omega are
#' similar to in the Ricker model.
#' @param fix If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.
#' @param immigration Logical specifying whether immigration is included in the model
#' @param iotaformula formula for the number of immigrants per site, per time step.
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  emf <- eFrameMNO(y=rem, numPrimary=1)
#'  mod <- remMNO(~1, ~1, ~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'

remMNO <- function(lamformula, gamformula, omformula, detformula,
    data, mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend", "ricker", "gompertz"),
    fix=c("none", "gamma", "omega"), immigration=FALSE, iotaformula = ~1,
    starts, method="BFGS", se=TRUE, ...) {

  #Check data source
  if(!is(data, "eFrameMNO"))
    stop("Data is not of class eFrameMNO.")

  #Check state model arguments
  mixture <- match.arg(mixture)
  dynamics <- match.arg(dynamics)

  if((identical(dynamics, "constant") || identical(dynamics, "notrend")) & immigration)
    stop("You can not include immigration in the constant or notrend models")

  if(identical(dynamics, "notrend") &
   !identical(lamformula, omformula))
    stop("lamformula and omformula must be identical for notrend model")

  fix <- match.arg(fix)

  D <- getDesign(data, lamformula, gamformula, omformula, detformula, iotaformula)
  y <- D$y

  M <- nrow(y)
  T <- data$numPrimary
  J <- ncol(y) / T
  piFun <- data$piFun

  y <- array(y, c(M, J, T))
  yt <- apply(y, c(1,3), function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
  })

  ytna <- apply(is.na(y), c(1,3), all)
  ytna <- matrix(ytna, nrow=M)
  ytna[] <- as.integer(ytna)

  first <- apply(!ytna, 1, function(x) min(which(x)))
  last  <- apply(!ytna, 1, function(x) max(which(x)))
  first1 <- which(first==1)[1]

  Xlam.offset <- D$Xlam.offset
  Xgam.offset <- D$Xgam.offset
  Xom.offset <- D$Xom.offset
  Xp.offset <- D$Xp.offset
  Xiota.offset <- D$Xiota.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
  if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
  if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
  if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*J)
  if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

  #K stuff
  if(missing(K)) {
    K <- max(y, na.rm=T) + 20
    warning("K was not specified and was set to ", K, ".")
  }
  if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
  k <- 0:K
  lk <- length(k)
  #Some k-related indices to avoid repeated calculations in likelihood
  lfac.k <- lgamma(k+1)
  kmyt <- array(0, c(lk, T, M))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- array(NA, c(M, T, lk)) #Indicator if given k is possible given y
  for(i in 1:M) {
    for(t in 1:T) {
      fin[i,t,] <- k - yt[i,t] >= 0
      if(sum(ytna[i,t])==0) {
        kmyt[,t,i] <- k - yt[i,t]
        lfac.kmyt[i,t, ] <- lgamma(kmyt[,t,i] + 1)
      }
    }
  }

  lamParms <- colnames(D$Xlam)
  gamParms <- colnames(D$Xgam)
  omParms <- colnames(D$Xom)
  detParms <- colnames(D$Xp)
  nAP <- ncol(D$Xlam)
  nGP <- ncol(D$Xgam)
  nOP <- ncol(D$Xom)
  nDP <-  ncol(D$Xp)

  nIP <- ifelse(immigration, ncol(D$Xiota), 0)
  iotaParms <- character(0)
  if(immigration) iotaParms <- colnames(D$Xiota)

  if(identical(fix, "gamma")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nGP > 1){
        stop("gamma covariates not allowed when fix==gamma")
    }else {
        nGP <- 0
        gamParms <- character(0)
    }
  } else if(identical(dynamics, "notrend")) {
    if(nGP > 1){
        stop("gamma covariates not allowed when dyamics==notrend")
    } else {
        nGP <- 0
        gamParms <- character(0)
    }
  }

  if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nOP > 1)
        stop("omega covariates not allowed when fix==omega")
    else {
        nOP <- 0
        omParms <- character(0)
    }
  } else if(identical(dynamics, "trend")) {
    if(nOP > 1)
        stop("omega covariates not allowed when dynamics='trend'")
    else {
        nOP <- 0
        omParms <- character(0)
    }
  }

  nP <- nAP + nGP + nOP + nDP + nIP + (mixture!="P")
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))

  nbParm <- character(0)
  if(identical(mixture,"NB"))
    nbParm <- "alpha"
  else if(identical(mixture, "ZIP"))
    nbParm <- "psi"

  paramNames <- c(lamParms, gamParms, omParms, detParms,
                 iotaParms, nbParm)

  #Create indices, all possible combinations of survivors and recruits,
  #finding all unique likelihood transitions
  I <- cbind(rep(k, times=lk), rep(k, each=lk))
  I1 <- I[I[,1] <= I[,2],]
  lik_trans <- .Call("get_lik_trans", I, I1, PACKAGE="unmarked")

  beta_ind <- matrix(NA, 6, 2)
  beta_ind[1,] <- c(1, nAP) #Abundance
  beta_ind[2,] <- c(1, nGP) + nAP #Gamma
  beta_ind[3,] <- c(1, nOP) + nAP + nGP #Omega
  beta_ind[4,] <- c(1, nDP) + nAP + nGP + nOP #Detection
  beta_ind[5,] <- c(1, nIP) + nAP + nGP + nOP + nDP #Iota
  beta_ind[6,] <- c(1, 1) + nAP + nGP + nOP + nDP + nIP #2nd abun param

  #Adjustments to objects to facilitate use in c++
  fin <- fin*1 #convert to numeric
  yperm <- aperm(y, c(1,3,2))
  yna <- is.na(yperm)*1

  nll <- function(parms) {
    .Call("nll_multmixOpen",
          yperm, yt,
          D$Xlam, D$Xgam, D$Xom, D$Xp, D$Xiota,
          parms, beta_ind - 1,
          Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
          ytna, yna,
          lk, mixture, first - 1, last - 1, first1 - 1, M, T, J,
          D$delta, dynamics, fix, D$go.dims, immigration,
          I, I1, lik_trans$Ib, lik_trans$Ip,
          piFun, lfac.k, kmyt, lfac.kmyt, fin,
          PACKAGE = "unmarked")
  }

  if(missing(starts)){
    starts <- rep(0, nP)
  }

  fm <- optim(starts, nll, method=method, hessian=se, ...)
  ests <- fm$par
  names(ests) <- paramNames
  covMat <- invertHessian(fm, nP, se)
  fmAIC <- 2*fm$value + 2*nP

  lamEstimates <- list(name = "Abundance", short.name = "lam",
                    estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]),
                    invlink = "exp", invlinkGrad = "exp")
  estimates <- list(lambda=lamEstimates)

  gamName <- switch(dynamics, constant = "gamConst", autoreg = "gamAR",
                              notrend = "", trend = "gamTrend",
                              ricker="gamRicker", gompertz = "gamGomp")
  if(!(identical(fix, "gamma") | identical(dynamics, "notrend"))){
    estimates$gamma <- list(name =
        ifelse(identical(dynamics, "constant") | identical(dynamics, "autoreg"),
        "Recruitment", "Growth Rate"), short.name = gamName,
        estimates = ests[(nAP+1) : (nAP+nGP)], covMat = as.matrix(covMat[(nAP+1) :
                           (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
  }

  if(!(identical(fix, "omega") | identical(dynamics, "trend"))) {
    if(identical(dynamics, "constant") | identical(dynamics, "autoreg") |
       identical(dynamics, "notrend")){
        estimates$omega <- list(name="Apparent Survival",
          short.name = "omega", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
          covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                                    (nAP+nGP+1) : (nAP+nGP+nOP)]),
          invlink = "logistic", invlinkGrad = "logistic.grad")
    } else if(identical(dynamics, "ricker")){
        estimates$omega <- list(name="Carrying Capacity",
          short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
          covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                            (nAP+nGP+1) : (nAP+nGP+nOP)]),
          invlink = "exp", invlinkGrad = "exp")
    } else{
      estimates$omega <- list(name="Carrying Capacity",
        short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                                  (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "exp", invlinkGrad = "exp")
    }
  }

  estimates$det <- list(name = "Detection", short.name = "p",
      estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
                        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
      invlink = "logistic", invlinkGrad = "logistic.grad")

  if(immigration) {
    estimates$iota <- list(name="Immigration", short.name = "iota",
      estimates = ests[(nAP+nGP+nOP+nDP+1) :(nAP+nGP+nOP+nDP+nIP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP),
                                (nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP)]),
      invlink = "exp", invlinkGrad = "exp")
  }

  if(identical(mixture, "NB")) {
    estimates$alpha <- list(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")
  }
  if(identical(mixture, "ZIP")) {
    estimates$psi <- list(name = "Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
  }

  efit <- list(fitType = "multmixOpen",
      call = match.call(), lamformula = lamformula, detformula=detformula,
      gamformula=gamformula, omformula=omformula, iotaformula=iotaformula,
      data = data, sitesRemoved=D$removed.sites, estimates = estimates, AIC = fmAIC,
      opt = fm, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
      dynamics = dynamics, fix = fix, immigration=immigration)

  class(efit) <- c('efitMNO','efit','list')

  return(efit)
}
