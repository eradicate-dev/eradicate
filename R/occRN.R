
#' occuRN
#'
#' @name occRN
#'
#' @description
#' \code{occuRN} fits the occupancy/abundance model of Royle & Nichols (2003).
#'
#' @usage occuRN(lamformula, detformula, data, K, starts, method="BFGS", se=TRUE, ...)
#'
#' @param lamformula formula for the latent abundance component.
#' @param detformula formula for the detection component.  Only
#'  site-level covariates are allowed for the detection component.
#'  This differs from the similar model in \code{unmarked}.
#' @param data A \code{eFrame} object containing the response (either (0/1) or counts)
#'  and site-level covariates. see \code{\link{eFrame}} for how to format
#'  the required data. count data will get trunctated to (0/1) by \code{occuRN()}.
#' @param K Integer upper index of integration for abundance. This should be
#'  set high enough so that it does not affect the parameter estimates. Note
#'  that computation time will increase with K
#' @param starts Initial values for parameters
#' @param method Optimisation method
#' @param se flag to return the standard error (hessian).
#'
#' @return a \code{efit} model object.
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  mod <- occuRN(~1, ~1, data=emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
occRN <- function(lamformula, detformula, data, K = 25, starts, method = "BFGS", se = TRUE, ...)
{
    if(!is(data, "eFrame")) stop("Data is not an eFrame")

    designMats <- getDesign(data, lamformula, detformula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))

  y <- truncateToBinary(y)

  J <- ncol(y)
  M <- nrow(y)

  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)

  nP <- nDP + nOP
  if(!missing(starts) && length(starts) != nP)
	   stop(paste("The number of starting values should be", nP))

  y.ji <- as.vector(y)
  navec <- is.na(y.ji)
  n <- 0:K

  nll <- function(parms, f = "Poisson") {

    ## compute individual level detection probabilities
    r.ij <- matrix(plogis(V %*% parms[(nOP + 1) : nP] + V.offset), M, J,
      byrow = TRUE)

    ## compute list of detection probabilities along N
    p.ij.list <- lapply(n, function(k) 1 - (1 - r.ij)^k)

    ## compute P(y_{ij} | N) (cell probabilities) along N
    cp.ij.list <- lapply(p.ij.list, function(pmat) pmat^y * (1-pmat)^(1-y))

    ## replace NA cell probabilities with 1.
    cp.ij.list <- lapply(cp.ij.list, function(cpmat) {
      cpmat[navec] <- 1
      cpmat
    })

    ## multiply across J to get P(y_i | N) along N
    cp.in <- sapply(cp.ij.list, rowProds)

    ## compute P(N = n | lambda_i) along i
    lambda.i <- exp(X %*% parms[1 : nOP] + X.offset)
    lambda.in <- sapply(n, function(x) dpois(x, lambda.i))

    ## integrate over P(y_i | N = n) * P(N = n | lambda_i) wrt n
    like.i <- rowSums(cp.in * lambda.in)

    -sum(log(like.i))
  }

	if(missing(starts)) starts <- rep(0, nP)
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
  names(ests)<- c(occParms, detParms)

  stateEstimates <- list(name = "Abundance",
      short.name = "lam",
      estimates = ests[1:nOP],
      covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = "exp",
      invlinkGrad = "exp")

  detEstimates <- list(name = "Detection", short.name = "p",
      estimates = ests[(nOP + 1) : nP],
      covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]),
      invlink = "logistic", invlinkGrad = "logistic.grad")

  estimates<- list(state=stateEstimates, det=detEstimates)

  efit <- list(fitType = "occuRN",call = match.call(),
               lamformula = lamformula, detformula=detformula, estimates=estimates,
               sitesRemoved = designMats$removed.sites, AIC = fmAIC, opt = opt,
               negLogLike = fm$value, nllFun = nll, data = data)
  class(efit) <- c('efit','list')
  return(efit)
}
