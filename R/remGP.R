#' remGP
#'
#' @name remGP
#'
#' @description
#' \code{remGP} fits the catch-effort model of Gould & Pollock (1997) to removal
#' data from a single site (or combined data from many sites).  At least 3 removal periods
#' are required to fit this model.  Additionally, the model also accepts index (count)
#' data collected in conjunction with the removal data.
#'
#' @usage remGP(catch, effort, index, starts, Nmax, se=TRUE, ...)
#'
#' @param catch the number of removals recorded for each period.
#' @param effort the overall effort expended in each period (i.e. trapnights)
#' @param index Index (i.e. count) data collected in conjunction with the removal
#' data.  The index data is assumed to be collected just before each removal period and
#' can consist of any relative index of abundance.
#' @param starts Initial values for parameters
#' @param Nmax Maximum search increment for abundance considered in the optimisation
#' @param se flag to return the standard error (hessian).
#' @alpha quantile for (1-alpha) level confidence intervals
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
remGP<- function (catch, effort, index=NULL, starts, Nmax=1000, se = TRUE, alpha=0.05, ...){

  if (length(catch) != length(effort))
    stop("unequal catch/effort vector lengths.")
  x <- as.data.frame(cbind(catch, effort))
  nobs<- length(x$catch)
  names(x) <- c("catch", "effort")
  x$samp <- seq_len(nobs)
  x$cpue <- x$catch/x$effort
  x$cumcatch<- cumsum(x$catch) - x$catch
  x$cumeffort<- cumsum(x$effort) - x$effort
  R <- sum(x$catch)

  if (length(x$catch) < 3)
    stop("ml method requires at least 3 observations!")
  if(missing(starts)) starts<- c(0.01, R)

    nll2 <- function(parm) {
      k <- parm
      p <- 1 - exp(-k * x$effort)
      qp<- rep(NA, nobs)
      qp[1]<- p[1]
      for (i in 2:nobs)
        qp[i] <- p[i] * (1-p[i])^(i-1)
      Q <- prod(1-p)
      pr <- x$catch * log(qp/(1 - Q))
      MP <- sum(pr)
      RF <- lgamma(R + 1)
      CF <- sum(lgamma(x$catch + 1))
      ll <- RF - CF + MP
      (-1)*ll
    }
    nll1<- function(parm, idx=FALSE) {
      N<- parm[1]
      p <- 1 - exp(-k * x$effort)
      Q <- prod(1-p)
      ll<- lgamma(N + 1) - (lgamma((N - R) + 1) + lgamma(R + 1)) + R*log(1-Q) + (N-R)*log(Q)
      if(idx){
        pm<- plogis(parm[2])
        Nr<- N - x$cumcatch
        lli<- sum(dpois(index, Nr*pm, log=TRUE))
      }
      else lli<- 0
      (-1)*(ll + lli)
    }
    # Estimate of catch coefficient lambda (k)
    upper <- -log(1e-6)/max(x$effort)
    m1 <- optim(starts[1], nll2, lower = 0, upper = upper, method="Brent", hessian=se)
    k <- m1$par
    var.k<- invertHessian(m1, 1, se)

    # conditional LL of abundance | lambda
    if(!is.null(index)){
      if(length(index) != nobs) stop("Index data must be same length as catch data")
      else {
        cat("Index data detected - adding index-removal estimation","\n\n")
        nP<- 2
        starts<- c(starts, 0)
        m2 <- optim(starts[2:3], nll1, idx=TRUE, method="BFGS", hessian=se)
        ests<- m2$par
        covMat<- invertHessian(m2, nP, se)
      }
    }
    else {
      cat("Gould and Pollock removal estimator","\n\n")
      nP<- 1
      upper <- R+Nmax
      m2 <- optim(starts[2], nll1, idx=FALSE, lower = R, upper = upper, method="Brent", hessian=se)
      ests<- m2$par
      covMat<- invertHessian(m2, nP, se)
    }

    typeNames<- c("state","catch")

    state <- list(name = "Abundance", short.name = "N",
                  estimates = ests[1],
                  covMat = as.matrix(covMat[1,1]),
                  invlink = "identLink",
                  invlinkGrad = "identLink.grad")

    catch <- list(name = "catchability", short.name = "lambda",
                estimates = k,
                covMat = var.k,
                invlink = "identLink",
                invlinkGrad = "identLink.grad")
    if(nP==2) {
      typeNames<- c(typeNames,"det")
      det<- list(name = "detection", short.name = "p",
           estimates = ests[2],
           covMat = as.matrix(covMat[2,2]),
           invlink = "logistic",
           invlinkGrad = "logistic.grad")
    }
    else
      det<- NULL


    efit <- list( fitType = "remGP", types=typeNames,
                  state=state, catch=catch, det=det, opt = list(m1=m1,m2=m2), data=x)
    class(efit) <- c('efitGP','efit','list')

    return(efit)
}






