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
remGP<- function (data, starts, Nmax=1000, se = TRUE, ...){
  if(!is(data, "eFrameGP"))
    stop("Data is not a eFrameGP")
  x<- data$counts
  nobs<- nrow(x)
  x$samp <- seq_len(nobs)
  x$cpue <- x$catch/x$effort
  x$cumcatch<- cumsum(x$catch) - x$catch
  x$cumeffort<- cumsum(x$effort) - x$effort
  R <- sum(x$catch)

  if (nobs < 3)
    stop("ml method requires at least 3 observations!")
  cstart<- -log(max(x$effort))
  if(!is.null(x$index)) istart<- -log(max(x$ieffort))
  else istart<- 0

  if(missing(starts)) starts<- c(cstart, log(R), istart)

    nll2 <- function(parm) {
      p <- 1 - exp(-exp(parm + log(x$effort)))
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
    nll1<- function(parm, k, idx=FALSE) {
      N<- exp(parm[1])
      p <- 1 - exp(-k * x$effort)
      Q <- prod(1-p)
      ll<- lgamma(N + 1) - (lgamma((N - R) + 1) + lgamma(R + 1)) + R*log(1-Q) + (N-R)*log(Q)
      if(idx){
        pm<- exp(parm[2])*x$ieffort
        Nr<- N - x$cumcatch
        llvec<- dpois(x$index, Nr*pm, log=TRUE)
        llvec[!is.finite(llvec)]<- 0
        lli<- sum(llvec)
      }
      else lli<- 0
      (-1)*(ll + lli)
    }
    # Estimate of catch coefficient lambda (k)
    m1 <- optim(starts[1], nll2, lower= -20, upper = 2, method="Brent", hessian=se)
    k <- exp(m1$par)
    var.k<- invertHessian(m1, 1, se)

    # conditional LL of abundance | lambda
    if(!is.null(x$index)){
        cat("Index data detected - adding index-removal estimation","\n\n")
        nP<- 2
        m2 <- optim(starts[2:3], nll1, k=k, idx=TRUE, method="L-BFGS-B", lower=c(log(R), -20),
                    upper=c(log(R+Nmax), 20), hessian=se)
        ests<- m2$par
        covMat<- invertHessian(m2, nP, se)
    }
    else {
      cat("Gould and Pollock removal estimator","\n\n")
      nP<- 1

      m2 <- optim(starts[2], nll1, k=k, idx=FALSE, lower = log(R), upper = log(R+Nmax),
                  method="Brent", hessian=se)
      ests<- m2$par
      covMat<- invertHessian(m2, nP, se)
    }

    typeNames<- c("state","catch")

    state <- list(name = "Abundance", short.name = "N",
                  estimates = ests[1],
                  covMat = as.matrix(covMat[1,1]),
                  invlink = "exp",
                  invlinkGrad = "exp")

    catch <- list(name = "catchability", short.name = "lambda",
                estimates = m1$par,
                covMat = var.k,
                invlink = "cloglog",
                invlinkGrad = "cloglog.grad")
    if(nP==2) {
      typeNames<- c(typeNames,"det")
      det<- list(name = "detection", short.name = "p",
           estimates = ests[2],
           covMat = as.matrix(covMat[2,2]),
           invlink = "cloglog",
           invlinkGrad = "cloglog.grad")
    }
    else
      det<- NULL


    efit <- list( fitType = "remGP", types=typeNames,
                  state=state, catch=catch, det=det, opt = list(m1=m1,m2=m2), data=x)
    class(efit) <- c('efitGP','efit','list')

    return(efit)
}






