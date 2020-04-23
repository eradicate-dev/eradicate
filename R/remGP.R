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
#' @param alpha quantile for (1-alpha) level confidence intervals
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
remGP<- function (data, starts, se = TRUE, ...){
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
  if(missing(starts)) {
    cstart<- -log(max(x$effort))
      if(!is.null(x$index)) {
        istart<- -log(max(x$ieffort))
        starts<- c(log(R+1), cstart, istart)
      }
      else {
        starts<- c(log(R+1), cstart)
      }
  }

    nll <- function(parm, idx=FALSE) {
      lambda<- exp(parm[1])
      p <- 1 - exp(-exp(parm[2] + log(x$effort)))
      pi<- removalPiFun(p)
      ll<- sum(dpois(x$catch, lambda*pi, log=TRUE))
      if(idx) {
        pm <- 1 - exp(-exp(parm[3] + log(x$ieffort)))
        pmi<- removalPiFun(pm)
        lli<- sum(dpois(x$index, lambda*pmi, log=TRUE))
      }
      else lli<- 0
      return((-1)*(ll+lli))
    }

    if(!is.null(x$index)) {
      nP<- 3
      lb<- c(log(R+1),-20,-20)
      ub<- c(log(R*10),2,2)
      m <- optim(starts, nll, idx=TRUE, method="L-BFGS-B", lower=lb, upper=ub, hessian=se)
    }
    else {
      nP<- 2
      lb<- c(log(R+1),-20)
      ub<- c(log(R*10),2)
      m <- optim(starts, nll, idx=FALSE, method="L-BFGS-B", lower=lb, upper=ub, hessian=se)
    }

    ests<- m$par
    covMat<- invertHessian(m, nP, se)

    state <- list(name = "Abundance", short.name = "N",
                  estimates = ests[1],
                  covMat = as.matrix(covMat[1,1]),
                  invlink = "exp",
                  invlinkGrad = "exp")

    catch <- list(name = "catchability", short.name = "lambda",
                estimates = ests[2],
                covMat = as.matrix(covMat[2,2]),
                invlink = "cloglog",
                invlinkGrad = "cloglog.grad")

    estimates<- list(state=state, catch=catch)
    if(nP==3) {
      det<- list(name = "detection", short.name = "p",
           estimates = ests[3],
           covMat = as.matrix(covMat[3,3]),
           invlink = "cloglog",
           invlinkGrad = "cloglog")
      estimates$det<- det
    }

    efit <- list(fitType = "remGP", estimates=estimates, opt = m, data=x)
    class(efit) <- c('efitGP','efit','list')

    return(efit)
}






