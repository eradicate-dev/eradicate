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
#' @usage remGP(data, starts, K, method = "BFGS", se=TRUE, ...)
#'
#' @param data \code{eFrameGP} object containing the catch (removal), effort and
#' optionally, the index data.  See \code{eFrameGP} for more details.
#' @param starts Initial values for parameters
#' @param K Integer representing upper bound for abundance for discrete integration
#' @param method optimsation method (see \code{?optim} for details)
#' @param se flag to return the standard error (hessian).

#' @return a \code{efit} model object.
#'
#' @examples
#'  y<- san_nic_rem$rem
#'  catch<- apply(y,2,sum)
#'  effort<- rep(nrow(y), length(catch))
#'  emf <- eFrameGP(catch, effort)
#'  mod <- remGP(emf)
#'  Nhat<- calcN(mod)
#'
#' @export
#'
remGP<- function (data, starts, K, method="Nelder-Mead", se = TRUE, ...){
  if(!is(data, "eFrameGP"))
    stop("Data is not a eFrameGP")
  dd<- remove_NA(data$counts)
  rr<- dd$x
  if(length(dd$narm) > 0) warning(paste0("a total of ",length(dd$narm)," rows with missing values removed"))
  rlist<- split(rr, rr$session)
  names(rlist)<- paste0("S",names(rlist))
  nP<- length(rlist)
  nobs<- sapply(rlist, nrow)
  if (any(nobs < 3))
    stop("ml method requires at least 3 observations!")
  rlist <- lapply(rlist, function(z) {z$cpue<- z$catch/z$effort
                                      z$cumcatch<- cumsum(z$catch) - z$catch
                                      z$cumeffort<- cumsum(z$effort) - z$effort
                                      return(z)
                                      })
 if(nP == 1) ans<- GPest(rlist[[1]], starts=starts, K=K, method=method, se=se)
  else {
    ans<- lapply(rlist, GPest, starts=starts, K=K, method=method, se=se)
    class(ans)<- c('efitGPlist')
  }
  return(ans)
}


GPest<- function(x, starts, K, method, se) {
  if(missing(K) || is.null(K)) K <- sum(x$catch) + 1000
  k <- 0:K
  if(missing(starts)) {
    cf <- coef(lm(cpue ~ cumcatch, data = x))
    cstart<- log(-cf[2])
    nstart<- log(-cf[1]/cf[2])
    if(!is.null(x$index)) {
      istart<- -log(max(x$ieffort, na.rm=TRUE))
      starts<- c(nstart, cstart, istart)
    }
    else {
      starts<- c(nstart, cstart)
    }
  }

  nll <- function(parm, x, idx=FALSE) {
    lambda<- exp(parm[1])
    p0<- plogis(parm[2])
    R<- sum(x$catch)
    p <- 1 - (1 - p0)^x$effort
    pi<- removalPiFun(p)
    pic<- pi/sum(pi)
    e<- dmultinom(x$catch, R, prob=pic)
    f<- dpois(R+k,lambda)
    g<- dbinom(R, R+k, sum(pi))
    ll<- log(sum(f*g)) + log(e)
    if(idx) {
      lk<- length(k)
      fin<- rep(NA, lk)
      pm <- exp(parm[3]) * x$ieffort
      for(i in 1:lk) {
        Nr <- (R+k[i]) - x$cumcatch
        fin[i]<- sum(dpois(x$index, Nr*pm), na.rm=TRUE)
      }
      lli<- log(sum(fin*f))
    }
    else lli<- 0
    return((-1)*(ll+lli))
  }

  if(!is.null(x$index)) {
    nP<- 3
    m <- optim(starts, nll, x=x, idx=TRUE, method=method, hessian=se)
  }
  else {
    nP<- 2
    m <- optim(starts, nll, x=x, idx=FALSE, method=method, hessian=se)
  }

  ests<- m$par
  covMat<- invertHessian(m, nP, se)

  state <- list(name = "Abundance", short.name = "N",
                estimates = ests[1],
                covMat = as.matrix(covMat[1,1]),
                invlink = "exp",
                invlinkGrad = "exp.grad")

  catch <- list(name = "catchability", short.name = "p",
                estimates = ests[2],
                covMat = as.matrix(covMat[2,2]),
                invlink = "logistic",
                invlinkGrad = "logistic.grad")

  estimates<- list(state=state, catch=catch)
  if(nP==3) {
    det<- list(name = "detection", short.name = "lambda",
               estimates = ests[3],
               covMat = as.matrix(covMat[3,3]),
               invlink = "exp",
               invlinkGrad = "exp.grad")
    estimates$det<- det
  }

  efit <- list(fitType = "remGP", estimates=estimates, opt = m, nllFun=nll, data=x)
  class(efit) <- c('efitGP','efit','list')
  return(efit)
}


remove_NA<- function(x){
  # Remove missing values from eFrameGP data
  c_na<- which(is.na(x$catch))
  e_na<- which(is.na(x$effort))
  inds<- c(c_na, e_na)
  if(length(inds) > 0) {
    x_narm<- x[-inds, ]
  }
  else {
    x_narm<- x
  }
  list(x=x_narm, narm=inds)
}

