#' raneffects
#'
#' @description
#' \code{raneffects} performs estimation of the marginal likelihood of N using
#' empirical Bayes methods.
#'
#' @param obj A fitted model object.
#' @param K maximum value of \code{N} to integrate over. The default rarely makes sense
#'
#' @return a \code{matrix} giving the marginal likelihood of N for each observation
#'
#' @export
#'
raneffects <- function(obj, ...){
  # method generic
  UseMethod("raneffects", obj)
}

#' postSamples
#'
#' @description
#' \code{postSamples} performs estimation of the posterior distribution of N using
#' empirical Bayes methods.
#'
#' @param obj An object of class 'raneffects'.
#'
#' @return a \code{matrix} giving the posterior samples of N
#'
#' @export
#'
postSamples <- function(obj, ...){
  # method generic
  UseMethod("postSamples", obj)
}

#-----------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitR<- function(obj, K, ...) {
  y <- obj$data$y
  srm <- obj$sitesRemoved
  if(length(srm) > 0)
    y <- y[-obj@sitesRemoved,]
  if(missing(K)) {
    warning("K was set to max(y)+50 by default")
    K <- max(y, na.rm=TRUE)+50
  }
  preds<- calcN(obj)
  lam <- preds$cellpreds$N
  R <- length(lam)
  cp <- calcP(obj)
  cp <- cbind(cp, 1-rowSums(cp))
  N <- 0:K
  post <- matrix(0, R, K+1)
  colnames(post) <- N
  for(i in 1:R) {
    f <- dpois(N, lam[i])
    g <- rep(1, K+1)
    if(any(is.na(y[i,])) | any(is.na(cp[i,])))
      next
    for(k in 1:(K+1)) {
      yi <- y[i,]
      ydot <- N[k] - sum(yi)
      if(ydot<0) {
        g[k] <- 0
        next
      }
      yi <- c(yi, ydot)
      g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,])
    }
    ml <- f*g
    if(all(ml == 0))
      post[i, ]<- g/sum(g)
    else
      post[i, ] <- ml / sum(ml)
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#' @rdname postSamples
#' @export
postSamples.raneffects<- function(obj, nsims=100, ...) {
  N <- dim(obj)[1]
  K <- dim(obj)[2]

  out <- matrix(NA, N, nsims)

  for (i in 1:N) {
    out[i, ] <- sample(0:(K-1), nsims, replace=TRUE, prob=obj[i, ])
  }
  out
}
