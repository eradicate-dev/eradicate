#' raneffects
#'
#' @description
#' \code{raneffects} performs estimation of the marginal likelihood of N using
#' empirical Bayes methods. A port of \code{ranef} in package \code{unmarked}.
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

#' blup
#'
#' @description
#' \code{blup} estimates the best linear unbiased predictor of the latent abundance
#' or occupancy state using empirical Bayes methods. A port of the similar function
#' in \code{unmarked}.
#'
#' @param obj A \code{raneffects} object.
#' #'
#' @return a \code{matrix} giving the marginal likelihood of N for each observation
#'
#' @export
#'
blup <- function(obj, ...){
  # method generic
  UseMethod("blup", obj)
}

#' postSamples
#'
#' @description
#' \code{postSamples} performs estimation of the posterior distribution of N using
#' empirical Bayes methods. A port of the similar function in \code{unmarked}.
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
    y <- y[-obj$sitesRemoved,]
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
  mix <- obj$mixture
  if(identical(mix, "NB"))
    alpha <- exp(coef(obj, type="alpha"))
  post <- array(0, c(R, K+1, 1))
  colnames(post) <- N
  for(i in 1:R) {
    switch(mix,
           P  = f <- dpois(N, lam[i]),
           NB = f <- dnbinom(N, mu=lam[i], size=alpha))
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
    post[i,,1] <- ml / sum(ml)
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}
#------------------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitGR<- function(obj, K, ...) {

  detformula <- as.formula(obj$detformula)
  phiformula <- as.formula(obj$phiformula)
  lamformula <- as.formula(obj$lamformula)
  data <- obj$data
  D <- getDesign(data, lamformula, phiformula, detformula)

  Xlam <- D$Xlam
  Xphi <- D$Xphi
  y <- D$y

  if(missing(K)) {
    K <- obj$K
  }

  Xlam.offset <- D$Xlam.offset
  Xphi.offset <- D$Xphi.offset

  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))

  beta.lam <- coef(obj, type="state")
  beta.phi <- coef(obj, type="phi")

  lambda <- exp(Xlam %*% beta.lam + Xlam.offset)
  if(is.null(beta.phi))
    phi <- rep(1, nrow(Xphi))
  else
    phi <- plogis(Xphi %*% beta.phi + Xphi.offset)

  cp<- calcP(obj)
  cp[is.na(y)] <- NA

  N <- 0:K

  M <- nrow(y)
  T <- data$numPrimary # Should be 1
  R <- ncol(y)
  J <- ncol(y) / T

  phi <- matrix(phi, M, byrow=TRUE)
  cpa <- array(cp, c(M,R,T))
  ya <- array(y, c(M,R,T))

  post <- array(0, c(M, K+1, 1))
  colnames(post) <- N
  mix <- obj$mixture

  if(identical(mix, "NB"))
    alpha <- exp(coef(obj, type="alpha"))
  for(i in 1:M) {
    switch(mix,
           P  = f <- dpois(N, lambda[i]),
           NB = f <- dnbinom(N, mu=lambda[i], size=alpha))
    g <- rep(1, K+1) # outside t loop

    for(t in 1:T) {
      if(all(is.na(ya[i,,t])) | is.na(phi[i,t]))
        next
      for(k in 1:(K+1)) {
        y.it <- ya[i,,t]
        ydot <- N[k] - sum(y.it, na.rm=TRUE)
        y.it <- c(y.it, ydot)

        if(ydot < 0) {
          g[k] <- 0
          next
        }
        cp.it <- cpa[i,,t]*phi[i,t]
        cp.it <- c(cp.it, 1-sum(cp.it, na.rm=TRUE))
        na.it <- is.na(cp.it)
        y.it[na.it] <- NA
        g[k] <- g[k]*dmultinom(y.it[!na.it], N[k], cp.it[!na.it])
      }
    }
    ml <- f*g
    post[i,,1] <- ml/sum(ml)
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#------------------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitGRM<- function(obj, K, ...) {

  detformula <- as.formula(obj$detformula)
  phiformula <- as.formula(obj$phiformula)
  lamformula <- as.formula(obj$lamformula)
  mdetformula<- as.formula(obj$mdetformula)
  data <- obj$data
  D <- getDesign(data, lamformula, phiformula, detformula, mdetformula)

  Xlam <- D$Xlam
  Xphi <- D$Xphi
  y <- D$y
  ym <- D$ym

  if(missing(K)) {
    K <- obj$K
  }

  Xlam.offset <- D$Xlam.offset
  Xphi.offset <- D$Xphi.offset

  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
  if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))

  beta.lam <- coef(obj, type="state")
  beta.phi <- coef(obj, type="phi")

  lambda <- exp(Xlam %*% beta.lam + Xlam.offset)
  if(is.null(beta.phi))
    phi <- rep(1, nrow(Xphi))
  else
    phi <- plogis(Xphi %*% beta.phi + Xphi.offset)

  allp<- calcP(obj)
  cp <- allp$cp
  cp[is.na(y)] <- NA

  cpm <- allp$cpm
  cpm[is.na(y)]<- NA

  N <- 0:K

  M <- nrow(y)
  T <- data$numPrimary # Should be 1
  R <- ncol(y)
  J <- ncol(y) / T

  phi <- matrix(phi, M, byrow=TRUE)
  cpa <- array(cp, c(M,R,T))
  cpma<- array(cpm, c(M,R,T))
  ya <- array(y, c(M,R,T))
  yma<- array(ym, c(M,R,T))

  cumy<- array(NA_real_, c(M,R,T))
  for(i in 1:M){
    for(t in 1:T) {
      if(all(is.na(ya[i,,t])))
        cumy[i,,t]<- NA
      else {
        for(r in 1:R){
          cumy[i,r,t]<- sum(ya[i,1:r,t], na.rm=TRUE) - ya[i,r,t]
        }
      }
    }
  }

  post <- array(0, c(M, K+1, 1))
  colnames(post) <- N
  mix <- obj$mixture

  if(identical(mix, "NB"))
    alpha <- exp(coef(obj, type="alpha"))
  for(i in 1:M) {
    switch(mix,
           P  = f <- dpois(N, lambda[i]),
           NB = f <- dnbinom(N, mu=lambda[i], size=alpha))
    g <- rep(1, K+1) # outside t loop
    h <- rep(1, K+1)
    for(t in 1:T) {
      if(all(is.na(ya[i,,t])) | all(is.na(yma[i,,t])) | is.na(phi[i,t]))
        next
      for(k in 1:(K+1)) {
        y.it <- ya[i,,t]
        ym.it <- yma[i,,t]
        ydot <- N[k] - sum(y.it, na.rm=TRUE)
        y.it <- c(y.it, ydot)
        Nr <- N[k] - cumy[i,,t]
        if(ydot < 0 | any(Nr < 0)) {
          g[k] <- 0
          h[k] <- 0
          next
        }
        cp.it <- cpa[i,,t]*phi[i,t]
        cp.it <- c(cp.it, 1-sum(cp.it, na.rm=TRUE))
        na.it <- is.na(cp.it)
        y.it[na.it] <- NA
        g[k] <- g[k]*dmultinom(y.it[!na.it], N[k], cp.it[!na.it])

        cpm.it<- cpma[i,,t]
        nam.it<- is.na(cpm.it)
        h[k]<- h[k]*exp(sum(dpois(ym.it[!nam.it], Nr*cpm.it[!nam.it], log=TRUE)))
      }
    }
    ml <- f*g*h
    post[i,,1] <- ml/sum(ml)
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#------------------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitMS<- function(obj, ...) {
# Dynamic (multi-season) occupancy model

  data <- obj$data
  M <- numSites(data)
  nY <- data$numPrimary
  J <- ncol(data$y)/nY
  psiParms <- coef(obj, 'state')
  detParms <- coef(obj, 'det')
  colParms <- coef(obj, 'col')
  extParms <- coef(obj, 'ext')
  D <- getDesign(data, obj$lamformula, obj$gamformula, obj$epsformula, obj$detformula)
  V.itj <- D$V
  X.it.gam <- D$X.gam
  X.it.eps <- D$X.eps
  W.i <- D$W

  y <- D$y
  y[y>1] <- 1
  ya <- array(y, c(M, J, nY))

  psiP <- plogis(W.i %*% psiParms)
  detP <- plogis(V.itj %*% detParms)
  colP <- plogis(X.it.gam  %*% colParms)
  extP <- plogis(X.it.eps %*% extParms)

  detP <- array(detP, c(J, nY, M))
  colP <- matrix(colP, M, nY, byrow = TRUE)
  extP <- matrix(extP, M, nY, byrow = TRUE)

  ## create transition matrices (phi^T)
  phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
  for(i in 1:M) {
    for(t in 1:(nY-1)) {
      phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                              1-extP[i,t]))
    }
  }

  ## first compute latent probs
  x <- array(NA, c(2, nY, M))
  x[1,1,] <- 1-psiP
  x[2,1,] <- psiP
  for(i in 1:M) {
    for(t in 2:nY) {
      x[,t,i] <- (phis[,,t-1,i] %*% x[,t-1,i])
    }
  }

  z <- 0:1
  post <- array(NA_real_, c(M, 2, nY))
  colnames(post) <- z

  for(i in 1:M) {
    for(t in 1:nY) {
      g <- rep(1, 2)
      for(j in 1:J) {
        if(is.na(ya[i,j,t]) | is.na(detP[j,t,i]))
          next
        g <- g * dbinom(ya[i,j,t], 1, z*detP[j,t,i])
      }
      tmp <- x[,t,i] * g
      post[i,,t] <- tmp/sum(tmp)
    }
  }

  class(post)<- c("raneffects",class(post))
  return(post)
}

#-----------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitMNS<- function(obj, K, ...) {
  data<- obj$data
  D<- getDesign(data, obj$lamformula, obj$detformula)
  y <- obj$data$y
  srm <- obj$sitesRemoved
  if(length(srm) > 0)
    y <- y[-obj$sitesRemoved,]
  if(missing(K)) {
    warning("K was set to max(y)+50 by default")
    K <- max(y, na.rm=TRUE)+50
  }
  M <- nrow(y)
  T <- data$numPrimary
  J <- ncol(y) / T
  ya<- array(t(y), c(J,T,M))
  ya<- aperm(ya, c(3, 1, 2))

  preds<- calcN(obj)
  lam <- preds$cellpreds$N
  lam <- matrix(lam, M)

  cp <- calcP(obj)
  cp <- cbind(cp, 1-rowSums(cp))
  cp<- array(t(cp), c(J+1, T, M))
  cp <- aperm(cp, c(3, 1, 2))

  N <- 0:K
  post <- array(0, c(M, K+1, T))
  colnames(post) <- N
  for(i in 1:M) {
    for(t in 1:T) {
      f <- dpois(N, lam[i,t])
      g <- rep(1, K+1)
      if(any(is.na(ya[i,,t])) | any(is.na(cp[i,,t])))
      next
    for(k in 1:(K+1)) {
      yi <- ya[i,,t]
      ydot <- N[k] - sum(yi)
      if(ydot<0) {
        g[k] <- 0
        next
      }
      yi <- c(yi, ydot)
      g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,,t])
    }
    ml <- f*g
    post[i,,t] <- ml / sum(ml)
    }
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#-----------------------------------
#' @rdname postSamples
#' @export
postSamples.raneffects<- function(obj, nsims=100, ...) {
  N <- dim(obj)[1]
  K <- dim(obj)[2]
  T <- dim(obj)[3]
  out <- array(NA, c(N, T, nsims))

  for (i in 1:N) {
    for(t in 1:T) {
    if(any(is.na(obj[i,,t]))) next
    out[i,t,] <- sample(0:(K-1), nsims, replace=TRUE, prob=obj[i,,t])
    }
  }
  out
}

#' @rdname blup
#' @export
blup.raneffects<- function(obj, ...) {
  re <- as.integer(colnames(obj))
  out <- apply(obj, c(1,3), function(x) sum(re*x))
  out <- drop(out)
  return(out)
}


