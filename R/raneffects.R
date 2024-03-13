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
  if(missing(K)) {
    warning("K was set to max(y)+50 by default")
    K <- max(y, na.rm=TRUE)+50
  }
  srm <- obj$sitesRemoved
  preds<- calcN(obj)
  lam <- preds$cellpreds$N
  if(length(srm) > 0) {
    y <- y[-obj$sitesRemoved,]
    lam <- lam[-obj$sitesRemoved]
  }
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
raneffects.efitGRM<- function(obj, K, ...) {

  detformula <- as.formula(obj$detformula)
  lamformula <- as.formula(obj$lamformula)
  mdetformula<- as.formula(obj$mdetformula)
  data <- obj$data
  D <- getDesign(data, lamformula, detformula, mdetformula)

  Xlam <- D$Xlam
  y <- D$y
  ym <- D$ym

  if(missing(K)) {
    K <- obj$K
  }

  Xlam.offset <- D$Xlam.offset

  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))

  beta.lam <- coef(obj, type="state")

  lambda <- exp(Xlam %*% beta.lam + Xlam.offset)

  allp<- calcP(obj)
  cp <- allp$cp
  cp[is.na(y)] <- NA

  cpm <- allp$cpm
  cpm[is.na(y)]<- NA

  N <- 0:K

  M <- nrow(y)
  R <- ncol(y)

  cumy<- matrix(NA_real_, M, R)
  for(i in 1:M){
      if(all(is.na(y[i,])))
        cumy[i,]<- NA
      else {
        for(r in 1:R){
          cumy[i,r]<- sum(y[i,1:r], na.rm=TRUE) - y[i,r]
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

    if(all(is.na(y[i,])) | all(is.na(ym[i,])))
        next
      for(k in 1:(K+1)) {
        y.i <- y[i,]
        ym.i <- ym[i,]
        ydot <- N[k] - sum(y.i, na.rm=TRUE)
        y.i <- c(y.i, ydot)
        Nr <- N[k] - cumy[i,]
        if(ydot < 0 | any(Nr < 0)) {
          g[k] <- 0
          h[k] <- 0
          next
        }
        cp.i <- cp[i,]
        cp.i <- c(cp.i, 1-sum(cp.i, na.rm=TRUE))
        na.i <- is.na(cp.i)
        y.i[na.i] <- NA
        g[k] <- g[k]*dmultinom(y.i[!na.i], N[k], cp.i[!na.i])

        cpm.i<- cpm[i,]
        nam.i<- is.na(cpm.i)
        h[k]<- h[k]*exp(sum(dpois(ym.i[!nam.i], Nr*cpm.i[!nam.i], log=TRUE)))
      }

    ml <- f*g*h
    post[i,,1] <- ml/sum(ml)
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#-----------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitMNS<- function(obj, K, ...) {
  data<- obj$data
  y <- obj$data$y
  srm <- obj$sitesRemoved
  if(length(srm) > 0)
    y <- y[-obj$sitesRemoved,]
  if(missing(K)) {
    warning("K was set to max(y)+50 by default")
    K <- max(y, na.rm=TRUE)+50
  }
  mix <- obj$mixture
  M <- data$numSites
  T <- data$numPrimary
  J <- data$numSecondary
  ya<- data$ya
  ya<- aperm(ya, c(1, 3, 2))

  preds<- calcN(obj)
  lam <- preds$cellpreds$N
  lam <- matrix(lam, M)

  cp <- calcP(obj)
  cp[is.na(y)] <- NA
  cp <- cbind(cp, 1-rowSums(cp, na.rm=TRUE))
  cp<- array(t(cp), c(J+1, T, M))
  cp <- aperm(cp, c(3, 1, 2))

  N <- 0:K
  post <- array(0, c(M, K+1, T))
  colnames(post) <- N
  for(i in 1:M) {
    for(t in 1:T) {
      if(identical(mix, "NB"))
        alpha <- exp(coef(obj, type="alpha"))
      switch(mix,
             P  = f <- dpois(N, lam[i,t]),
             NB = f <- dnbinom(N, mu=lam[i,t], size=alpha))
      g <- rep(1, K+1)
      if(all(is.na(ya[i,,t])) | all(is.na(cp[i,,t])))
        next
      for(k in 1:(K+1)) {
        yi <- ya[i,,t]
        ydot <- N[k] - sum(yi)
        if(ydot<0) {
          g[k] <- 0
          next
        }
        yi <- c(yi, ydot)
        cp.i<- cp[i,,t]
        na.i<- is.na(cp.i)
        g[k] <- g[k] * dmultinom(yi[!na.i], size=N[k], prob=cp.i[!na.i])
      }
      ml <- f*g
      post[i,,t] <- ml / sum(ml)
    }
  }
  class(post)<- c("raneffects",class(post))
  return(post)
}

#-----------------------------------------
#' @rdname raneffects
#' @export
raneffects.efitGRMS<- function(obj, K, ...) {
  data<- obj$data
  y <- obj$data$y
  srm <- obj$sitesRemoved
  if(length(srm) > 0)
    y <- y[-obj$sitesRemoved,]
  if(missing(K)) {
    warning("K was set to max(y)+50 by default")
    K <- max(y, na.rm=TRUE)+50
  }
  mix <- obj$mixture
  M <- data$numSites
  T <- data$numPrimary
  J <- data$numSecondary
  ya<- data$ya
  ya<- aperm(ya, c(1, 3, 2))

  yma<- data$yma
  yma<- aperm(yma, c(1,3,2))

  preds<- calcN(obj)
  lam <- preds$cellpreds$N
  lam <- matrix(lam, M)

  allp<- calcP(obj)

  cp <- allp$cp
  cp[is.na(y)] <- NA
  cp <- cbind(cp, 1-rowSums(cp, na.rm=TRUE))
  cp<- array(t(cp), c(J+1, T, M))
  cp <- aperm(cp, c(3, 1, 2))

  cpm <- allp$cpm
  cpm[is.na(y)]<- NA
  cpm<- array(t(cpm), c(J, T, M))
  cpm <- aperm(cpm, c(3, 1, 2))

  cumy<- array(NA_real_, c(M, J, T))
  for(i in 1:M){
    for(t in 1:T) {
      if(all(is.na(ya[i,,t])))
        cumy[i,,t]<- NA
      else {
        for(j in 1:J){
          cumy[i,j,t]<- sum(ya[i,1:j,t], na.rm=TRUE) - ya[i,j,t]
        }
      }
    }
  }

  N <- 0:K
  post <- array(0, c(M, K+1, T))
  colnames(post) <- N
  for(i in 1:M) {
    for(t in 1:T) {
      if(identical(mix, "NB"))
        alpha <- exp(coef(obj, type="alpha"))
      switch(mix,
             P  = f <- dpois(N, lam[i,t]),
             NB = f <- dnbinom(N, mu=lam[i,t], size=alpha))
      g <- rep(1, K+1)
      h <- rep(1, K+1)
      if(all(is.na(ya[i,,t])) | all(is.na(cp[i,,t])))
        next
      for(k in 1:(K+1)) {
        yi <- ya[i,,t]
        ydot <- N[k] - sum(yi)
        Nr <- N[k] - cumy[i,,t]
        if(ydot < 0 | any(Nr < 0)) {
          g[k] <- 0
          h[k] <- 0
          next
        }
        cp.i<- cp[i,,t]
        na.i<- is.na(cp.i)
        yi <- c(yi, ydot)
        yi[na.i]<- NA
        g[k] <- g[k] * dmultinom(yi[!na.i], size=N[k], prob=cp.i[!na.i])
        # now index data
        yma.i<- yma[i,,t]
        cpm.i<- cpm[i,,t]
        nam.i<- is.na(cpm.i)
        h[k]<- h[k]*exp(sum(dpois(yma.i[!nam.i], Nr*cpm.i[!nam.i], log=TRUE)))
      }
      ml <- f*g*h
      post[i,,t] <- ml / sum(ml)
    }
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


