
# Process design matrices
getDesign <- function(obj, ...){
  # method generic
  UseMethod("getDesign", obj)
}

handleNA <- function(obj, ...){
  # method generic
  UseMethod("handleNA", obj)
}


getDesign.eFrame<- function(emf, lamformula, detformula, na.rm=TRUE) {

    detformula <- as.formula(detformula)
    stateformula <- as.formula(lamformula)

    M <- numSites(emf)
    R <- numY(emf)


    if(is.null(siteCovs(emf))) {
        siteCovs <- data.frame(placeHolder = rep(1, M))
    } else {
        siteCovs <- siteCovs(emf)
    }

    ## Compute detection design matrix
    if(is.null(obsCovs(emf))) {
      obsCovs <- data.frame(obsNum = as.factor(rep(1:R, M)))
    } else {
      obsCovs <- obsCovs(emf)
      obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
    }
    ## Record future column names
    colNames <- c(colnames(obsCovs), colnames(siteCovs))

    ## add site Covariates at observation-level
    obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
    colnames(obsCovs) <- colNames

    ## Compute design matrices
    X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
    X <- model.matrix(stateformula, X.mf)
    X.offset <- as.vector(model.offset(X.mf))
    if (!is.null(X.offset)) {
        X.offset[is.na(X.offset)] <- 0
    }

    V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    V <- model.matrix(detformula, V.mf)
    V.offset <- as.vector(model.offset(V.mf))
    if (!is.null(V.offset)) {
        V.offset[is.na(V.offset)] <- 0
    }

    if (na.rm) {
        out <- handleNA(emf, X, X.offset, V, V.offset)
        y <- out$y
        X <- out$X
        X.offset <- out$X.offset
        V <- out$V
        V.offset <- out$V.offset
        removed.sites <- out$removed.sites
    } else {
        y=getY(emf)
        removed.sites=integer(0)
    }

    return(list(y = y, X = X, X.offset = X.offset, V = V,
                V.offset = V.offset, removed.sites = removed.sites))
}

#-------------------------------------------------------------------------
# REST model

getDesign.eFrameREST<- function(emf, lamformula, na.rm=TRUE) {
  # REST model assumes cameras sample viewshed of area A with perfect detection
  # Future version will model effective viewshed using distance sampling techniques

  stateformula <- as.formula(lamformula)

  M <- numSites(emf)
  R <- numY(emf)

  ## Compute design matrix
  if(is.null(siteCovs(emf))) {
    siteCovs <- data.frame(placeHolder = rep(1, M))
  } else {
    siteCovs <- siteCovs(emf)
  }
  X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
  X <- model.matrix(stateformula, X.mf)
  X.offset <- as.vector(model.offset(X.mf))
  if (!is.null(X.offset)) {
    X.offset[is.na(X.offset)] <- 0
  }

  if (na.rm) {
    out <- handleNA(emf, X, X.offset)
    y <- out$y
    X <- out$X
    effort<- out$effort
    X.offset <- out$X.offset
    removed.sites <- out$removed.sites
  } else {
    y<- getY(emf)
    removed.sites=integer(0)
  }

  return(list(y = y, X = X, X.offset = X.offset, effort=effort, stay=emf$stay,
              cens=emf$cens, area=emf$area, removed.sites = removed.sites))
}


#---------------------------------------------------------
# DM for generalized removal estimator

getDesign.eFrameGR<- function(emf, lamformula, detformula, na.rm = TRUE) {

    detformula <- as.formula(detformula)
    lamformula <- as.formula(lamformula)

    M <- numSites(emf)
    R <- numY(emf)

    # Compute site-level design matrix for lambda
    if(is.null(siteCovs(emf))) {
      siteCovs <- data.frame(placeHolder = rep(1, M))
    } else siteCovs <- siteCovs(emf)

    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset)) Xlam.offset[is.na(Xlam.offset)] <- 0

    ## Compute detection design matrix
    if(is.null(obsCovs(emf))) {
      obsCovs <- data.frame(obsNum = as.factor(rep(1:R, M)))
    } else {
      obsCovs <- obsCovs(emf)
      obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
    }

    # add sitecovs to obscovs
    cnames <- c(colnames(obsCovs), colnames(siteCovs))
    obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
    colnames(obsCovs) <- cnames

    Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    Xdet <- model.matrix(detformula, Xdet.mf)
    Xdet.offset <- as.vector(model.offset(Xdet.mf))
    if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0

    if(na.rm)
      out <- handleNA(emf, Xlam, Xlam.offset, Xdet, Xdet.offset)
    else
      out <- list(y=getY(emf), Xlam=Xlam, Xlam.offset = Xlam.offset,
                  Xdet=Xdet, Xdet.offset =Xdet.offset,
                  removed.sites=integer(0))

    return(list(y = out$y, Xlam = out$Xlam,
                Xdet = out$Xdet,
                Xlam.offset = out$Xlam.offset,
                Xdet.offset = out$Xdet.offset,
                removed.sites = out$removed.sites))
  }

#------------------------
getDesign.eFrameGRM<- function(emf, lamformula, detformula, mdetformula, na.rm = TRUE) {

  detformula <- as.formula(detformula)
  lamformula <- as.formula(lamformula)
  mdetformula <- as.formula(mdetformula)

  M <- numSites(emf)
  R <- numY(emf)

  # Compute site-level design matrix for lambda
  if(is.null(siteCovs(emf))) {
    siteCovs <- data.frame(placeHolder = rep(1, M))
  } else siteCovs <- siteCovs(emf)
  Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
  Xlam <- model.matrix(lamformula, Xlam.mf)
  Xlam.offset <- as.vector(model.offset(Xlam.mf))
  if(!is.null(Xlam.offset)) Xlam.offset[is.na(Xlam.offset)] <- 0

  ## Compute detection design matrix
  if(is.null(obsCovs(emf))) {
    obsCovs <- data.frame(obsNum = as.factor(rep(1:R, M)))
  } else {
    obsCovs <- obsCovs(emf)
    obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
  }

  # add sitecovs to obscovs
  cnames <- c(colnames(obsCovs), colnames(siteCovs))
  obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
  colnames(obsCovs) <- cnames

  Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
  Xdet <- model.matrix(detformula, Xdet.mf)
  Xdet.offset <- as.vector(model.offset(Xdet.mf))
  if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0

  Xdetm.mf <- model.frame(mdetformula, obsCovs, na.action = NULL)
  Xdetm <- model.matrix(mdetformula, Xdetm.mf)
  Xdetm.offset <- as.vector(model.offset(Xdetm.mf))
  if(!is.null(Xdetm.offset)) Xdetm.offset[is.na(Xdetm.offset)] <- 0


  if(na.rm)
    out <- handleNA(emf, Xlam, Xlam.offset, Xdet,
                    Xdet.offset, Xdetm, Xdetm.offset)
  else
    out <- list(y=emf$y, ym=emf$ym, Xlam=Xlam, Xlam.offset = Xlam.offset,
                Xdet=Xdet, Xdetm=Xdetm, Xdetm.offset = Xdetm.offset,
                removed.sites=integer(0))

  return(list(y = out$y, ym=out$ym, Xlam = out$Xlam,
              Xdet = out$Xdet, Xdetm=out$Xdetm,
              Xlam.offset = out$Xlam.offset,
              Xdet.offset = out$Xdet.offset,
              Xdetm.offset = out$Xdetm.offset,
              removed.sites = out$removed.sites))
}

#---------------------------------------------
# Method for for occuMS
getDesign.eFrameMS<- function(emf, lamformula, gamformula, epsformula, detformula, na.rm = TRUE) {

    detformula <- as.formula(detformula)
    lamformula <- as.formula(lamformula)
    gamformula <- as.formula(gamformula)
    epsformula <- as.formula(epsformula)

    M <- numSites(emf)
    R <- numY(emf)
    nY <- emf$numPrimary
    J<- R / nY


    ## Compute default design matrix for gamma/epsilon
    seasCovs <- data.frame(.season = as.factor(rep(1:nY, M)))

    ## add siteCovs in so they can be used as well
    if(!is.null(emf$siteCovs)) {
      sC <- emf$siteCovs[rep(1:M, each = nY),,drop=FALSE]
      seasCovs <- cbind(seasCovs, sC)
    }

    ## Compute site-level design matrix for psi
    if(is.null(siteCovs(emf)))
      siteCovs <- data.frame(placeHolder = rep(1, M))
    else
      siteCovs <- siteCovs(emf)

    W.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    if(!is.null(model.offset(W.mf)))
      stop("offsets not currently allowed in occuMS", call.=FALSE)
    W <- model.matrix(lamformula, W.mf)

    ## Compute detection design matrix
    if(is.null(obsCovs(emf))) {
      obsCovs <- data.frame(obsNum = as.factor(rep(1:R, M)))
    } else {
      obsCovs <- obsCovs(emf)
      obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
    }

    ## add site and season covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(seasCovs))
    obsCovs <- cbind(obsCovs, seasCovs[rep(1:(M*nY), each = J),])
    colnames(obsCovs) <- cnames

    V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    if(!is.null(model.offset(V.mf)))
      stop("offsets not currently allowed in occuMS", call.=FALSE)
    V <- model.matrix(detformula, V.mf)

    ## in order to drop factor levels that only appear in last year,
    ## replace last year with NAs and use drop=TRUE
    seasCovs[seq(nY,M*nY,by=nY),] <- NA
    seasCovs <- as.data.frame(lapply(seasCovs, function(x) {x[,drop = TRUE]}))

    X.mf.gam <- model.frame(gamformula, seasCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.gam)))
      stop("offsets not currently allowed in occuMS", call.=FALSE)
    X.gam <- model.matrix(gamformula, X.mf.gam)
    X.mf.eps <- model.frame(epsformula, seasCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.eps)))
      stop("offsets not currently allowed in occuMS", call.=FALSE)
    X.eps <- model.matrix(epsformula, X.mf.eps)

    if(na.rm)
      out <- handleNA(emf, W, X.gam, X.eps, V)
    else
      out <- list(y=getY(emf), X.gam=X.gam, X.eps=X.eps, W=W, V=V,
                  removed.sites=integer(0))

    return(list(y = out$y, X.eps = out$X.eps, X.gam = out$X.gam, W = out$W,
                V = out$V, removed.sites = out$removed.sites))
}

#---------------------------------------------
# Stacked data for trend analysis

getDesign.eFrameMNS<- function(emf, lamformula, detformula, na.rm = TRUE) {

  detformula <- as.formula(detformula)
  stateformula <- as.formula(lamformula)

  M <- emf$numSites
  T <- emf$numPrimary
  J <- emf$numSecondary
  delta<- emf$delta
  y <- emf$y

  ## Compute default design matrix for seasonal strata and numeric trend
  season <- data.frame(.season = as.factor(rep(1:T, each = M)))
  trend<- data.frame(.trend = rep(delta, each = M))

  if(is.null(siteCovs(emf))) {
    siteCovs <- data.frame(placeHolder = rep(1, M))
  } else {
    siteCovs <- siteCovs(emf)
  }

  # Combine above
  seasCovs <- cbind(season, trend, siteCovs[rep(1:M, T),,drop=FALSE])


  ## Compute detection design matrix
  if(is.null(obsCovs(emf))) {
    obsCovs <- data.frame(obsNum = as.factor(rep(1:J, M*T)))
  } else {
    obsCovs <- obsCovs(emf)
    obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:J, M*T)))
  }

  ## add site and season covariates, which contain siteCovs
  cnames <- c(colnames(obsCovs), colnames(seasCovs))
  obsCovs <- cbind(obsCovs, seasCovs[rep(1:(M*T), each = J),])
  colnames(obsCovs) <- cnames

   ## Compute design matrices
  X.mf <- model.frame(stateformula, seasCovs, na.action = NULL)
  X <- model.matrix(stateformula, X.mf)
  X.offset <- as.vector(model.offset(X.mf))
  if (!is.null(X.offset)) {
    X.offset[is.na(X.offset)] <- 0
  }

  V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
  V <- model.matrix(detformula, V.mf)
  V.offset <- as.vector(model.offset(V.mf))
  if (!is.null(V.offset)) {
    V.offset[is.na(V.offset)] <- 0
  }

  if (na.rm) {
    out <- handleNA(emf, X, X.offset, V, V.offset)
    y <- out$y
    X <- out$X
    X.offset <- out$X.offset
    V <- out$V
    V.offset <- out$V.offset
    removed.sites <- out$removed.sites
  } else {
    removed.sites=integer(0)
  }

  return(list(y = y, X = X, X.offset = X.offset, V = V,
              V.offset = V.offset, removed.sites = removed.sites))
}

#------------------------

getDesign.eFrameGRMS<- function(emf, lamformula, detformula, mdetformula, na.rm = TRUE) {

  detformula <- as.formula(detformula)
  lamformula <- as.formula(lamformula)
  mdetformula <- as.formula(mdetformula)

  M <- emf$numSites
  T <- emf$numPrimary
  J <- emf$numSecondary
  delta<- emf$delta
  y <- emf$y
  ym <- emf$ym

  ## Compute default design matrix for seasonal strata and numeric trend
  season <- data.frame(.season = as.factor(rep(1:T, each = M)))
  trend<- data.frame(.trend = rep(delta, each = M))


  # site-level covariates
  if(is.null(siteCovs(emf))) {
    siteCovs <- data.frame(placeHolder = rep(1, M))
  } else siteCovs <- siteCovs(emf)

  # Combine above
  seasCovs <- cbind(season, trend, siteCovs[rep(1:M, T),,drop=FALSE])

  ## Compute detection design matrix
  if(is.null(obsCovs(emf))) {
    obsCovs <- data.frame(obsNum = as.factor(rep(1:J, M*T)))
  } else {
    obsCovs <- obsCovs(emf)
    obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:J, M*T)))
  }

  ## add site and season covariates, which contain siteCovs
  cnames <- c(colnames(obsCovs), colnames(seasCovs))
  obsCovs <- cbind(obsCovs, seasCovs[rep(1:(M*T), each = J),])
  colnames(obsCovs) <- cnames


  Xlam.mf <- model.frame(lamformula, seasCovs, na.action = NULL)
  Xlam <- model.matrix(lamformula, Xlam.mf)
  Xlam.offset <- as.vector(model.offset(Xlam.mf))
  if(!is.null(Xlam.offset)) Xlam.offset[is.na(Xlam.offset)] <- 0


  Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
  Xdet <- model.matrix(detformula, Xdet.mf)
  Xdet.offset <- as.vector(model.offset(Xdet.mf))
  if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0

  Xdetm.mf <- model.frame(mdetformula, obsCovs, na.action = NULL)
  Xdetm <- model.matrix(mdetformula, Xdetm.mf)
  Xdetm.offset <- as.vector(model.offset(Xdetm.mf))
  if(!is.null(Xdetm.offset)) Xdetm.offset[is.na(Xdetm.offset)] <- 0


  if(na.rm){
      out <- handleNA(emf, Xlam, Xlam.offset, Xdet,
                    Xdet.offset, Xdetm, Xdetm.offset)
      y <- out$y
      ym <- out$ym
      Xlam <- out$Xlam
      Xdet <- out$Xdet
      Xdetm <- out$Xdetm
      Xlam.offset <- out$Xlam.offset
      Xdet.offset <- out$Xdet.offset
      Xdetm.offset <- out$Xdetm.offset
      removed.sites <- out$removed.sites
    }
  else {
    removed.sites=integer(0)
  }

  return(list(y = y, ym=ym, Xlam = Xlam,
              Xdet = Xdet, Xdetm=Xdetm,
              Xlam.offset = Xlam.offset,
              Xdet.offset = Xdet.offset,
              Xdetm.offset = Xdetm.offset,
              removed.sites = removed.sites))
}
#---------------------------------------------
# Method for prediction
getDesign.efit<- function(obj, siteCovs, na.rm=TRUE) {
  # Create design matrix for prediction

  stateformula <- as.formula(obj$lamformula)
  X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
  X <- model.matrix(stateformula, X.mf)
  X.offset <- as.vector(model.offset(X.mf))
  if (!is.null(X.offset)) {
    X.offset[is.na(X.offset)] <- 0
  }

  colNames <- colnames(siteCovs)

  if (na.rm) {
    sites.to.remove <- apply(siteCovs, 1, function(x) any(is.na(x)))
    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
      X <- X[!sites.to.remove, ,drop = FALSE]
      X.offset <- X.offset[!sites.to.remove]
      warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
    }
    removed.sites <- which(sites.to.remove)
    retained.sites <- which(!sites.to.remove)
  } else {
    removed.sites<- integer(0)
    retained.sites<- 1:nrow(X)
  }

  return(list(X = X, X.offset = X.offset, retained.sites = retained.sites))
}
#-------------------------------------------------------
# Handle missing data methods
#--------------------------------------------------------
handleNA.eFrame<- function(emf, X, X.offset, V, V.offset) {

  J <- numY(emf)
  M <- numSites(emf)

  X.long <- X[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  V.long <- V[rep(1:M, each = J),]
  V.long.na <- is.na(V.long)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, J, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove, ,drop = FALSE]
    X <- X[!sites.to.remove, ,drop = FALSE]
    X.offset <- X.offset[!sites.to.remove]
    V <- V[!sites.to.remove, ,drop = FALSE]
    V.offset <- V.offset[!sites.to.remove]
    warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
  }

  list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
       removed.sites = which(sites.to.remove))
}

#-------------------------------------
handleNA.eFrameGR<- function(emf, X, X.offset, V, V.offset) {

  J <- numY(emf)
  M <- numSites(emf)

  X.long <- X[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  V.long <- V[rep(1:M, each = J),]
  V.long.na <- is.na(V.long)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, J, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove, ,drop = FALSE]
    X <- X[!sites.to.remove, ,drop = FALSE]
    X.offset <- X.offset[!sites.to.remove]
    V <- V[!sites.to.remove, ,drop = FALSE]
    V.offset <- V.offset[!sites.to.remove]
    warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
  }

  list(y = y, Xlam = X, Xlam.offset = X.offset, Xdet = V, Xdet.offset = V.offset,
       removed.sites = which(sites.to.remove))
}

#-------------------------------------
handleNA.eFrameGRM<- function(emf, Xlam, Xlam.offset, Xdet, Xdet.offset,
                              Xdetm, Xdetm.offset) {

  J <- numY(emf)
  M <- numSites(emf)

  X.long <- Xlam[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  V.long <- Xdet[rep(1:M, each = J),]
  V.long.na <- is.na(V.long)

  W.long <- Xdetm[rep(1:M, each = J),]
  W.long.na <- is.na(W.long)

  covs.na <- apply(cbind(X.long.na, V.long.na, W.long.na), 1, any)

  y.long <- as.vector(t(emf$y))
  y.long.na <- is.na(y.long)

  ym.long <- as.vector(t(emf$ym))

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    ym.long[y.new.na]<- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, J, byrow = TRUE)
  ym <- matrix(ym.long, M, J, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove,, drop = FALSE]
    ym<- ym[!sites.to.remove,, drop = FALSE]
    Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
    Xlam.offset <- Xlam.offset[!sites.to.remove]
    Xdet <- Xdet[!sites.to.remove[rep(1:M, each = J)],,
                 drop=FALSE]
    Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=J)]]
    Xdetm <- Xdetm[!sites.to.remove[rep(1:M, each = J)],,
                 drop=FALSE]
    Xdetm.offset <- Xdetm.offset[!sites.to.remove[rep(1:M, each=J)]]
    warning(paste(num.to.remove,
                  "sites have been discarded because of missing data."), call.=FALSE)
  }
  list(y = y, ym=ym, Xlam = Xlam, Xlam.offset = Xlam.offset,
       Xdet = Xdet, Xdet.offset = Xdet.offset,
       Xdetm = Xdetm, Xdetm.offset = Xdetm.offset, removed.sites = which(sites.to.remove))
}


#-------------------------
handleNA.eFrameREST<- function(emf, X, X.offset) {

  J <- numY(emf)
  M <- numSites(emf)

  X.long <- X[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  e.long <- as.vector(t(emf$effort))

  if(ncol(X) == 1) X.long.na<- as.matrix(X.long.na, nrow=1)
  covs.na <- apply(X.long.na, 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    e.long[y.new.na] <- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, J, byrow = TRUE)
  effort<- matrix(e.long, M, J, byrow=TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove, ,drop = FALSE]
    effort <- effort[!sites.to.remove, ,drop = FALSE]
    X <- X[!sites.to.remove, ,drop = FALSE]
    X.offset <- X.offset[!sites.to.remove]
    warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
  }

  list(y = y, X = X, X.offset = X.offset, effort=effort, removed.sites = which(sites.to.remove))
}

#-------------------------

handleNA.eFrameMS<- function(emf, W, X.gam, X.eps, V) {

  M <- numSites(emf)
  R <- numY(emf)
  nY <- emf$numPrimary
  J<- R / nY

  obsToY <- diag(R)

  ## treat both X's #######no: and W together
  #    X <- cbind(X.gam, X.eps, W[rep(1:M, each = nY), ])
  X <- cbind(X.gam, X.eps)

  X.na <- is.na(X)
  X.na[seq(nY,M*nY,by=nY),] <- FALSE  ## final years are unimportant.
  ## not true for W covs!!!
  W.expand <- W[rep(1:M, each=nY),,drop=FALSE]
  W.na <- is.na(W.expand)
  X.na <- cbind(X.na, W.na) # NAs in siteCovs results in removed site

  X.long.na <- X.na[rep(1:(M*nY), each = J),]

  V.long.na <- apply(V, 2, function(x) {
    x.mat <- matrix(x, M, R, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% obsToY
    x.long <- as.vector(t(x.mat))
    x.long > 0
  })
  V.long.na <- apply(V.long.na, 1, any)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  # It doesn't make sense to combine X.gam/eps with W here b/c
  # a X.eps does not map correctly to y
  covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because correspoding
            covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, R, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove, ,drop = FALSE]
    X.gam <- X.gam[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
    X.eps <- X.eps[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
    W <- W[!sites.to.remove,, drop = FALSE] # !!! Recent bug fix
    V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
    warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
  }
  list(y = y, X.gam = X.gam, X.eps = X.eps, W = W, V = V,
       removed.sites = which(sites.to.remove))
}

#--------
# Stacked multinomial models

handleNA.eFrameMNS<- function(emf, X, X.offset, V, V.offset) {

  M <- emf$numSites
  T <- emf$numPrimary
  J <- emf$numSecondary
  y <- emf$y

  X.long <- X[rep(1:(M*T), each = J),]
  X.long.na <- is.na(X.long)

  V.long <- V[rep(1:(M*T), each = J),]
  V.long.na <- is.na(V.long)

  y.long <- as.vector(t(y))
  y.long.na <- is.na(y.long)

  covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, (M*T), J, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove, ,drop = FALSE]
    X <- X[!sites.to.remove, ,drop = FALSE]
    X.offset <- X.offset[!sites.to.remove]
    V <- V[!sites.to.remove, ,drop = FALSE]
    V.offset <- V.offset[!sites.to.remove]
    warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
  }

  list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
       removed.sites = which(sites.to.remove))
}

#-------------------------------------
handleNA.eFrameGRMS<- function(emf, Xlam, Xlam.offset, Xdet, Xdet.offset,
                               Xdetm, Xdetm.offset) {

  M <- emf$numSites
  T <- emf$numPrimary
  J <- emf$numSecondary
  y <- emf$y
  ym<- emf$ym

  X.long <- Xlam[rep(1:(M*T), each = J),]
  X.long.na <- is.na(X.long)

  V.long <- Xdet[rep(1:(M*T), each = J),]
  V.long.na <- is.na(V.long)

  W.long <- Xdetm[rep(1:(M*T), each = J),]
  W.long.na <- is.na(W.long)

  covs.na <- apply(cbind(X.long.na, V.long.na, W.long.na), 1, any)

  y.long <- as.vector(t(y))
  y.long.na <- is.na(y.long)

  ym.long <- as.vector(t(ym))

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    ym.long[y.new.na]<- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, (M*T), J, byrow = TRUE)
  ym <- matrix(ym.long, (M*T), J, byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove,, drop = FALSE]
    ym<- ym[!sites.to.remove,, drop = FALSE]
    Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
    Xlam.offset <- Xlam.offset[!sites.to.remove]
    Xdet <- Xdet[!sites.to.remove[rep(1:M, each = J)],,
                 drop=FALSE]
    Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=J)]]
    Xdetm <- Xdetm[!sites.to.remove[rep(1:M, each = J)],,
                   drop=FALSE]
    Xdetm.offset <- Xdetm.offset[!sites.to.remove[rep(1:M, each=J)]]
    warning(paste(num.to.remove,
                  "sites have been discarded because of missing data."), call.=FALSE)
  }
  list(y = y, ym=ym, Xlam = Xlam, Xlam.offset = Xlam.offset,
       Xdet = Xdet, Xdet.offset = Xdet.offset,
       Xdetm = Xdetm, Xdetm.offset = Xdetm.offset, removed.sites = which(sites.to.remove))
}
