
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

getDesign.eFrameGR<- function(emf, lamformula, phiformula, detformula, na.rm = TRUE) {

    M <- numSites(emf)
    T <- emf$numPrimary
    R <- numY(emf) # 2*T for double observer sampling
    # 1*T for distance sampling
    # nPasses*T for removal sampling

    ## Compute phi design matrices
    if(is.null(emf$primaryCovs)) {
      primaryCovs <- data.frame(placeHolder = rep(1, M*T))
    } else primaryCovs <- emf$primaryCovs


    # add siteCovs in so they can be used as well
    if(!is.null(emf$siteCovs)) {
      sC <- emf$siteCovs[rep(1:M, each = T),,drop=FALSE]
      primaryCovs <- cbind(primaryCovs, sC)
    }

    Xphi.mf <- model.frame(phiformula, primaryCovs, na.action = NULL)
    Xphi <- model.matrix(phiformula, Xphi.mf)
    Xphi.offset <- as.vector(model.offset(Xphi.mf))
    if(!is.null(Xphi.offset)) Xphi.offset[is.na(Xphi.offset)] <- 0

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

    # add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(primaryCovs))
    obsCovs <- cbind(obsCovs, primaryCovs[rep(1:(M*T), each = R/T),])
    colnames(obsCovs) <- cnames

    Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    Xdet <- model.matrix(detformula, Xdet.mf)
    Xdet.offset <- as.vector(model.offset(Xdet.mf))
    if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0

    if(na.rm)
      out <- handleNA(emf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet,
                      Xdet.offset)
    else
      out <- list(y=getY(emf), Xlam=Xlam, Xlam.offset = Xlam.offset,
                  Xphi=Xphi, Xphi.offset = Xphi.offset, Xdet=Xdet,
                  removed.sites=integer(0))

    return(list(y = out$y, Xlam = out$Xlam, Xphi = out$Xphi,
                Xdet = out$Xdet,
                Xlam.offset = out$Xlam.offset,
                Xphi.offset = out$Xphi.offset,
                Xdet.offset = out$Xdet.offset,
                removed.sites = out$removed.sites))
  }

#------------------------
getDesign.eFrameGRM<- function(emf, lamformula, phiformula, detformula, mdetformula, na.rm = TRUE) {

  M <- numSites(emf)
  T <- emf$numPrimary
  R <- numY(emf) # 2*T for double observer sampling
  # 1*T for distance sampling
  # nPasses*T for removal sampling

  ## Compute phi design matrices
  if(is.null(emf$primaryCovs)) {
    primaryCovs <- data.frame(placeHolder = rep(1, M*T))
  } else primaryCovs <- emf$primaryCovs

  # add siteCovs in so they can be used as well
  if(!is.null(emf$siteCovs)) {
    sC <- emf$siteCovs[rep(1:M, each = T),,drop=FALSE]
    primaryCovs <- cbind(primaryCovs, sC)
  }

  Xphi.mf <- model.frame(phiformula, primaryCovs, na.action = NULL)
  Xphi <- model.matrix(phiformula, Xphi.mf)
  Xphi.offset <- as.vector(model.offset(Xphi.mf))
  if(!is.null(Xphi.offset)) Xphi.offset[is.na(Xphi.offset)] <- 0

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

  # add site and yearlysite covariates, which contain siteCovs
  cnames <- c(colnames(obsCovs), colnames(primaryCovs))
  obsCovs <- cbind(obsCovs, primaryCovs[rep(1:(M*T), each = R/T),])
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
    out <- handleNA(emf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet,
                    Xdet.offset, Xdetm, Xdetm.offset)
  else
    out <- list(y=emf$y, ym=emf$ym, Xlam=Xlam, Xlam.offset = Xlam.offset,
                Xphi=Xphi, Xphi.offset = Xphi.offset, Xdet=Xdet, Xdetm=Xdetm,
                Xdetm.offset = Xdetm.offset, removed.sites=integer(0))

  return(list(y = out$y, ym=out$ym, Xlam = out$Xlam, Xphi = out$Xphi,
              Xdet = out$Xdet, Xdetm=out$Xdetm,
              Xlam.offset = out$Xlam.offset,
              Xphi.offset = out$Xphi.offset,
              Xdet.offset = out$Xdet.offset,
              Xdetm.offset = out$Xdetm.offset,
              removed.sites = out$removed.sites))
}
#---------------------------------------------
# Method for for occuMS
getDesign.eFrameMS<- function(emf, lamformula, gamformula, epsformula, detformula, na.rm = TRUE) {

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
# Multinomial open models
getDesign.eFrameMNO <- function(emf, lamformula, gamformula, omformula, detformula,
                                iotaformula, na.rm = TRUE) {

  y <- emf$y
  M <- nrow(y)
  T <- emf$numPrimary
  J <- ncol(y) / T
  delta <- emf$primaryPeriod

  if(is.null(emf$primaryCovs)) {
    primaryCovs <- data.frame(placeHolder = rep(1, M*T))
  } else primaryCovs <- emf$primaryCovs

  # add siteCovs in so they can be used as well
  if(!is.null(emf$siteCovs)) {
    sC <- emf$siteCovs[rep(1:M, each = T),,drop=FALSE]
    primaryCovs <- cbind(primaryCovs, sC)
  }

  if(is.null(siteCovs(emf)))
    siteCovs <- data.frame(placeHolder = rep(1, M))
  else
    siteCovs <- siteCovs(emf)

  Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
  Xlam <- model.matrix(lamformula, Xlam.mf)
  Xlam.offset <- as.vector(model.offset(Xlam.mf))
  if(!is.null(Xlam.offset))
    Xlam.offset[is.na(Xlam.offset)] <- 0

  # add observation level vars
  if(is.null(obsCovs(emf)))
    obsCovs <- data.frame(obsNum = as.factor(rep(1:(J*T), M)))
  else{
    obsCovs <- obsCovs(emf)
  }
  # add site and primarycovs covariates, which contain siteCovs
  cnames <- c(colnames(obsCovs), colnames(primaryCovs))
  obsCovs <- cbind(obsCovs, primaryCovs[rep(1:(M*T), each = J),])
  colnames(obsCovs) <- cnames

  # Ignore last year of data
  transCovs <- primaryCovs[-seq(T, M*T, by=T),,drop=FALSE]
  for(i in 1:ncol(transCovs))
    if(is.factor(transCovs[,i]))
      transCovs[,i] <- factor(transCovs[,i]) # drop unused levels

  Xiota.mf <- model.frame(iotaformula, transCovs, na.action = NULL)
  Xiota <- model.matrix(iotaformula, Xiota.mf)
  Xiota.offset <- as.vector(model.offset(Xiota.mf))
  if(!is.null(Xiota.offset))
    Xiota.offset[is.na(Xiota.offset)] <- 0
  Xp.mf <- model.frame(detformula, obsCovs, na.action = NULL)
  Xp <- model.matrix(detformula, Xp.mf)
  Xp.offset <- as.vector(model.offset(Xp.mf))
  if(!is.null(Xp.offset))
    Xp.offset[is.na(Xp.offset)] <- 0
  Xgam.mf <- model.frame(gamformula, transCovs, na.action = NULL)
  Xgam <- model.matrix(gamformula, Xgam.mf)
  Xgam.offset <- as.vector(model.offset(Xgam.mf))
  if(!is.null(Xgam.offset))
    Xgam.offset[is.na(Xgam.offset)] <- 0
  Xom.mf <- model.frame(omformula, transCovs, na.action = NULL)
  Xom <- model.matrix(omformula, Xom.mf)
  Xom.offset <- as.vector(model.offset(Xom.mf))
  if(!is.null(Xom.offset))
    Xom.offset[is.na(Xom.offset)] <- 0

  # determine if gamma, omega, and iota are scalar, vector, or matrix valued
  # Runtime is much faster for scalars and vectors
  Xgo <- cbind(Xgam, Xom, Xiota)
  getGOdims <- function(x) {
    xm <- matrix(x, M, T-1, byrow=TRUE)
    col.table <- apply(xm, 2, table)
    row.table <- apply(xm, 1, table)
    if(is.vector(col.table) & !is.list(col.table)) {
      return("rowvec")
    } else if(is.vector(row.table) & !is.list(row.table)) {
      return("colvec")
    } else
      return("matrix")
  }
  if(isTRUE(all.equal(gamformula, ~1)) & isTRUE(all.equal(omformula, ~1)) &
     isTRUE(all.equal(iotaformula, ~1)))
    go.dims <- "scalar"
  else {
    go.dims.vec <- apply(Xgo, 2, getGOdims)
    if(all(go.dims.vec == "rowvec"))
      go.dims <- "rowvec"
    else if(all(go.dims.vec == "colvec"))
      go.dims <- "matrix" ##"colvec"  ## NOTE: Temporary fix to the problem reported with time-only-varying covariates
    else
      go.dims <- "matrix"
  }

  if(na.rm)
    out <- handleNA(emf, Xlam, Xgam, Xom, Xp, Xiota,
                    Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
                    delta)
  else {   # delta needs to be formatted first
    ya <- array(y, c(M, J, T))
    yna <- apply(is.na(ya), c(1,3), all)
    delta <- formatDelta(delta, yna)
    out <- list(y=y, Xlam=Xlam, Xgam=Xgam, Xom=Xom, Xp=Xp, Xiota=Xiota,
                Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                Xiota.offset=Xiota.offset,
                delta=delta, removed.sites=integer(0))
  }

  return(list(y = out$y, Xlam = out$Xlam, Xgam = out$Xgam,
              Xom = out$Xom, Xp = out$Xp, Xiota = out$Xiota,
              Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
              Xom.offset=Xom.offset, Xp.offset=Xp.offset,
              Xiota.offset=Xiota.offset, delta = out$delta,
              removed.sites = out$removed.sites, go.dims = go.dims))
}

#--------
# Stacked data for trend analysis

getDesign.eFrameMNS<- function(emf, lamformula, detformula, na.rm = TRUE) {

  M <- numSites(emf)
  T <- emf$numPrimary
  R <- numY(emf) # 2*T for double observer sampling
  J <- R/T
  delta<- emf$delta
  y <- stack.data(emf$y, T)

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

  detformula <- as.formula(detformula)
  stateformula <- as.formula(lamformula)

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
handleNA.eFrameGR<- function(emf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, Xdet.offset) {

  M <- numSites(emf)
  T <- emf$numPrimary
  R <- numY(emf)
  J <- R/T

  obsToY <- diag(J)
  obsToY[upper.tri(obsToY)] <- 1
  obsToY <- kronecker(diag(T), obsToY)

  # treat Xphi and Xlam together
  X <- cbind(Xphi, Xlam[rep(1:M, each = T), ])

  X.na <- is.na(X)
  X.long.na <- X.na[rep(1:(M*T), each = J),]

  Xdet.long.na <- apply(Xdet, 2, function(x) {
    x.mat <- matrix(x, M, R, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% obsToY
    x.long <- as.vector(t(x.mat))
    x.long > 0
  })

  Xdet.long.na <- apply(Xdet.long.na, 1, any)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  covs.na <- apply(cbind(X.long.na, Xdet.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because
            correspoding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, numY(emf), byrow = TRUE)
  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove,, drop = FALSE]
    Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
    Xlam.offset <- Xlam.offset[!sites.to.remove]
    Xphi <- Xphi[!sites.to.remove[rep(1:M, each = T)],, drop = FALSE]
    Xphi.offset <- Xphi.offset[!sites.to.remove[rep(1:M, each = T)]]
    Xdet <- Xdet[!sites.to.remove[rep(1:M, each = R)],,
                 drop=FALSE]
    Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=R)]]
    warning(paste(num.to.remove,
                  "sites have been discarded because of missing data."), call.=FALSE)
  }
  list(y = y, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi,
       Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
       removed.sites = which(sites.to.remove))
}

#-------------------------------------
handleNA.eFrameGRM<- function(emf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, Xdet.offset,
                              Xdetm, Xdetm.offset) {

  M <- numSites(emf)
  T <- emf$numPrimary
  R <- numY(emf)
  J <- R/T

  obsToY <- diag(J)
  obsToY[upper.tri(obsToY)] <- 1
  obsToY <- kronecker(diag(T), obsToY)

  # treat Xphi and Xlam together
  X <- cbind(Xphi, Xlam[rep(1:M, each = T), ])

  X.na <- is.na(X)
  X.long.na <- X.na[rep(1:(M*T), each = J),]

  Xdet.long.na <- apply(Xdet, 2, function(x) {
    x.mat <- matrix(x, M, R, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% obsToY
    x.long <- as.vector(t(x.mat))
    x.long > 0
  })

  Xdetm.long.na <- apply(Xdetm, 2, function(x) {
    x.mat <- matrix(x, M, R, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% obsToY
    x.long <- as.vector(t(x.mat))
    x.long > 0
  })

  Xdet.long.na <- apply(Xdet.long.na, 1, any)
  Xdetm.long.na <- apply(Xdetm.long.na, 1, any)

  y.long <- as.vector(t(emf$y))
  y.long.na <- is.na(y.long)

  ym.long <- as.vector(t(emf$ym))

  covs.na <- apply(cbind(X.long.na, Xdet.long.na, Xdetm.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    ym.long[y.new.na]<- NA
    warning("Some observations have been discarded because correspoding
            covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, numY(emf), byrow = TRUE)
  ym<- matrix(ym.long, M, numY(emf), byrow=TRUE)

  sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!sites.to.remove,, drop = FALSE]
    ym<- ym[!sites.to.remove,, drop = FALSE]
    Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
    Xlam.offset <- Xlam.offset[!sites.to.remove]
    Xphi <- Xphi[!sites.to.remove[rep(1:M, each = T)],, drop = FALSE]
    Xphi.offset <- Xphi.offset[!sites.to.remove[rep(1:M, each = T)]]
    Xdet <- Xdet[!sites.to.remove[rep(1:M, each = R)],,
                 drop=FALSE]
    Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=R)]]
    Xdetm <- Xdetm[!sites.to.remove[rep(1:M, each = R)],,
                 drop=FALSE]
    Xdetm.offset <- Xdetm.offset[!sites.to.remove[rep(1:M, each=R)]]
    warning(paste(num.to.remove,
                  "sites have been discarded because of missing data."), call.=FALSE)
  }
  list(y = y, ym=ym, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi,
       Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
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

#---------------------------------------
#multinomial open

handleNA.eFrameMNO<- function(emf, Xlam, Xgam, Xom, Xp, Xiota, Xlam.offset, Xgam.offset,
                   Xom.offset, Xp.offset, Xiota.offset, delta) {

  M <- numSites(emf)
  T <- emf$numPrimary
  y <- getY(emf)
  J <- ncol(y) / T
  R <- numY(emf)

  obsToY <- diag(J)
  obsToY[upper.tri(obsToY)] <- 1
  obsToY <- kronecker(diag(T), obsToY)

  Xlam.long <- Xlam[rep(1:M, each = J*T),]
  Xlam.long.na <- is.na(Xlam.long)

  long.na <- function(x) {
    x.mat <- matrix(x, M, R, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% obsToY
    x.long <- as.vector(t(x.mat))
    x.long > 0
  }

  o2y2 <- diag(T)
  o2y2 <- o2y2[-T, -T]

  long.na2 <- function(x) {
    x.mat <- matrix(x, M, T-1, byrow = TRUE)
    x.mat <- is.na(x.mat)
    x.mat <- x.mat %*% o2y2
    x.long <- as.vector(t(x.mat))
    x.long > 0
  }

  Xp.long.na <- apply(Xp, 2, long.na)
  Xp.long.na <- apply(Xp.long.na, 1, any)

  Xgam.long.na <- apply(Xgam, 2, long.na2)
  Xgam.long.na <- apply(Xgam.long.na, 1, any)
  Xom.long.na <- apply(Xom, 2, long.na2)
  Xom.long.na <- apply(Xom.long.na, 1, any)

  y.long <- as.vector(t(y))
  y.long.na <- is.na(y.long)

  #  delta.long <- as.vector(t(delta))
  #	delta.long.na <- is.na(delta.long)

  covs.na <- apply(cbind(Xlam.long.na, Xp.long.na), 1, any)
  covs.na2 <- apply(cbind(Xgam.long.na, Xom.long.na), 1, any)
  covs.na3 <- rep(covs.na2, each=J)
  # If gamma[1, 1] is NA, remove y[1, 2]
  #common <- 1:(M*J*(T-1))
  ignore <- rep(seq(1, M*J*T, by=J*T), each=J) + 0:(J-1)
  covs.na[-ignore] <- covs.na[-ignore] | covs.na3

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    warning("Some observations have been discarded because corresponding covariates were missing.",
            call. = FALSE)
  }

  y.wide <- matrix(y.long, nrow=M, ncol=J*T, byrow = TRUE)
  #	delta <- matrix(delta.long, nrow=M, ncol=T, byrow = TRUE)
  sites.to.remove <- apply(y.wide, 1, function(x) all(is.na(x)))
  ya <- array(y.wide, c(M, J, T))
  yna <- apply(is.na(ya), c(1,3), all)
  delta <- formatDelta(delta, yna)

  num.to.remove <- sum(sites.to.remove)
  if(num.to.remove > 0) {
    y.wide <- y.wide[!sites.to.remove, ,drop = FALSE]
    Xlam <- Xlam[!sites.to.remove, ,drop = FALSE]
    Xlam.offset <- Xlam.offset[!sites.to.remove]
    Xgam <- Xgam[!sites.to.remove[rep(1:M, each = T-1)],,
                 drop = FALSE]
    Xgam.offset <- Xgam.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                               drop = FALSE]
    Xom <- Xom[!sites.to.remove[rep(1:M, each = T-1)],,
               drop = FALSE]
    Xom.offset <- Xom.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                             drop = FALSE]
    Xp <- Xp[!sites.to.remove[rep(1:M, each = J*T)],,
             drop = FALSE]
    Xp.offset <- Xp.offset[!sites.to.remove[rep(1:M, each = J*T)],,
                           drop = FALSE]
    Xiota <- Xiota[!sites.to.remove[rep(1:M, each = T-1)],,
                   drop = FALSE]
    Xiota.offset <- Xiota.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                                 drop = FALSE]
    delta <- delta[!sites.to.remove, ,drop =FALSE]
    warning(paste(num.to.remove, "sites have been discarded because of missing data."), call.=FALSE)
  }

  list(y = y.wide, Xlam = Xlam, Xgam = Xgam, Xom = Xom, Xp = Xp, Xiota = Xiota,
       Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
       Xom.offset=Xom.offset, Xp.offset=Xp.offset,
       Xiota.offset=Xiota.offset, delta = delta,
       removed.sites = which(sites.to.remove))
}

#--------
# Stacked multinomial
handleNA.eFrameMNS<- function(emf, X, X.offset, V, V.offset) {

  T <- emf$numPrimary
  y <- stack.data(emf$y, T)
  J <- ncol(y)
  M <- nrow(y)

  X.long <- X[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  V.long <- V[rep(1:M, each = J),]
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

