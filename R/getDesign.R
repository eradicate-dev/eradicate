
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
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    } else {
        siteCovs <- siteCovs(emf)
        obsCovs <- siteCovs[rep(1:M, each = R),]
    }

    ## Record future column names
    colNames <- colnames(siteCovs)
    colnames(obsCovs) <- colNames

    ## add observation number
    obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

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

  ## Record future column names for obsCovs
  colNames <- colnames(siteCovs)

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

    # add observation level vars
    obsCovs <- data.frame(placeHolder = rep(1, M*R))
      # add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(primaryCovs))
    obsCovs <- cbind(obsCovs, primaryCovs[rep(1:(M*T), each = R/T),])
    colnames(obsCovs) <- cnames
     # add observation number if not present
    if(!("obsNum" %in% names(obsCovs)))
      obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

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

  # add observation level vars
  obsCovs <- data.frame(placeHolder = rep(1, M*R))
  # add site and primarycovs covariates, which contain siteCovs
  cnames <- c(colnames(obsCovs), colnames(primaryCovs))
  obsCovs <- cbind(obsCovs, primaryCovs[rep(1:(M*T), each = R/T),])
  colnames(obsCovs) <- cnames
  # add observation number if not present
  if(!("obsNum" %in% names(obsCovs)))
    obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

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
  } else {
    removed.sites=integer(0)
  }

  return(list(X = X, X.offset = X.offset, removed.sites = removed.sites))
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
  ym.long.na <- is.na(ym.long)

  covs.na <- apply(cbind(X.long.na, Xdet.long.na, Xdetm.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na & !ym.long.na

  if(sum(y.new.na) > 0) {
    y.long[y.new.na] <- NA
    ym.long[y.new.na]<- NA
    warning("Some observations have been discarded because correspoding
            covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, numY(emf), byrow = TRUE)
  ym<- matrix(ym.long, M, numY(emf), byrow=TRUE)

  y.to.remove <- apply(y, 1, function(x) all(is.na(x)))
  ym.to.remove <- apply(ym, 1, function(x) all(is.na(x)))
  all.sites.to.remove<- as.logical(y.to.remove + ym.to.remove)

  num.to.remove <- sum(all.sites.to.remove)
  if(num.to.remove > 0) {
    y <- y[!y.to.remove,, drop = FALSE]
    ym<- ym[!all.sites.to.remove,, drop = FALSE]
    Xlam <- Xlam[!y.to.remove,, drop = FALSE]
    Xlam.offset <- Xlam.offset[!y.to.remove]
    Xphi <- Xphi[!y.to.remove[rep(1:M, each = T)],, drop = FALSE]
    Xphi.offset <- Xphi.offset[!y.to.remove[rep(1:M, each = T)]]
    Xdet <- Xdet[!y.to.remove[rep(1:M, each = R)],,
                 drop=FALSE]
    Xdet.offset <- Xdet.offset[!y.to.remove[rep(1:M, each=R)]]
    Xdetm <- Xdetm[!all.sites.to.remove[rep(1:M, each = R)],,
                 drop=FALSE]
    Xdetm.offset <- Xdetm.offset[!all.sites.to.remove[rep(1:M, each=R)]]
    warning(paste(num.to.remove,
                  "sites have been discarded because of missing data."), call.=FALSE)
  }
  list(y = y, ym=ym, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi,
       Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
       Xdetm = Xdetm, Xdetm.offset = Xdetm.offset, removed.sites = which(all.sites.to.remove))
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

