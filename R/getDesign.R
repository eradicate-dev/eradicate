
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
# Removal Monitor data for eFrameRM
# just need to handle extra monitoring data

getDesign.eFrameRM<- function(emf, lamformula, detformula, na.rm=TRUE) {

      detformula <- as.formula(detformula)
      stateformula <- as.formula(lamformula)

      M <- numSites(emf)
      R <- numY(emf)
      Z <- emf$Z
      ym<- emf$ym

      ## Compute removal design matrix
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

      ## Compute monitoring detection design matrix (intercept only)
      obsCovs <- data.frame(placeHolder = rep(1, M*R))
      Vm.mf <- model.frame(as.formula(~1), obsCovs, na.action = NULL)
      Vm <- model.matrix(as.formula(~1), Vm.mf)

      if (na.rm) {
        out <- handleNA(emf, X, X.offset, V, V.offset, Vm)
        y <- out$y
        ym<- out$ym
        X <- out$X
        X.offset <- out$X.offset
        V <- out$V
        V.offset <- out$V.offset
        Xm<- out$Xm
        Xm.offset<- out$Xm.offset
        Vm<- Vm
        y.sites<- out$y.sites
        ym.sites<- out$ym.sites
        removed.sites <- out$removed.sites
      } else {
        y<- getY(emf)
        Xm<- X
        Xm.offset<- X.offset
        y.sites<- seq_len(M)
        ym.sites<- y.sites
        removed.sites=integer(0)
      }

      return(list(y = y, ym=ym, X = X, X.offset = X.offset, V = V,
                  V.offset = V.offset, Xm=Xm, Xm.offset=Xm.offset, Vm=Vm, Z=Z,
                  y.sites=y.sites, ym.sites=ym.sites, removed.sites = removed.sites))
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
#--------------------------------------------------------
handleNA.eFrameRM<- function(emf, X, X.offset, V, V.offset, Vm) {

  J <- numY(emf)
  M <- numSites(emf)
  sites<- seq_len(M)

  X.long <- X[rep(1:M, each = J),]
  X.long.na <- is.na(X.long)

  V.long <- V[rep(1:M, each = J),]
  V.long.na <- is.na(V.long)

  y.long <- as.vector(t(getY(emf)))
  y.long.na <- is.na(y.long)

  ym.long <- as.vector(t(emf$ym))
  ym.long.na <- is.na(ym.long)

  covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

  ## are any NA in covs not in y already?
  y.new.na <- covs.na & !y.long.na
  ym.new.na <- covs.na & !ym.long.na

  if(sum(y.new.na) > 0 | sum(ym.new.na) > 0) {
    y.long[y.new.na] <- NA
    ym.long[ym.new.na] <- NA
    warning("Some observations have been discarded because
                corresponding covariates were missing.", call. = FALSE)
  }

  y <- matrix(y.long, M, J, byrow = TRUE)
  y.remove <- apply(y, 1, function(x) all(is.na(x)))

  ym <- matrix(ym.long, M, J, byrow = TRUE)
  ym.remove <- apply(ym, 1, function(x) all(is.na(x)))

  y <- y[!y.remove, ,drop = FALSE]
  Xr <- X[!y.remove, ,drop = FALSE]
  Xr.offset <- X.offset[!y.remove]
  V <- V[!y.remove, ,drop = FALSE]
  V.offset <- V.offset[!y.remove]
  # now do ym
  ym <- ym[!ym.remove, ,drop = FALSE]
  Xm <- X[!ym.remove, ,drop = FALSE]
  Xm.offset <- X.offset[!ym.remove]
  Vm <- Vm[!ym.remove, ,drop = FALSE]

  if(sum(y.remove) > 0){
    y.sites.remove<- which(y.remove)
    y.sites<- sites[-y.sites.remove]
    y.removed.sites<- y.sites.remove
    warning(paste(sum(y.remove)," removal sites have been discarded because of missing data."),call. = FALSE)
  } else {
    y.sites<- sites
    y.removed.sites<- integer(0)
  }

  if(sum(ym.remove) > 0) {
    ym.sites.remove<- which(ym.remove)
    ym.sites<- sites[-ym.sites.remove]
    ym.removed.sites<- ym.sites.remove
    warning(paste(sum(ym.remove),"monitoring sites have been discarded because of missing data."),call. = FALSE)
  } else {
    ym.sites<- sites
    ym.removed.sites<- integer(0)
  }
   removed.sites<- unique(c(y.removed.sites,ym.removed.sites))

  list(y = y, ym=ym, X = Xr, X.offset = Xr.offset, V = V, V.offset = V.offset,
       Xm=Xm, Xm.offset = Xm.offset, Vm=Vm, y.sites = y.sites, ym.sites = ym.sites,
       removed.sites = removed.sites)
}

#---------
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

