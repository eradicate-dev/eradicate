#======================================================
# data objects with similar structure to unmarkedFrame
#======================================================

#' eFrame
#'
#' \code{eFrame} creates an eFrame data object for use with Nmixture
#' or occupancy type models
#'
#' @param y An MxJ matrix of the observed measured data, where M is the
#'    number of sites and J is the maximum number of observations per site.
#' @param SiteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#'
#' @return a \code{eFrame} holding data containing the response and
#'  covariates required for each model
#'
#' @examples
#'  counts<- san_nic_pre$counts
#'  emf <- eFrame(y=counts)
#'  summary(emf)
#'
#' @export
#'
eFrame <- function(y, siteCovs = NULL) {
    if("data.frame" %in% class(y))
        y <- as.matrix(y)
    M <- nrow(y)
    J <- ncol(y)
    if(!is.null(siteCovs))
      if(nrow(siteCovs) != M) stop("siteCov Data does not have same size number of sites as y")
    emf <- list(y = y, siteCovs = siteCovs)
    class(emf) <- c('eFrame', 'list')
    return(emf)
}

#-------------------------------
#' eFrameREST
#'
#' \code{eFrameREST} creates an eFrame data object for use with the REST
#' model applied to camera trap data using the MLE given in:
#' Nakashima,Y., Fukasawa, K., and Samejima H. (2018). J.App.Ecol,55,735-744.
#'
#' @param y An MxJ matrix of the observed number of encounters at a camera over
#' a specified period (usually 24h but could be less - see \code{active_hours})
#' where M is the number of sites and J is the maximum number of observations per site.
#' @param stay A vector of the observed staying times, which are the difference
#'  between the exit and entry times for all (or a sample) of individuals
#' @param cens A vector indicating whether both entry and exit times were
#'  observed for each of the staying times (\code{cens}=1) or only the entry time
#'  was observed with the exit time unknown (\code{cens}=0).  In the latter case,
#'  the staying time is the maximum time that the individual was observed to be in
#'  front of the camera.
#' @param area A scalar indicating The area of the camera viewshed in square meters.
#'  This area is assumed to have a detection probability of 1.0 (effective sampling area).
#' @param active_hours A scalar indicating the number of hours in each 24 hour
#'  period that the species is active
#' @inheritParams eFrame
#'
#' @return a \code{eFrameREST} holding data containing the response and
#'  covariates required for the REST model
#'
#' @examples
#' counts<- rest$y
#' stay<- rest$stay
#' cens<- rest$cens
#' A<- rest$area # in km2
#' active<- rest$active_hours
#'
#' emf<- eFrameREST(counts, stay, cens, A, active, siteCovs = site.data)
#' summary(emf)
#'
#' @export
#'
eFrameREST <- function(y, stay, cens, area, active_hours, siteCovs = NULL) {

  if(length(stay) != length(cens)) stop("length of staying times != length of censor variable")
  if(any(is.na(stay))) stop("Staying times must have no missing values")
  if(!(all(cens %in% c(0,1)))) stop("censoring variable must only contain 0/1 values and no missing")
  R<- ncol(y)
  M<- nrow(y)
  tt<- active_hours * 60 * 60 # effort per day in secs
  effort<- matrix(tt, M, R)
  emf <- eFrame(y, siteCovs)
  emf$area<- area
  emf$effort<- effort
  emf$stay<- stay
  emf$cens<- cens
  class(emf) <- c("eFrameREST",class(emf))
  emf
}

#' eFrameR
#'
#' \code{eFrameR} creates an eFrameR data object for use with removal models.
#' This model is designed to be used with data over a single primary period.
#'
#' @param y An MxJ matrix of the observed removal data, where M is the
#'    number of sites and J is the maximum number of removal (primary)
#'    periods per site. Each primary period can consist of k secondary
#'    periods but this is not used here.
#' @param type sampling method. Currently supports "removal" for removal
#' sampling and "double" for double observer sampling
#' @param SiteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#'
#' @return a \code{eFrameR} holding data containing the response and
#'  covariates required for removal models
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  emf<- eFrameR(rem, type="removal")
#'  summary(emf)
#'
#' @export
#'
eFrameR <- function(y, type, siteCovs = NULL) {
  if(!missing(type)) {
    switch(type,
           removal = piFun <- "removalPiFun",
           double = piFun <- "doublePiFun")
  } else stop("Removal type required")

  emf <- eFrame(y, siteCovs)
  emf$piFun<- piFun
  emf$samplingMethod<- type
  class(emf) <- c("eFrameR",class(emf))
  emf
}

#' eFrameGR
#'
#' \code{eFrameGR} creates an eFrameGR data object for use with generalized removal
#' models using the robust design where sampling occurs over a number of primary and
#' secondary periods.
#'
#'
#' @param y An NxJ matrix of the additional monitoring data, where N is the
#'    number of monitored sites and J is the maximum number of primary periods
#'    per site as for \code{y}.
#' @param numPrimary the number of primary periods. For each primary period, the
#' population is assumed to be closed.
#' @param Z integer indicating the amount of monitoring effort occurring for each
#' monitoring device used during each primary period. Ignored if \code{numPrimary}
#' is greater than 1.
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate.
#' @param primaryCovs A \code{data.frame} of covariates that vary at the
#'    site x primary period level.
#' @return a \code{eFrameGR} holding data containing the response and
#'  covariates required for removal models
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  y1<- san_nic_rem$y1 # detections from additional monitoring
#'  mtraps<- san_nic_rem$cells
#'  nights<- san_nic_rem$nights
#'
#'  emf<-eFrameRM(rem, y1, mtraps, nights, type="removal")
#'  summary(emf)
#'
#' @export
#'
eFrameGR<- function(y, numPrimary, siteCovs = NULL, primaryCovs = NULL, type) {

  if(!missing(type)) {
    switch(type,
           removal = piFun <- "removalPiFun",
           double = piFun <- "doublePiFun")
  } else stop("Removal type required")

  emf <- eFrame(y, siteCovs)
  emf$piFun<- piFun
  emf$samplingMethod<- type
  emf$numPrimary <- numPrimary
  emf$primaryCovs <- covsToDF(primaryCovs, "primaryCovs", numPrimary, nrow(y))
  class(emf) <- c("eFrameGR",class(emf))
  emf
}


#' eFrameGRM
#'
#' \code{eFrameGRM} creates an eFrameGRM data object for use with generalized removal
#' models with additional monitoring using the robust design where sampling occurs
#' over a number of primary and secondary periods.  The additional monitoring
#' usually consists of an index of abundance and hence, inference follows the
#' combined removal + index/removal method (Chen, Pollock & Hoenig (1998) Biometrics,
#' 54, 815-827).
#'
#' @param y An MxJ matrix of the observed removal data, where M is the
#'    number of sites and J is the maximum number of removal (primary)
#'    periods per site. Each primary period can consist of k secondary
#'    periods but this is not used here.
#' @param ym An MxJ matrix of the additional monitoring (index) data.
#' @param numPrimary the number of primary periods. For each primary period, the
#' population is assumed to be closed.
#' @param Z integer indicating the amount of monitoring effort occurring for each
#' monitoring device used during each primary period. Ignored if \code{numPrimary}
#' is greater than 1.
#' @param SiteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param primaryCovs A \code{data.frame} of covariates that vary at the
#'    site x primary period level.
#' @return a \code{eFrameRGP} holding data containing the response and
#'  covariates required for removal models
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  ym<- san_nic_rem$ym # detections from additional monitoring
#'
#'  emf<-eFrameGRM(rem, ym, numPrimary=1, type="removal")
#'  summary(emf)
#'
#' @export
#'
eFrameGRM<- function(y, ym, numPrimary, siteCovs = NULL, primaryCovs = NULL, type) {

  if(!missing(type)) {
    switch(type,
           removal = piFun <- "removalPiFun",
           double = piFun <- "doublePiFun")
  } else stop("Removal type required")

  emf <- eFrame(y, siteCovs)
  emf$ym<- ym
  emf$piFun<- piFun
  emf$samplingMethod<- type
  emf$numPrimary <- numPrimary
  emf$primaryCovs <- covsToDF(primaryCovs, "primaryCovs", numPrimary, nrow(y))
  class(emf) <- c("eFrameGRM",class(emf))
  emf
}

#' eFrameGP
#'
#' \code{eFrameGP} creates an eFrameGP data object for use with the single site
#' removal estimator \code{remGP()}.
#'
#' @param catch A vector of removals for each period.
#' @param effort A vector of removal effort employed during each period (i.e. trapnights).
#' @param index Optional vector of relative abundance indices for each period.
#'
#' @return a \code{eFrameGP} holding data suitable for use in \code{remGP}
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  ym<- san_nic_rem$ym # detections from additional monitoring
#'  catch<- apply(rem,2,sum)
#'  effort<- rep(nrow(rem), length(catch))
#'  index<- apply(ym,2,sum)
#'
#'  emf<-eFrameGP(catch, effort, index)
#'  summary(emf)
#'
#' @export
#'
eFrameGP<- function(catch, effort, index=NULL, ieffort=NULL) {
  if(length(catch) != length(effort))
    stop("catch and effort vectors must be same length")
  x <- data.frame(catch=catch, effort=effort)
  nobs<- nrow(x)
  if(!is.null(index) & length(index) == nobs) {
    x$index<- index
    if(!is.null(ieffort) & length(ieffort) == nobs)
      x$ieffort<- ieffort
    else stop("index supplied with no or faulty effort data")
  }
  emf<- list(counts=x)
  class(emf) <- c("eFrameGP")
  emf
}

############################ EXTRACTORS ##################################

siteCovs<- function(object) return(object$siteCovs)

numSites<- function(object) nrow(object$y)

numY<- function(object) ncol(object$y)

getY<- function(object) object$y

covsToDF <- function(covs, name, obsNum, numSites){
# Convert covs provided as list of matrices/dfs to data frame
  if(!inherits(covs, "list")) return(covs)
  lapply(covs, function(x){
    if(!inherits(x, c("matrix", "data.frame")))
      stop(paste("At least one element of", name, "is not a matrix or data frame."))
    if(ncol(x) != obsNum | nrow(x) != numSites)
      stop(paste("At least one element of", name, "has incorrect number of dimensions."))
  })
  data.frame(lapply(covs, function(x) as.vector(t(x))))
}

################################### PRINT/SUMMARY METHODS ######################
#' print.eFrame
#'
#' \code{print} method for eFrame objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrame} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#' ## uses san_nic_pre
#'  emf <- eFrame(y=counts)
#'  summary(emf)
#'
#' @export
#'
print.eFrame<- function(object,...) {
  #S3 method for eFrame
  cat("eFrame Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of observations per site:",numY(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of observations per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#' summary.eFrame
#'
#' \code{summary} method for eFrame objects.
#'
#' @param object An \code{eFrame} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#' ## uses san_nic_pre
#'  emf <- eFrame(y=counts)
#'  summary(emf)
#'
#' @export
#'
summary.eFrame<- function(object,...) {
  #S3 method for eFrame
  cat("eFrame Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of observations per site:",numY(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of observations per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#-----------------------------------------
#' print.eFrameR
#'
#' \code{print} method for eFrameR objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameR} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#' ## uses san_nic_rem
#'  emf <- eFrameR(y=rem, type="removal")
#'  summary(emf)
#'
#' @export
#'
print.eFrameR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameR Object\n\n")
  cat("sampling method: ",object$samplingMethod,"\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y,na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y,na.rm=TRUE),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#' summary.eFrameR
#'
#' \code{summary} method for eFrameR objects.
#'
#' @param object An \code{eFrameR} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#' ## uses san_nic_rem
#'  emf <- eFrame(y=rem, type="removal")
#'  summary(emf)
#'
#' @export
#'
summary.eFrameR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameR Object\n\n")
  cat("sampling method: ",object$samplingMethod,"\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y,na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y,na.rm=TRUE),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#-----------------------------------------
#' print.eFrameGR
#'
#' \code{print} method for eFrameGR objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameGR} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses san_nic_rem
#'  emf <- eFrameRM(y=rem, y1=y1, cells=cells, Z=nights)
#'  summary(emf)
#'
#' @export
#'
print.eFrameGR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGR Object\n\n")
  cat(nrow(object$y), " removal sites\n\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",numY(object)/object$numPrimary,"\n\n")
  cat("Mean effort per secondary period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y, na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#' summary.eFrameGRM
#'
#' \code{summary} method for eFrameRM objects.
#'
#' @param object An \code{eFrameRM} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses san_nic_rem
#'  emf <- eFrameRM(y=rem, y1=y1, cells=cells, Z=nights)
#'  summary(emf)
#'
#' @export
#'
summary.eFrameGR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGR Object\n\n")
  cat(nrow(object$y), " removal sites\n\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",numY(object)/object$numPrimary,"\n\n")
  cat("Mean effort per secondary period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y, na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}

#-----------------------------------------
#' print.eFrameGRM
#'
#' \code{print} method for eFrameRM objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameRM} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses san_nic_rem
#'  emf <- eFrameRM(y=rem, y1=y1, cells=cells, Z=nights)
#'  summary(emf)
#'
#' @export
#'
print.eFrameGRM<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGRM Object\n\n")
  cat(nrow(object$y), " removal sites\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",numY(object)/object$numPrimary,"\n\n")
  cat("Mean effort per secondary period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y, na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(nrow(object$ym), "Additional monitoring sites\n")
  cat("Tabulation of additional monitoring observations:")
  print(table(object$ym, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}
#' summary.eFrameGRM
#'
#' \code{summary} method for eFrameGRM objects.
#'
#' @param object An \code{eFrameRM} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses san_nic_rem
#'  emf <- eFrameRM(y=rem, y1=y1, cells=cells, Z=nights)
#'  summary(emf)
#'
#' @export
#'
summary.eFrameGRM<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGRM Object\n\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",numY(object)/object$numPrimary,"\n\n")
  cat("Mean effort per secondary period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y,na.rm=TRUE)))
  cat("\n")
  cat("Total removed:",sum(object$y,na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(nrow(object$y1), "monitoring sites\n")
  cat("Tabulation of monitoring observations:")
  print(table(object$y1, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
}

#-------------------------------------
#' print.eFrameREST
#'
#' \code{print} method for eFrameREST objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameREST} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses rest
#'  emf <- eFrameREST(y, stay, cens, area, active_hours)
#'  summary(emf)
#'
#' @export
#'
print.eFrameREST<- function(object,...) {
  #S3 method for eFrameREST
  cat("eFrameREST Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of observations per site:",numY(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of observations per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  cat("\nNumber of Staying times: ",length(object$stay),"\n")
  cens.obs<- table(object$cens)
  cat("\nNumber of censored observations: ",cens.obs[1],"\n")
  cat("\nStaying time summary:\n")
    print(summary(object$stay))

}
#' summary.eFrameREST
#'
#' \code{summary} method for eFrameREST objects.
#'
#' @param object An \code{eFrameREST} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  ## uses rest
#'  emf <- eFrameREST(y, stay, cens, area, active_hours)
#'  summary(emf)
#'
#' @export
#'
summary.eFrameREST<- function(object,...) {
   #S3 method for eFrameREST
  cat("eFrameREST Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of observations per site:",numY(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of observations per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  cat("\nNumber of Staying times: ",length(object$stay),"\n")
  cens.obs<- table(object$cens)
  cat("\nNumber of censored observations: ",cens.obs[1],"\n")
  cat("\nStaying time summary:\n")
    print(summary(object$stay))
}
