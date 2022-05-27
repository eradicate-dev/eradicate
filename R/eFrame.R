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
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
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
eFrame <- function(y, siteCovs = NULL, obsCovs = NULL) {
    if("data.frame" %in% class(y))
        y <- as.matrix(y)
    M <- nrow(y)
    R <- ncol(y)
    if(!is.null(siteCovs))
      if(nrow(siteCovs) != M) stop("siteCov Data does not have same size number of sites as y")
    if(!is.null(obsCovs)) {
      obsCovs <- covsToDF(obsCovs, "obsCovs", R, M)
    }
    emf <- list(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
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
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
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
eFrameR <- function(y, siteCovs = NULL, obsCovs = NULL) {
  emf <- eFrame(y, siteCovs, obsCovs)
  emf$piFun<- "removalPiFun"
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
#' @param y An MxJ matrix of the observed removal data, where M is the
#'    number of sites and J is the maximum number of removal (primary)
#'    periods per site. Each primary period can consist of k secondary
#'    periods but this is not used here.
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate.
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
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
eFrameGR<- function(y, siteCovs = NULL, obsCovs = NULL) {
  emf <- eFrame(y, siteCovs, obsCovs)
  emf$piFun<- "removalPiFun"
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
#' @param SiteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
#' @return a \code{eFrameGRM} holding data containing the response and
#'  covariates required for removal models
#'
#' @examples
#'  rem<- san_nic_rem$rem
#'  ym<- san_nic_rem$ym # detections from additional monitoring
#'
#'  emf<-eFrameGRM(rem, ym, type="removal")
#'  summary(emf)
#'
#' @export
#'
eFrameGRM<- function(y, ym, siteCovs = NULL, obsCovs = NULL) {
  emf <- eFrame(y, siteCovs, obsCovs)
  if((ncol(ym) != ncol(y)) | (nrow(ym) != nrow(y))) stop("ym has different dimensions to y")
  emf$ym<- ym
  emf$piFun<- "removalPiFun"
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
#' @param ieffort Optional vector of effort employed for the relative index
#' @param session Optional vector indicating multiple removal 'sessions'. The population for
#' a particular session is considered to be closed and removals for each session are
#' analysed independently. Analysis of catch/effort data for multiple sessions returns a list
#' with length equal to the number of sessions.
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
eFrameGP<- function(catch, effort, session=NULL, index=NULL, ieffort=NULL) {
  if(length(catch) != length(effort))
    stop("catch and effort vectors must be same length")
  if(is.null(session)) session<- rep(1, length(catch))
  if(!is.null(session) & (length(session) != length(catch)))
    stop("session and catch vectors must be same length")
  x <- data.frame(catch=catch, effort=effort, session=session)
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

#' eFrameMS
#'
#' \code{eFrameMS} creates an eFrameMS data object for use with multi-season
#' occupancy models \code{occMS} where sampling occurs over a number of
#' primary and secondary periods.
#'
#' @param df A \code{data.frame} of the observed detection/non-detection data for each site
#' in rows and secondary periods in columns, indexed by session (primary period). values
#' greater than 0 will be truncated to 1. The \code{df} must contain a column
#'  \code{session} with at least two unique values.
#' @param obsPerSeason scalar of the maximum number of observations (replicates)
#' per season.
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
#' @return a \code{eFrameMS} holding data containing the response and
#'  covariates required for \code{occuMS}
#'
#' @examples
#' ### NEED TO SUPPLY GOOD EXAMPLE DATA
#'  rem<- san_nic_rem
#'  #'
#'  emf<-eFrameMS(y=rem, obsPerSeason=XXX, type="removal")
#'  summary(emf)
#'
#' @export
#'
eFrameMS<- function(df, siteCovs = NULL, obsCovs =  NULL) {
  yl<- unstack.data(df)
  M <- yl$M; T <- yl$T; J <- yl$J
  y<- yl$y
  y <- truncateToBinary(y)
  emf <- eFrame(y, siteCovs, obsCovs)
  emf$numPrimary <- T
  class(emf) <- c("eFrameMS",class(emf))
  emf
}

#' eFrameMNS
#'
#' \code{eFrameMNS} creates an eFrameMNS 'stacked' data object for use with closed population
#' multinomial removal models using the robust design where sampling occurs over
#' a number of primary and secondary periods. Data for each primary period is 'stacked'
#' into rows with an indicator variable added to identify each primary period. The data can
#' then be analysed using closed population removal models to estimate the trend in abundance
#' between primary periods.
#'
#' @param rem A \code{data.frame} of the observed removal data for each site
#' in rows and secondary periods in columns, indexed by session (primary period).
#' The data.frame \code{rem} must contain a column \code{session} with at least two unique values.
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
#' @param delta A vector with elements giving the time units between primary periods for each site
#' beginning with 1 for the first primary period. A default value of 1 is used to
#' indicate equal time intervals between primary periods.
#' @return a \code{eFrameMNS} holding stacked data for each primary period.
#'
#' @examples
#'  rem<- san_nic_open$removals
#'
#'  emf<-eFrameMNS(rem)
#'  summary(emf)
#'
#' @export
#'
eFrameMNS<- function(rem, siteCovs = NULL, obsCovs = NULL, delta = NULL) {

  rd<- get.dims.data(rem)
  M <- rd$M
  T <- rd$T
  J <- rd$J
  y<- rd$y

  if (J %% 1 != 0) stop("Unequal number of secondary periods")
  if (J < 2) stop("less than 2 primary periods present")
  if(missing(delta))
    delta <- rep(1.0, T)
  if(!is.vector(delta) || length(delta) != T)
    stop("delta should be a vector of length T")
  if(any(delta < 0, na.rm=TRUE))
    stop("Negative delta values are not allowed.")
  if(any(is.na(delta)))
    stop("Missing values are not allowed in delta.")
  if(!is.null(siteCovs))
    if(nrow(siteCovs) != M) stop("siteCov Data does not have same size number of sites as y")
  if(!is.null(obsCovs)) {
    obsCovs <- covsToDF(obsCovs, "obsCovs", J, M)
  }
  ya <- array(y, c(M, J, T))
  num.removed <- apply(ya, 3, sum, na.rm=TRUE)
  emf <- list(y=y, siteCovs=siteCovs, obsCovs=obsCovs)
  emf$ya<- ya
  emf$piFun<- "removalPiFun"
  emf$numSites<- M
  emf$numPrimary <- T
  emf$numSecondary<- J
  emf$delta<- cumsum(delta)
  emf$num.removed <- num.removed
  class(emf) <- c('eFrameMNS','eFrame','list')
  emf
}

#' eFrameGRMS
#'
#' \code{eFrameGRMS} creates an eFrameGRMS 'stacked' data object for use with closed population
#' generalized removal models + index data using the robust design where sampling occurs over
#' a number of primary and secondary periods. Data for each primary period is 'stacked'
#' into rows with an indicator variable added to identify each primary period. The data can
#' then be analysed using closed population removal models to estimate the trend in abundance
#' between primary periods.
#'
#' @param rem A \code{data.frame} of the observed removal data for each site
#' in rows and secondary periods in columns, indexed by session (primary period).
#' The \code{data.frame} must contain a column \code{session} with at least two unique values.
#' @param idx A \code{data.frame} of the observed index data for each site
#' in rows and secondary periods in columns, indexed by session (primary period).
#' The \code{data.frame} must contain a column \code{session} with at least two unique values.
#' @param siteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#' @param obsCovs A list of matrices or data.frames of variables varying within sites.
#' Each matrix or data.frame must be of dimension MxJ.
#' @param delta A vector with elements giving the time units between primary periods for each site
#' beginning with 1 for the first primary period. A default value of 1 is used to
#' indicate equal time intervals between primary periods.
#' @return a \code{eFrameGRMS} holding stacked data for each primary period.
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  ym<- san_nic_open$index # detections from additional monitoring
#'
#'  emf<-eFrameGRMS(rem, index)
#'  summary(emf)
#'
#' @export
#'
eFrameGRMS<- function(rem, idx, siteCovs = NULL, obsCovs = NULL, delta = NULL) {

  rd<- get.dims.data(rem)
  M <- rd$M
  T <- rd$T
  J <- rd$J
  y<- rd$y
  ri<- get.dims.data(idx)
  ym<- ri$y

  if(M != ri$M) stop("index data has different number of sites to removal data")
  if(T != ri$T) stop("index data has different number of sessions to removal data")
  if(J != ri$J) stop("index data has different number of secondary periods to removal data")
  if (J %% 1 != 0) stop("Unequal number of secondary periods")
  if (J < 2) stop("less than 2 primary periods present")
  if(missing(delta))
    delta <- rep(1.0, T)
  if(!is.vector(delta) || length(delta) != T)
    stop("delta should be a vector of length T")
  if(any(delta < 0, na.rm=TRUE))
    stop("Negative delta values are not allowed.")
  if(any(is.na(delta)))
    stop("Missing values are not allowed in delta.")
  ya <- array(y, c(M, J, T))
  num.removed <- apply(ya, 3, sum, na.rm=TRUE)
  obsCovs <- covsToDF(obsCovs, "obsCovs", J*T, M)
  emf <- list(y=y, siteCovs=siteCovs, obsCovs=obsCovs)
  emf$ym<- ym
  emf$ya<- ya
  emf$piFun<- "removalPiFun"
  emf$numSites<- M
  emf$numPrimary <- T
  emf$numSecondary<- J
  emf$delta<- cumsum(delta)
  emf$num.removed <- num.removed
  class(emf) <- c("eFrameGRMS",'eFrame','list')
  emf
}

#' eFrameDS
#'
#' \code{eFrameDS} creates an eFrame data object for data collected by
#' distance sampling, specifically from remote cameras.  Distance sampling
#' is used to estimate the area of a camera viewshed by fitting a detection
#' function and estimating the effective sampled area.
#'
#' @param distance A vector of distances (or distance bin number) for each
#' detected individual. Actual distances are the distance (m) to the midpoint
#' of each distance bin.
#' @param size the number of individuals recorded for each distance measurement
#' (i.e. groupsize).
#' @param siteID A vector indicating the camera ID for each distance measurement.
#' @param cutpoints vector of bin cutpoints indicating the distance to the end of
#' each bin. cutpoints should begin at zero and end with w.
#' @param w Truncation distance or maximum distance from camera that will be considered
#' in the analysis. All distances or bins further than this will be discarded.
#' @param bin_nums A logical indicating whether \code{distance} records
#' distances or bin numbers.  Bin numbers are assumed to be numbered from 1.
#' @return a \code{eFrameDS} holding data containing the data suitable for estimating a
#' camera detection function
#'
#' @examples
#'  ddata<- HogDeer$ddata
#'  cutpoints<- HogDeer$cutpoints
#'  emf<- eFrameDS(ddata$distance, ddata$size, ddata$camID, cutpoints, w=12.5, bin_num=TRUE)
#'
#' @export
#'
eFrameDS <- function(distance, size, siteID, cutpoints, w, bin_nums=FALSE) {
  M<- length(distance)
  if(any(is.na(distance))) stop("distance vector should not have missing values")
  if(bin_nums) distance<- convert_bin_to_dist(distance, cutpoints)
  if((cutpoints[1] != 0) & (cutpoints[M] != w))
      stop("cutpoints should start at 0 and end at w")
  if(length(size) != M) stop("size vector should be same length as distance")
  if(length(siteID) != M) stop("siteID vector should be same length as distance")
  ddata <- data.frame(distance=distance,size=size,siteID=siteID)
  emf<- list(ddata=ddata, cutpoints=cutpoints, w=w)
  class(emf) <- c('eFrameDS', 'list')
  return(emf)
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  cat("Maximum number of secondary periods:",numY(object),"\n\n")
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  cat("Maximum number of secondary periods:",numY(object),"\n\n")
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  cat("Maximum number of secondary periods:",numY(object),"\n\n")
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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
  cat("Maximum number of secondary periods:",numY(object),"\n\n")
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
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
  }
}

#-----------------------------------------
#' print.eFrameMNS
#'
#' \code{print} method for eFrameMNS objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameMNS} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  emf <- eFrameMNS(rem)
#'  print(emf)
#'
#' @export
#'
print.eFrameMNS<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameMNS Object\n\n")
  cat(object$numSites, " removal sites\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",object$numSecondary,"\n\n")
  cat("Number of removals per session:","\n")
  print(data.frame(Session=1:object$numPrimary,Removed=object$num.removed))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(object$ya, 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(object$numSites, "Additional monitoring sites\n")
  cat("Tabulation of additional monitoring observations:")
  print(table(object$ym, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
  }
}
#' summary.eFrameMNS
#'
#' \code{summary} method for eFrameMNS objects.
#'
#' @param object An \code{eFrameMNS} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  emf <- eFrameMNS(rem)
#'  summary(emf)
#'
#' @export
#'
summary.eFrameMNS<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameMNS Object\n\n")
  cat(object$numSites, " removal sites\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",object$numSecondary,"\n\n")
  cat("Number of removals per session:","\n")
  print(data.frame(Session=1:object$numPrimary,Removed=object$num.removed))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(object$ya, 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(object$numSites, "Additional monitoring sites\n")
  cat("Tabulation of additional monitoring observations:")
  print(table(object$ym, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
  }
}

#-----------------------------------------
#' print.eFrameGRMS
#'
#' \code{print} method for eFrameGRMS objects. Basically the same as \code{summary.eFrame}
#'
#' @param object An \code{eFrameGRMS} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  idx<- san_nic_open$index
#'  emf <- eFrameGRMS(rem, idx)
#'  print(emf)
#'
#' @export
#'
print.eFrameGRMS<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGRMS Object\n\n")
  cat(object$numSites, " removal sites\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",object$numSecondary,"\n\n")
  cat("Number of removals per session:","\n")
  print(data.frame(Session=1:object$numPrimary,Removed=object$num.removed))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(object$ya, 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(object$numSites, "Additional monitoring sites\n")
  cat("Tabulation of additional monitoring observations:")
  print(table(object$ym, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
  }
}
#' summary.eFrameGRMS
#'
#' \code{summary} method for eFrameGRMS objects.
#'
#' @param object An \code{eFrameGRMS} object.
#'
#' @return a \code{list} containing various summaries of the data
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  idx<- san_nic_open$index
#'  emf <- eFrameGRMS(rem, idx)
#'  summary(emf)
#'
#' @export
#'
summary.eFrameGRMS<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameGRMS Object\n\n")
  cat(object$numSites, " removal sites\n")
  cat("Number of primary periods:",object$numPrimary,"\n\n")
  cat("Maximum number of secondary periods:",object$numSecondary,"\n\n")
  cat("Number of removals per session:","\n")
  print(data.frame(Session=1:object$numPrimary,Removed=object$num.removed))
  cat("\n")
  cat("Total removed:",sum(object$y, na.rm=TRUE),"\n\n")
  cat("Sites with at least one detection:",
      sum(apply(object$ya, 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  cat("Tabulation of y observations:")
  print(table(object$y, exclude=NULL))
  cat("\n")
  cat(object$numSites, "Additional monitoring sites\n")
  cat("Tabulation of additional monitoring observations:")
  print(table(object$ym, exclude=NULL))
  if(!is.null(object$siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object$siteCovs))
  }
  if(!is.null(object$obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object$obsCovs))
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


############################ EXTRACTORS ##################################

siteCovs<- function(object) return(object$siteCovs)

obsCovs<- function(object, matrices=FALSE) {
  M<- numSites(object)
  R<- numY(object)
  if(matrices) {
    value <- list()
    for(i in seq(length=length(object$obsCovs))){
      value[[i]] <- matrix(object$obsCovs[,i], M, R, byrow = TRUE)
    }
    names(value) <- names(object$obsCovs)
  } else {
    value <- object$obsCovs
  }
  return(value)
}

numSites<- function(object) {
  if(!is.null(object$numSites)) return(object$numSites)
  else return(nrow(object$y))
}

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

#-----------------------------------------
# bracket methods
#-----------------------------------------

#' [.eFrame
#'
#' @description Site extractor methods for eFrame objects.
#'
#' @param emf A eFrame object.
#'
#' @export
#'
`[.eFrame` <- function(x, i) {
  M <- numSites(x)
  y <- getY(x)
  J <- ncol(y)

  if(length(i) == 0) return(x)
  if(any(i < 0) && any(i > 0))
    stop("i must be all positive or all negative indices.")
  if(all(i < 0)) { # if i is negative, then convert to positive
    i <- (1:M)[i]
  }
  y <- getY(x)[i,]
  if (length(i) == 1) {
    y <- t(y)
  }
  siteCovs <- siteCovs(x)
  obsCovs <- obsCovs(x)
  if (!is.null(siteCovs)) {
    siteCovs <- siteCovs(x)[i, , drop = FALSE]
  }
  if (!is.null(obsCovs)) {
    .site <- rep(1:M, each = J)
    obsCovs <- obsCovs[which(.site %in% i),]
  }
  emf <- x
  emf$y <- y
  emf$siteCovs <- siteCovs
  emf$obsCovs <- obsCovs
  emf
}

#----------
#' [.eFrameMS
#'
#' @description Site extractor methods for eFrame objects.
#'
#' @param emf A eFrame object.
#'
#' @export
#'
`[.eFrameMS` <- function(x, i) {
  M <- numSites(x)
  J <- x$obsPerSeason
  y <- getY(x)
  T <- ncol(y) / J

  if(length(i) == 0) return(x)
  if(any(i < 0) && any(i > 0))
    stop("i must be all positive or all negative indices.")
  if(all(i < 0)) { # if i is negative, then convert to positive
    i <- (1:M)[i]
  }
  y <- getY(x)[i,]
  if (length(i) == 1) {
    y <- t(y)
  }
  siteCovs <- siteCovs(x)
  obsCovs <- obsCovs(x)

  if (!is.null(siteCovs)) {
    siteCovs <- siteCovs(x)[i, , drop = FALSE]
  }
  if (!is.null(obsCovs)) {
    .site <- rep(1:M, each = J)
    obsCovs <- obsCovs[which(.site %in% i),]
  }
  emf <- x
  emf$y <- y
  emf$siteCovs <- siteCovs
  emf$obsCovs <- obsCovs
  emf
}

