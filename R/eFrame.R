#======================================================
# data objects with similar structure to unmarkedFrame
#======================================================

#' eFrame
#'
#' \code{eFrame} creates an eFrame data object for use with Nmixture
#' or occupancy type models
#'
#'
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
#'  emf <- eFrame(y=counts, siteCovs=site.df)
#'  summary(emf)
#'
#' @export
#'
eFrame <- function(y, siteCovs = NULL) {
    if(("data.frame" %in% class(y)) | ("cast_matrix" %in% class(y)))
        y <- as.matrix(y)
    M <- nrow(y)
    J <- ncol(y)
    if(!is.null(siteCovs))
      if(nrow(siteCovs) != M) stop("siteCov Data does not have same size number of sites as y")
    emf <- list(y = y, siteCovs = siteCovs)
    class(emf) <- c('eFrame', 'list')
    return(emf)
}

#' eFrameR
#'
#' \code{eFrameR} creates an eFrameR data object for use with removal models.
#'
#'
#'
#' @param y An MxJ matrix of the observed removal data, where M is the
#'    number of sites and J is the maximum number of removal (primary)
#'    periods per site. Each primary period can consist of k secondary
#'    periods but this is not used here.
#' @param type sampling method
#' @param SiteCovs A \code{data.frame} of covariates that vary at the
#'    site level. This should have M rows and one column per covariate
#'
#' @return a \code{eFrameR} holding data containing the response and
#'  covariates required for removal models
#'
#' @export
#'
eFrameR <- function(y, type, siteCovs = NULL) {
    if(!missing(type)) {
        switch(type,
            removal = piFun <- "removalPiFun",
            double = piFun <- "doublePiFun",
            depDouble = piFun <- "depDoublePiFun")
    } else stop("Removal type required")

    emf <- eFrame(y, siteCovs)
    emf$piFun<- piFun
    emf$samplingMethod<- type
    class(emf) <- c("eFrameR",class(emf))
    emf
}

#' eFrameRM
#'
#' \code{eFrameRM} creates an eFrameRM data object for use with removal models
#' that include additional monitoring at some or all of the sites.
#'
#'
#' @inheritParams eFrameR
#' @param y1 An NxJ matrix of the additional monitoring data, where N is the
#'    number of monitored sites and J is the maximum number of primary periods
#'    per site as for \code{y}.
#' @param cells vector indictaing which of the \code{y} sites were subject
#'  to monitoring
#' @param Z integer indicating the number of secondary periods per
#'  primary period for \code{y1}
#'
#'
#' @return a \code{eFrameRM} holding data containing the response and
#'  covariates required for removal models
#'
#' @export
#'
eFrameRM <- function(y, y1, cells, Z, type, siteCovs = NULL) {

  if(!missing(type)) {
    switch(type,
           removal = piFun <- "removalPiFun",
           double = piFun <- "doublePiFun",
           depDouble = piFun <- "depDoublePiFun")
  } else stop("Removal type required")

  emf <- eFrame(y, siteCovs)
  if(length(cells) != nrow(y1)) stop("number of cells needs to be same as monitored sites")
  emf$y1<- y1
  emf$cells<- cells
  emf$Z<- Z
  emf$piFun <- piFun
  emf$samplingMethod <- type
  class(emf) <- c("eFrameRM",class(emf))
  emf
}

#-------------------------------
#' eFrameREST
#'
#' \code{eFrameREST} creates an eFrame data object for use with the REST
#' model applied to camera trap data using the MLE given in:
#' Nakashima,Y., Fukasawa, K., and Samejima H. (2018). J.App.Ecol,55,735-744.
#'
#' @param y An MxJ matrix of the observed measured data, where M is the
#'    number of sites and J is the maximum number of observations per site.
#' @param stay A vector of the observed staying times, which are the difference
#'  between the exit and entry times for all (or a sample) of individuals
#' @param cens A vector indicating whether both entry and exit times were
#'  observed for each of the staying times (cens=1) or only the entry time
#'  was observed with the exit time unknown (cens=0).  In the latter case,
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
#' @export
#'
eFrameREST <- function(y, stay, area, cens, active_hours, siteCovs = NULL) {

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

############################ EXTRACTORS ##################################

siteCovs<- function(object) return(object$siteCovs)

obsCovs<- function(object) return(object$obsCovs)

numSites<- function(object) nrow(object$y)

numY<- function(object) ncol(object$y)

getY<- function(object) object$y


################################### PRINT/SUMMARY METHODS ######################
#' @rdname print
#' @export
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
#' @rdname summary
#' @export
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
#' @rdname print
#' @export
print.eFrameR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameR Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y)))
  cat("\n")
  cat("Total removed:",sum(y),"\n")
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
#' @rdname summary
#' @export
summary.eFrameR<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameR Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y)))
  cat("\n")
  cat("Total removed:",sum(y),"\n")
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
#' @rdname print
#' @export
print.eFrameRM<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameRM Object\n\n")
  cat(nrow(object$y), " removal sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n")
  cat("Mean replicate counts per period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y)))
  cat("\n")
  cat("Total removed:",sum(y),"\n\n")
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
#' @rdname summary
#' @export
summary.eFrameRM<- function(object,...) {
  #S3 method for eFrame
  cat("eFrameRM Object\n\n")
  cat(nrow(object$y), "sites\n")
  cat("Maximum number of periods per site:",numY(object),"\n")
  cat("Mean replicate counts per period:",mean(object$Z),"\n\n")
  cat("Number of removals per period:","\n")
  print(data.frame(Period=1:numY(object),Removed=colSums(object$y)))
  cat("\n")
  cat("Total removed:",sum(y),"\n\n")
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
#' @rdname print
#' @export
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
#' @rdname summary
#' @export
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
