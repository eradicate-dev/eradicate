% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eFrame.R
\name{eFrameGRM}
\alias{eFrameGRM}
\title{eFrameGRM}
\usage{
eFrameGRM(y, ym, siteCovs = NULL, obsCovs = NULL)
}
\arguments{
\item{y}{An MxJ matrix of the observed removal data, where M is the
number of sites and J is the maximum number of removal (primary)
periods per site. Each primary period can consist of k secondary
periods but this is not used here.}

\item{ym}{An MxJ matrix of the additional monitoring (index) data.}

\item{obsCovs}{A list of matrices or data.frames of variables varying within sites.
Each matrix or data.frame must be of dimension MxJ.}

\item{SiteCovs}{A \code{data.frame} of covariates that vary at the
site level. This should have M rows and one column per covariate}
}
\value{
a \code{eFrameGRM} holding data containing the response and
 covariates required for removal models
}
\description{
\code{eFrameGRM} creates an eFrameGRM data object for use with generalized removal
models with additional monitoring using the robust design where sampling occurs
over a number of primary and secondary periods.  The additional monitoring
usually consists of an index of abundance and hence, inference follows the
combined removal + index/removal method (Chen, Pollock & Hoenig (1998) Biometrics,
54, 815-827).
}
\examples{
 rem<- san_nic_rem$rem
 ym<- san_nic_rem$ym # detections from additional monitoring

 emf<-eFrameGRM(rem, ym, type="removal")
 summary(emf)

}
