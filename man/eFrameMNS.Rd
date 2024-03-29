% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eFrame.R
\name{eFrameMNS}
\alias{eFrameMNS}
\title{eFrameMNS}
\usage{
eFrameMNS(rem, siteCovs = NULL, obsCovs = NULL, delta = NULL)
}
\arguments{
\item{rem}{A \code{data.frame} of the observed removal data for each site
in rows and secondary periods in columns, indexed by session (primary period).
The data.frame \code{rem} must contain a column \code{session} with at least two unique values.
The dimensions of \code{rem} are thus \code{MT} rows and \code{J+1} columns
where \code{T} is the number of sessions.}

\item{siteCovs}{A \code{data.frame} of covariates that vary at the
site level. This should have M rows and one column per covariate}

\item{obsCovs}{A list of matrices or data.frames of variables varying within sites.
Each matrix or data.frame must be of dimension \code{MT x J}.}

\item{delta}{A vector with elements giving the time units between primary periods for each site
beginning with 1 for the first primary period. A default value of 1 is used to
indicate equal time intervals between primary periods.}
}
\value{
a \code{eFrameMNS} holding stacked data for each primary period.
}
\description{
\code{eFrameMNS} creates an eFrameMNS 'stacked' data object for use with closed population
multinomial removal models using the robust design where sampling occurs over
a number of primary and secondary periods. Data for each primary period is 'stacked'
into rows with an indicator variable added to identify each primary period. The data can
then be analysed using closed population removal models to estimate the trend in abundance
between primary periods.
}
\examples{
 rem<- san_nic_open$removals

 emf<-eFrameMNS(rem)
 summary(emf)

}
