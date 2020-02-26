
#' San_nic_pre
#'
#' Cat monitoring on San Nicolas Island
#'
#' A simulated example dataset of pre-implementation monitoring of a cat population on San
#' Nicolas Island, California.  The data consist of the counts of cats recorded on 55
#' systematically spaced cameras on the island set for 30 nights.
#'
#' @format A list consisting of 2 elements
#'  \describe{
#'  \item{counts}{A matrix of counts with one row for each camera}
#'  \item{traps}{a data.frame of camera coordinates}
#'  }
#'
"san_nic_pre"

#' San_nic_rem
#'
#' Cat removal and monitoring on San Nicolas Island
#'
#' A simulated example dataset of removal sampling and additional camera monitoring
#' on San Nicolas Island, California.  The data consist of the counts of cats removed
#' at each of 100 leg-hold traps spaced randomly on the island over 20 separate removal
#' (primary) periods. For each primary period, traps were set for 10 consecutive nights
#' (secondary period).
#'
#'  The additional monitoring data consist of cat detections at cameras set at 30 of the
#'  leg-hold trap locations for a similar number of primary and secondary periods.
#'
#' @format A list consisting of 2 elements
#'  \describe{
#'  \item{rem}{An \code{M} x \code{J} matrix of counts of individuals removed at each trap
#'  location for each primary period where M is the number of traps and J is the number of
#'  primary periods}
#'  \item{y1}{An \code{N} x \code{J} matrix of detections recorded over the same period as
#'   \code{rem} where \code{N} is the number of cameras and \code{J} is the number of primary periods}
#'  \item{cells}{A vector of length \code{N} giving the indices of each of the the removal traps
#'  where the additional monitoring was undertaken.  Hence the additional monitoring data
#'  \code{y1} was collected at \code{rem[cells,]}}
#'  \item{nights}{The number of secondary periods of sampling for each of the \code{J} primary
#'  periods}
#'  \item{traps}{a data.frame with \code{M} rows of trap locations}
#'  }
#'
"san_nic_rem"

#' rest
#'
#' Example camera monitoring data collected for use in the REST model
#'
#' A simulated example dataset of camera detections collected for use with the
#' REST model of (Nakashima et al. 2018).  The data consist of a the number of individuals
#' detected at each camera as well as the time taken for each individual to cross the
#' camera viewshed (staying time).  It is assumed that detection within the camera
#' viewshed is perfect.
#'
#'
#' @format A list consisting of 2 elements
#'  \describe{
#'  \item{y}{An \code{M} x \code{J} matrix of counts of individuals detected at each camera
#'  where M is the number of cameras and J is the maximum number of days the camera was
#'  operating}
#'  \item{stay}{A vector of the observed staying times, which are the difference between
#'  the exit and entry times for all (or a sample) of individuals}
#'  \item{cens}{A vector indicating whether both entry and exit times were
#'  observed for each of the staying times (\code{cens}=1) or only the entry time
#'  was observed with the exit time unknown (\code{cens}=0).  In the latter case,
#'  the staying time is the maximum time that the individual was observed to be in
#'  front of the camera}
#'  \item{area}{The area of the camera viewshed (km2)}
#'  \item{active_hours}{The number of hours in a single 24 hour period that individuals of
#'  the target species were active for}
#'  \item{traps}{a data.frame with \code{M} rows of trap locations}
#'  }
#'
"rest"


#' @name san_nic_habitat.tif
#' @title Simulated habitat (relative risk) on San Nicolas Island.
#' @description A raster file at 100 meter resolution with simulated habitat values (either
#' \code{0} or \code{1}). Can be resampled for use with any of the models.
#' @examples
#'  library(raster)
#'  habitat<- raster(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))
#'
NULL

#' @name shape
#'
#' @title San Nicholas island.
#' @description ESRI shapefile of San Nicolas Island.
#' @examples
#'  library(sf)
#'  region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
NULL


