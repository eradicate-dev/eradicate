
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

#' san_nic_rest
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
"san_nic_rest"

#' HogDeer
#'
#' Example camera trapping data set containing camera encounters of hog deer
#' (Axis porcinus).
#'
#' The dataset consists of distance sampling data from camera traps as well as data
#' on residence times in front of the camera (staying time) for each encounter.
#' Distance sampling data is used for estimating the effective sampling area of the
#' camera viewshed using \code{dist_func}. The effective camera viewshed area and staying
#' time data are required for use in the \code{REST} model. Distance bin numbers for
#' images of hog deer were derived from markers placed in the camera viewshed at
#' 2.5, 5.0, 7.5 and 10 m from the camera. The right truncation distance was set at 12.5 m.
#'
#'
#' @format A list consisting of 3 elements
#'  \describe{
#'  \item{sites}{A \code{data.frame} containing information for each camera, including
#'  the first day (\code{sday}) and last day (\code{eday}) of operation for each camera}
#'  \item{encounters}{A \code{data.frame} containing columns \code{day} giving the day
#'  of encounter, \code{distance}, giving the distance  bin numbers. Distance bins are
#' numbered with 1 indicating an image in the 1st bin (i.e. 0 - 2.5 m) and 5
#' indicating an image in the 5th bin (i.e. 10 - 12.5m). These are converted to distances
#' (i.e. distance to the midpoint of each bin in meters) automatically by \code{eFrameDS}.
#'   Group size of the encounter is given in \code{count} and \code{cam} gives the camera ID.
#'  \code{stay} gives the estimate of residence time (staying time in seconds) and
#'   \code{censor} is an variable indicating a fully observed (\code{censor} = 1) or right
#'   censored (\code{censor} = 0) staying time.}
#'  \item{cutpoints}{A vector of the distance categories (bins) used.}
#'  }
#'
"HogDeer"

#' @name san_nic_habitat
#' @title Simulated habitat (relative risk) on San Nicolas Island.
#' @description A raster file at 100 meter resolution with simulated habitat values (either
#' \code{0} or \code{1}). Can be resampled for use with any of the models.
#' @examples
#'  library(raster)
#'  habitat<- raster(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))
#'
NULL

#' @name san_nic_region
#' @title San Nicholas island.
#' @description ESRI shapefile of San Nicolas Island.
#' @examples
#'  library(sf)
#'  region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
NULL


