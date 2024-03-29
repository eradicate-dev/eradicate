% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_data.R
\docType{data}
\name{san_nic_open}
\alias{san_nic_open}
\title{San_nic_open}
\format{
A list consisting of 2 elements
 \describe{
 \item{removal}{A \code{MT} x \code{J} data.frame of counts of individuals removed at each trap
 location for each primary period where `M` is the number of traps, `T' is the number of primary
 periods and `J` is the number of secondary periods (nights). The data.frame must have an additional
 column named \code{session} indicating the primary period for each observation. The remaining columns
 contain the number removed during each of the \code{J} secondary periods}
 \item{traps}{a data.frame with \code{M} rows of trap locations}
 }
}
\usage{
san_nic_open
}
\description{
Simulated rabbit removal on San Nicolas Island
}
\details{
A simulated example dataset of removal sampling of rabbits from an open population
on San Nicolas Island, California.  The data consist of the counts of rabbits removed
from each of 100 traps spaced randomly on the island over 8 separate removal
(primary) periods. For each primary period, traps were set for 4 consecutive nights
(secondary period).  Primary period trapping was undertaken every 3 months over a two year
period and hence, it would be unreasonable to expect the population to be closed. Hence,
these data require analysis with an open population removal model (i.e. `remMNO()`).
}
\keyword{datasets}
