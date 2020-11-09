## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(eradicate)
library(sf)
library(raster)

## -----------------------------------------------------------------------------

region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
habitat<- raster(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))


## ----fig.height=5, fig.width=7------------------------------------------------
traps<- san_nic_open$traps
y<- san_nic_open$removals
head(y)

plot(habitat, axes = FALSE)
plot(st_geometry(region), add=TRUE)
points(traps, pch=16)

## -----------------------------------------------------------------------------
habvals<- extract(habitat, traps, buffer=500)
pgrass<- sapply(habvals, function(x) mean(x, na.rm=T))
site.data<- cbind(traps, pgrass)


## -----------------------------------------------------------------------------
emf<- eFrameMNO(y=y, numPrimary = 8, siteCovs = site.data)
summary(emf)


## -----------------------------------------------------------------------------

fit <- remMNO(~pgrass, ~1, ~1, ~1, K=100, data=emf, dynamics="constant")
summary(fit)


## -----------------------------------------------------------------------------

calcN(fit)$Nhat


## -----------------------------------------------------------------------------
set.seed(1)
re<- raneffects(fit)
ppd<- postSamples(re, 1000)
Nhat<- round(apply(ppd, c(1,2), mean))

apply(Nhat,2,sum)  # Total population size for each trapping session (= primary period)


