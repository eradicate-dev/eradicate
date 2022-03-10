## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(eradicate)
library(sf)
library(terra)

## -----------------------------------------------------------------------------

region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
habitat<- rast(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))


## ----fig.height=5, fig.width=7------------------------------------------------
traps<- san_nic_open$traps
y<- san_nic_open$removals
head(y)

plot(habitat, axes = FALSE)
plot(st_geometry(region), add=TRUE)
points(traps, pch=16)

## -----------------------------------------------------------------------------
traps_sf<- st_as_sf(traps, coords=c(1,2), crs=st_crs(region))
traps_buff<- st_buffer(traps_sf, dist=500)
pgrass<- terra::extract(habitat, vect(traps_buff), fun=mean, na.rm=TRUE)
names(pgrass)<- c("id","pgrass")
site.data<- cbind(pgrass, traps)


## -----------------------------------------------------------------------------
emf<- eFrameMNO(y=y, numPrimary = 8, siteCovs = site.data)
summary(emf)


## -----------------------------------------------------------------------------

fit <- remMNO(~pgrass, ~1, ~1, ~1, K=100, data=emf, dynamics="constant")
summary(fit)


## -----------------------------------------------------------------------------

ests<- calcN(fit)
#Initial population size
ests$Nhat


## -----------------------------------------------------------------------------
# Population size at the start of each session (.season)
ests$Nseason

# Residual population size
ests$Nresid


