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
emf<- eFrameMNS(y, siteCovs = site.data)
summary(emf)


## -----------------------------------------------------------------------------

fit1 <- remMNS(~pgrass + .season, ~1, data=emf)
fit2 <- remMNS(~pgrass + .season, ~.trend, data=emf)
fit3 <- remMNS(~pgrass + .trend, ~1, data=emf)
fit4 <- remMNS(~pgrass + .trend, ~.trend, data=emf)

AIC(fit1)
AIC(fit2)
AIC(fit3)
AIC(fit4)


## -----------------------------------------------------------------------------

summary(fit1)

summary(fit3)


## -----------------------------------------------------------------------------

ests1<- calcN(fit1)
#Initial population size
ests1$Nhat

ests3<- calcN(fit3)
ests3$Nhat

## -----------------------------------------------------------------------------
# Residual population size
ests1$Nresid

ests3$Nresid


## -----------------------------------------------------------------------------

idx<- san_nic_open$index
head(idx)


## -----------------------------------------------------------------------------

emf<- eFrameGRMS(y, idx, siteCovs = site.data)
mod<- remGRMS(~pgrass + .season, ~1, ~1, data=emf)
summary(mod)


## -----------------------------------------------------------------------------

ests<- calcN(mod)

ests$Nhat

ests$Nresid

## -----------------------------------------------------------------------------

emf <- eFrameMS(y, siteCovs=site.data)

# order of terms is lamformula (occupancy), gamformula (colonisation), epsformula (extinction)
# and detformula (detection)

occ1 <- occMS(~pgrass, ~1, ~1, ~1, emf)  
occ2 <- occMS(~pgrass, ~.season, ~1, ~1, emf)
occ3 <- occMS(~pgrass, ~1, ~.season, ~1, emf)
occ4 <- occMS(~pgrass, ~.season, ~.season, ~1, emf)

AIC(occ1)
AIC(occ2)
AIC(occ3)
AIC(occ4)


## -----------------------------------------------------------------------------

summary(occ4)


## -----------------------------------------------------------------------------
summary(occ3)


## -----------------------------------------------------------------------------
ests<- calcN(occ3)

ests$Nhat


