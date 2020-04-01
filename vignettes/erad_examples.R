## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(eradicate)
library(sf)
library(raster)

## ---- fig.height=5, fig.width=7-----------------------------------------------
region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
habitat<- raster(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))

plot(habitat, axes=F)
plot(st_geometry(region), add=T)

## ----fig.height=5, fig.width=7------------------------------------------------
traps<- san_nic_pre$traps
counts<- san_nic_pre$counts
dim(counts)

plot(habitat, axes = FALSE)
plot(st_geometry(region), add=TRUE)
points(traps, pch=16)

## -----------------------------------------------------------------------------
habvals<- extract(habitat, traps, buffer=500)
pgrass<- sapply(habvals, function(x) mean(x, na.rm=T))
site.data<- cbind(traps, pgrass)


## -----------------------------------------------------------------------------
emf<- eFrame(counts, siteCovs = site.data)
summary(emf)

# Occupancy model
m1<- occuM(~pgrass, ~1, data=emf)
summary(m1)

# Royle/Nichols model
m2<- occuRN(~pgrass, ~1, K=50, data=emf)
summary(m2)

m3<- nmix(~pgrass, ~1, K=50, data=emf)
summary(m3)

## -----------------------------------------------------------------------------
calcN(m1)
calcN(m2)
calcN(m3)

## -----------------------------------------------------------------------------
counts<- san_nic_rest$y
stay<- san_nic_rest$stay
cens<- san_nic_rest$cens
A<- san_nic_rest$area # in km2
active<- san_nic_rest$active_hours

emf<- eFrameREST(counts, stay, cens, A, active, siteCovs = site.data)

m4<- REST(~pgrass, data=emf)
summary(m4)

calcN(m4)


## -----------------------------------------------------------------------------
rem<- san_nic_rem$rem
ym<- san_nic_rem$ym # detections from additional monitoring
nights<- san_nic_rem$nights
traps<- san_nic_rem$traps

habvals<- raster::extract(habitat, traps, buffer=500)
pgrass<- sapply(habvals, function(x) mean(x, na.rm=T))
site.data<- cbind(traps, pgrass)

# Poisson abundance 
emf<- eFrameR(rem, type="removal", siteCovs = site.data)
r1<- remPois(~pgrass, ~1, data=emf)
summary(r1)
calcN(r1)

# Generalized Poisson
emf<- eFrameGR(rem, numPrimary=1, type="removal", siteCovs = site.data)
r2<- remGR(~pgrass, ~1, ~1, data=emf)
summary(r2)
calcN(r2)

emf<- eFrameGRM(rem, ym, numPrimary=1, type="removal", siteCovs = site.data)
r3<- remGRM(~pgrass, ~1, ~1, ~1, data=emf)
summary(r3)
calcN(r3)

## -----------------------------------------------------------------------------
catch<- apply(rem,2,sum)
effort<- rep(nrow(rem), length(catch))
# extra monitoring/effort data
index<- apply(ym,2,sum)
ieffort<- rep(nrow(ym), length(index))

emf<- eFrameGP(catch, effort, index, ieffort)
mod<- remGP(emf)
summary(mod)
calcN(mod)


