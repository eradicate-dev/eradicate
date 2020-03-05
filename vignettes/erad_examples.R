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
counts<- rest$y
stay<- rest$stay
cens<- rest$cens
A<- rest$area # in km2
active<- rest$active_hours

emf<- eFrameREST(counts, stay, cens, A, active, siteCovs = site.data)

m4<- REST(~pgrass, data=emf)
summary(m4)

calcN(m4)


## -----------------------------------------------------------------------------
rem<- san_nic_rem$rem
y1<- san_nic_rem$y1 # detections from additional monitoring
mtraps<- san_nic_rem$cells
nights<- san_nic_rem$nights
traps<- san_nic_rem$traps

habvals<- raster::extract(habitat, traps, buffer=500)
pgrass<- sapply(habvals, function(x) mean(x, na.rm=T))
site.data<- cbind(traps, pgrass)


emf<-eFrameRM(rem, y1, mtraps, nights, type="removal", siteCovs = site.data)
r1<- remPoisM(~ pgrass,  ~ 1, data=emf, K=50)
summary(r1)
calcN(r1)


emf<- eFrameR(rem, type="removal", siteCovs = site.data)
r2<- remPois(~pgrass, ~1, data=emf)
summary(r2)
calcN(r2)


