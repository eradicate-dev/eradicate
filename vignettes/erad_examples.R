## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(eradicate)
library(sf)
library(terra)

## ---- fig.height=5, fig.width=7-----------------------------------------------
region<- read_sf(system.file("extdata", "shape/san_nic_region.shp", package="eradicate"))
habitat<- rast(system.file("extdata", "san_nic_habitat.tif", package="eradicate"))

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
traps_sf<- st_as_sf(traps, coords=c(1,2), crs=st_crs(region))
traps_buff<- st_buffer(traps_sf, dist=500)
pgrass<- terra::extract(habitat, vect(traps_buff), fun=mean, na.rm=TRUE)
names(pgrass)<- c("id","pgrass")
site.data<- cbind(pgrass, traps)


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
nhat1<- calcN(m1)
nhat2<- calcN(m2)
nhat3<- calcN(m3)

nhat1$Occ
nhat2$Nhat
nhat3$Nhat


## -----------------------------------------------------------------------------
counts<- san_nic_rest$y
stay<- san_nic_rest$stay
cens<- san_nic_rest$cens
A<- san_nic_rest$area # in km2
active<- san_nic_rest$active_hours

emf<- eFrameREST(counts, stay, cens, A, active, siteCovs = site.data)

m4<- REST(~pgrass, data=emf)
summary(m4)

calcN(m4)$Nhat


## -----------------------------------------------------------------------------
rem<- san_nic_rem$rem
ym<- san_nic_rem$ym # detections from additional monitoring
nights<- san_nic_rem$nights
traps<- san_nic_rem$traps

traps_sf<- st_as_sf(traps, coords=c(1,2), crs=st_crs(region))
traps_buff<- st_buffer(traps_sf, dist=500)
pgrass<- terra::extract(habitat, vect(traps_buff), fun=mean, na.rm=TRUE)
names(pgrass)<- c("id","pgrass")
site.data<- cbind(pgrass, traps)

# Poisson abundance 
emf<- eFrameR(rem, siteCovs = site.data)
r1<- remMN(~pgrass, ~1, data=emf)
summary(r1)
nr1<- calcN(r1)
nr1$Nhat
nr1$Nresid

# Generalized Poisson
emf<- eFrameGR(rem, numPrimary=1, siteCovs = site.data)
r2<- remGR(~pgrass, ~1, ~1, data=emf)
summary(r2)
nr2<- calcN(r2)
nr2$Nhat
nr2$Nresid

# Including additional camera monitoring data in ym 
emf<- eFrameGRM(rem, ym, numPrimary=1, siteCovs = site.data)
r3<- remGRM(~pgrass, ~1, ~1, ~1, data=emf)
summary(r3)
nr3<- calcN(r3)
nr3$Nhat
nr3$Nresid


## -----------------------------------------------------------------------------
catch<- apply(rem,2,sum)
effort<- rep(nrow(rem), length(catch))
# extra monitoring/effort data
index<- apply(ym,2,sum)
ieffort<- rep(nrow(ym), length(index))

emf<- eFrameGP(catch, effort, index=index, ieffort=ieffort)
mod<- remGP(emf)
summary(mod)
nce<- calcN(mod)
nce$Nhat
nce$Nresid


