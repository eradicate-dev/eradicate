## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(eradicate)
library(sf)
library(raster)


## ---- fig.height=6, fig.width=8-----------------------------------------------
region<- read_sf(system.file("extdata", "shape/gipps_lakes.shp", package="eradicate"))
habitat<- raster(system.file("extdata", "hog_habitat.tif", package="eradicate"))

sites<- HogDeer$sites
encounters<- HogDeer$encounters


plot(habitat, axes=F)
plot(st_geometry(region), add=TRUE)
points(sites$Easting, sites$Northing, pch=16, col="red", cex=1.2)


## -----------------------------------------------------------------------------

head(sites)

head(encounters)


## ----  fig.height=5, fig.width=7----------------------------------------------
cutpoints<- HogDeer$cutpoints
w<- cutpoints[length(cutpoints)] # truncation distance
emf<- with(encounters, {
  eFrameDS(distance = dist_bin, size = count, siteID = cam, cutpoints = cutpoints, w=w, bin_nums=TRUE)
})

hist(emf$ddata$distance, breaks=cutpoints, main="Distance categories of hog deer encounters",
     xlab="Distance (m)")


## -----------------------------------------------------------------------------
dfuncs<- dist_func(emf, theta = 40)
dfuncs$fit.table

esa<- dfuncs$fit.table$esa[1] # Select esa estimate for hn based on AIC

## ---- fig.height=5, fig.width=7-----------------------------------------------
plot(dfuncs$fits$hn, main="Fitted Detection function", xlab="Distance (m)",
     showpoints=TRUE, lwd=2, xlim=c(0, w))


## -----------------------------------------------------------------------------
# First construct the matrix of encounters having dimensions M x J where M is the number of sites
# and J is the maximum days operation 
y<- make_encounters(sites, encounters)
stay<- encounters$stay
censor<- encounters$censor
active_hours<- 24  # Assume hog deer active at any time of day
site.data<- subset(sites, select=EVC)
area<- esa/10^6  # Effective sampling area in km^2


emf<- eFrameREST(y=y, stay=stay, cens=censor, area=area, active_hours = 24, siteCovs = site.data)

mod<- REST(~EVC, data=emf)
summary(mod)

## -----------------------------------------------------------------------------
newdat<- data.frame(EVC = factor(1:5, labels = c("coastal_scrub","heathlands",
                                                 "forest_woodlands","herb_woodlands","wetlands")))

calcN(mod, newdata=newdat)$cellpreds


## ---- fig.height=6, fig.width=8-----------------------------------------------
#extract habitat EVC values from raster layer

newdat<- as.data.frame(habitat)
names(newdat)<- "EVC"
newdat$EVC<- factor(newdat$EVC) # convert to factor
inds<- which(!is.na(newdat$EVC)) # indicator for all cells in habitat containg EVC values

preds<- calcN(mod, newdata=newdat, off.set = 0.01)
preds$Nhat  # Total abundance

# make copy of habitat layer
habs<- habitat
habs[inds]<- preds$cellpreds$N

plot(habs, col=rev(heat.colors(10)),legend.args = list(text = 'Density/ha',cex=0.8))
plot(st_geometry(region), add=TRUE)

