## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(sf)
library(terra)
library(rTARDIS)

## -----------------------------------------------------------------------------
# build a global raster with a constant cell value
rst <- rast(nrow = 180, ncol = 360, vals = 1)

# convert to a geoglist, then a tardis object
rtd <- build_tardis(rast_to_geoglist(rst))

## -----------------------------------------------------------------------------
# calculate cost surface and plot
cs <- cost_surface(rtd)
plot(cs)

## -----------------------------------------------------------------------------
# we supply our point as a two column, lon-lat matrix
orig <- point_check(rtd, cbind(0, 0))

# create a cumulative cost surface from this origin point
cc <- cumulative_cost(rtd, origin = orig)

# plot the surface
plot(cc)

## -----------------------------------------------------------------------------
# load and mask data
gal <- galapagos()[[5]]
gal_l <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# build the geoglist and tardis
rst <- rast_to_geoglist(gal, gal_l)
rtd <- build_tardis(rst)

## -----------------------------------------------------------------------------
# resolve our land point and get our surface plots
pt <- point_check(rtd, rbind(c(-90.75, 0.58)))
cc <- cumulative_cost(rtd, origin = pt)
cs <- cost_surface(rtd)
plot(cs)
plot(cc)

## -----------------------------------------------------------------------------
# mst
rst <- link_islands(rst, klink = 1)
rtd1 <- build_tardis(rst)
cc1 <- cumulative_cost(rtd1, origin = pt)

# voronoi
rst2 <- link_islands(rst, klink = NULL)
rtd2 <- build_tardis(rst2)
cc2 <- cumulative_cost(rtd2, origin = pt)

# plot solutions together
par(mfrow = c(1, 2))
plot(cc1, zlim = c(0, 4e5), pal = sf.colors(20))
plot(rst$links, add = T)
plot(cc2, zlim = c(0, 4e5), pal = sf.colors(20))
plot(rst2$links, add = T)

## -----------------------------------------------------------------------------
# example

## -----------------------------------------------------------------------------
# weighting function for elevation-adjusted great circle distance (tardis default)
wfun <- function(origin, dest) {
  return(sqrt(origin$hdist^2 + abs(origin$vdist)^2))
}

# weighting function for the average of some additional value supplied
wfun <- function(origin, dest) {
  return((origin$value + dest$value) / 2)
}

# weighting function where base costs increase between layers according to a random multiplier
wfun <- function(origin, dest) {
  base_cost <- (origin$value + dest$value) / 2
  multiplier <- rgamma(n = length(base_cost), shape = origin$layer[1])
  return(base_cost * multiplier)
}

## -----------------------------------------------------------------------------
wfun <- function(origin, dest) {
  
  # calculate true geographic distance
  gdist <- sqrt(origin$hdist^2 + abs(origin$vdist)^2)
  
  # get the indices of the mask links
  mlink <- which(origin$type == 1)
  
  # multiply mask link costs by 100
  gdist[mlink] <- gdist[mlink] * 100
  
  # return costs
  return(gdist)
}

## -----------------------------------------------------------------------------
# apply the weighting
rtd <- weight_tardis(rtd, name = "mask100", wfun = wfun)

# view the scheme
head(rtd$edges[,"mask100"])

## -----------------------------------------------------------------------------
wfun <- function(origin, dest) {
  
  # get the slope incurred along an edge (rise over run)
  slope <- origin$vdist / origin$hdist
  
  # Tobler's hiking function
  speed <- 6 * exp(-3.5 * abs(slope + 0.05))
  
  # We want higher speeds to reflect a lower cost, so we return the reciprocal
  return(1 / speed)
}

## -----------------------------------------------------------------------------
wfun <- function(origin, dest) {
  
  # get the slope incurred along an edge (rise over run)
  slope <- origin$vdist / origin$hdist
  
  # Tobler's hiking function
  speed <- 6 * exp(-3.5 * abs(slope + 0.05))
  
  # multiply by normalised horizontal distance
  speed_norm <- speed * (origin$hdist / max(origin$hdist))
  
  # return as the reciprocal (resistance rather than conductance)
  return(1 / speed_norm)
}

## -----------------------------------------------------------------------------
# apply the weighting
rtd <- weight_tardis(rtd, name = "tobler", wfun = wfun)

# calculate cost surface
cc <- cost_surface(rtd, weights = "tobler")

# plot
plot(rst)
plot(cc)

## -----------------------------------------------------------------------------
# put our geoglist into a base R list
wvars <- list(elev = rst)

# a very simple cost function that penalises movement directly according to elevation
wfun <- function(origin, dest) {
  return(origin$elev)
}

# weight the tardis
rtd <- weight_tardis(rtd, name = "elev", vars = wvars, wfun = wfun)

