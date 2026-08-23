## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
if(!require("chronosphere")) {
  install.packages("chronosphere",
                   repos = "http://cran.us.r-project.org")
}
if(!require("via")) {
  install.packages("via",
                   repos = "http://cran.us.r-project.org")
}

## ----setup, echo=F------------------------------------------------------------
library(terra)
library(rTARDIS)
library(via)
library(chronosphere)

## -----------------------------------------------------------------------------
library(chronosphere)
datasets()[which(datasets()$defaultClass == "RasterArray"),4]

## ----include=F----------------------------------------------------------------
dem <- fetch("paleomap", "dem")

## -----------------------------------------------------------------------------
tms <- dem@index
dem <- dem@stack
unloadNamespace("via")
dem

## ----include=FALSE------------------------------------------------------------
# load the dataset: slices 49 to 46 correspond to 385, 380, 375, and 370 million years ago
dvn <- dem[[78:75]]

# create the land-sea mask
dvn_l <- classify(dvn, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# convert to geoglist
hexes <- rast_to_geoglist(dvn, dvn_l, as.hex = T, hex = 3)

# plot
plot(hexes, 2)

## ----echo=F, fig.width=9.5, fig.height=4.5, dev="svg"-------------------------
# plot
plot(hexes, 2)

## ----include=FALSE------------------------------------------------------------
# link islands
hexes <- link_islands(hexes, klink = 1)

# plot
plot(hexes, 2)

## ----echo=F, fig.width=9.5, fig.height=4.5, dev="svg"-------------------------
# plot
plot(hexes, 2)

## -----------------------------------------------------------------------------
# construct our time vector
tms <- c(387.5, 382.5, 377.5, 372.5, 367.5)

# get rotations for the PALEOMAP model
rots <- get_rotations(hexes, tms, "PALEOMAP")

# interrogate the output
length(rots)
head(rots[[1]])

## -----------------------------------------------------------------------------
# build tardis
htd <- build_tardis(hexes, tms, rotations = rots)

