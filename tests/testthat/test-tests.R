## SETUP ##

# load required libraries
library(terra)
library(TARDIS)
library(testthat)

# set up objects
gal <- galapagos()
gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
hexes <- link_islands(hexes, klink = 1)


## FUNCTIONS ##

# build_tardis
htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))

# cost_surface

# cretaceous

# cumulative_cost

# extract_geoglist

# galapagos

# get_grid

# get_rotations

# instantiate_tardis

# isochrone

# least_cost

# link_delauney

# link_islands

# min_span

# optim_route

# point_check

# project_geoglist

# rand_walk

# rast_to_geoglist

# regions

# slice_tardis

# tardis_to_matrix

# weight_tardis
