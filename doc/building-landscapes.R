## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=F------------------------------------------------------------
library(terra)
library(rTARDIS)

## -----------------------------------------------------------------------------
# convert a toy dataset to a geoglist and view its structure
rasts <- rast_to_geoglist(galapagos())
rasts

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# plot the first layer of the geoglist above
plot(rasts, 1)

## ----fig.width=9.5, fig.height=5, dev="svg"-----------------------------------
# load a dataset of the Galapagos archipelago and select the present day layer
gal <- galapagos()[[5]]

# create a land-sea mask from the archipelago raster set
gal_l <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# plot original data and the land-sea mask
par(mfrow = c(1, 2))
plot(gal)
plot(gal_l)

## ----fig.width=9.5, fig.height=5, dev="svg"-----------------------------------
# create a shallow sea mask at -200 to 0 metres (the photic zone) from the archipelago raster set
gal_m <- classify(gal, matrix(c(-Inf, -200, NA, -200, 0, 1, 0, Inf, NA), ncol = 3, byrow = T), right = F)

# plot original data and the land-sea mask
par(mfrow = c(1, 2))
plot(gal)
plot(gal_m)

## ----fig.width=9.5, fig.height=5, dev="svg"-----------------------------------
# create a land mask
land <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# get all cells around the boundaries of land cell areas and find their indices
coast <- boundaries(land)
coast <- which(coast[] == 1)

# add those cell indices to our previous shallow sea mask
gal_m[coast] <- 1

# plot original data and the land-sea mask
par(mfrow = c(1, 2))
plot(gal)
plot(gal_m)

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
rst <- rast_to_geoglist(gal, mask = gal_l)
plot(rst)

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# detect regions (bounds.only = F means that all island cells are returned, rather than just their margins)
rg <- regions(rst, bounds.only = F)

# get the voronoi domains for those regions
dm <- domains(rst)

# plot together, cropping to the voronoi domain extent
plot(rst, xlim = ext(dm)[1:2], ylim = ext(dm)[3:4])
plot(dm, add = T)

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# minimum spanning tree linkage, then plot
rst <- link_islands(rst, klink = 1)
plot(rst, xlim = ext(dm)[1:2], ylim = ext(dm)[3:4])

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# link then plot
rst <- link_islands(rst, klink = 2)
plot(rst, xlim = ext(dm)[1:2], ylim = ext(dm)[3:4])

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# link then plot
rst <- link_islands(rst, klink = NULL)
plot(rst, xlim = ext(dm)[1:2], ylim = ext(dm)[3:4], lcol = "gold", lwd = 2)
plot(dm, add = T, border = "grey")

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# link then plot
rst <- link_delaunay(rst)
plot(rst)

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# link then plot
regs <- regions(rst, use.links = T, bounds.only = F)
plot(regs)

## ----include=FALSE------------------------------------------------------------
# load the dataset and select the first time slice
cret <- cretaceous()[[1]]

# create the land-sea mask
cret_l <- classify(cret, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# convert to geoglist with h = 2, and h = 3
h2 <- rast_to_geoglist(cret, cret_l, as.hex = T, hex = 2)
h3 <- rast_to_geoglist(cret, cret_l, as.hex = T, hex = 3)

# plot to compare (zlim set to)
plot(h2)
plot(h3)

## ----fig.width=9.5, fig.height=4.5, dev="svg"---------------------------------
plot(h2)
plot(h3)

## -----------------------------------------------------------------------------
# voronoi linkage (klink is NULL by default)
h3 <- link_islands(h3)

## ----fig.width=9.5, fig.height=10, dev="svg"----------------------------------
# project to Equal Earth (you can look up proj strings for other map types online)
heq <- project_geoglist(h3, crs = "+proj=eqearth")

# plot
plot(heq, bg = "aliceblue")

## -----------------------------------------------------------------------------
rtd <- build_tardis(rst)

## -----------------------------------------------------------------------------
head(rtd$edges)

## -----------------------------------------------------------------------------
# get a random sample of cell IDs
smp <- sample(rtd$edges[,1], 10)

# get their temporal position (i.e., layer 1 of 1)
smp %/% rtd$gdat["ncell"] + 1

# get their spatial position (raster cell number). As we expect, these values are unchanged as they come from the first (and only) layer in the geoglist
smp %% rtd$gdat["ncell"]

## -----------------------------------------------------------------------------
head(rtd$tgraph$dict)

## ----echo=FALSE, fig.width=9.5, fig.height=4.5, dev="svg"---------------------
# load data
gal <- galapagos()

# create a land-sea mask from the archipelago raster set
gal_l <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

# plot original data and the land-sea mask
par(mfrow = c(1, 2))
plot(gal[[1]])
plot(gal[[3]])

## -----------------------------------------------------------------------------
# create our time vector
tms <- c(seq(2.25, 0, -0.5), 0)

# build a tardis graph from the masked landscape
rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))

# view the temporal information
rtd$tdat
rtd$tlink

# view the link types
unique(rtd$edges[,"type"])

## -----------------------------------------------------------------------------
# consistent seed for the vignette output
set.seed(1)

# sample cells randomly
smp <- sample(rtd$edges[,1], 100)

# get their temporal position (multiple layers represented)
tmp <- smp %/% rtd$gdat["ncell"] + 1

# get their spatial position (raster cell number).
spt <- smp %% rtd$gdat["ncell"]

# convert spatial position to coordinates
crds <- xyFromCell(gal, spt)

# print the first 10 as examples
cbind(smp, tmp, spt, crds)[1:10,]

