## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=F------------------------------------------------------------
library(sf)
library(h3jsr)
library(terra)
library(rTARDIS)

## ----echo=FALSE, fig.width=9.5, fig.height=4.5, dev="svg"---------------------
foo <- as.polygons(rast(ncol = 46, nrow = 23))
foo$area <- expanse(foo) / max(expanse(foo))
foo2 <- st_sf(cell_to_polygon(get_children(get_res0(), res = 1)))
foo2 <- suppressWarnings(st_wrap_dateline(foo2, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")))
foo2 <- vect(foo2)
foo2$area <- expanse(foo2) / max(expanse(foo2))
plot(foo, values = foo$area, range = c(0,1), breaks = 20, legend = F, cex = 0.7)
legend_cont("topleft", legend = c(0, 1))
plot(foo2, values = foo2$area, range = c(0,1), breaks = 20, legend = F, xlim = c(-180, 180), ylim = c(-90, 90))
legend_cont("topleft", legend = c(0, 1))

## ----echo=F, fig.width=9.5, dev="svg", include=F------------------------------
par(mfrow = c(1, 2))
plot(st_as_sf(as.polygons(rast(ncol = 3, nrow = 3))), col = c(2, 4, 2, 4, 5, 4, 2, 4, 2))
plot(cell_to_polygon(get_ring(suppressMessages(point_to_cell(c(10, -40), 2)))), col = 2)
plot(cell_to_polygon(suppressMessages(point_to_cell(c(10, -40), 2))), add = T, col =5)

## ----echo=F, fig.width=9.5, fig.height=4.5, dev="svg"-------------------------
foo <- st_sf(cell_to_polygon(polygon_to_cells(
  st_sfc(
    st_polygon(list(cbind(c(-10, -10, -9, -9, -10),
                          c(-10, -9, -9, -10, -10)))), 
    crs = "EPSG:4326"), res = 6)))
foo$id <- 1:nrow(foo)
foo2 <- rast(xmin = -10, xmax = -9, ymin = -10, ymax = -9, ncol = 20, nrow = 20)
foo2$id <- 1:ncell(foo2)
par(mfrow = c(1, 2))
plot(as.polygons(foo2), value = foo2$id, breaks = 20, legend = F)
legend_cont("topleft", legend = foo2)
plot(vect(foo), value = foo$id, breaks = 20, legend = F)

## ----echo=F, fig.width=9.5, fig.height=5, dev="svg"---------------------------
# load and plot galapagos dataset
gal <- galapagos()
gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
par(mfrow = c(1, 2))
plot(gal[[1]])
plot(gal_m[[1]], legend = F)

## ----include=F----------------------------------------------------------------
# load and plot galapagos dataset
rasts <- rast_to_geoglist(gal[[1]], gal_m[[1]])
rasts <- link_islands(rasts, klink = 1)

## ----echo=F, fig.width=9.5, fig.height=10, dev="svg"--------------------------
plot(rasts, legend = F)

## ----echo=F, fig.width=9.5, fig.height=5, dev="svg"---------------------------
par(mfrow = c(1, 2))
plot(gal_m[[1]], legend = F)
plot(gal_m[[2]], legend = F)

## ----echo=F, fig.width=9.5, fig.height=10, dev="svg"--------------------------
rtd <- build_tardis(rasts)
pts <- point_check(rtd, rbind(c(-88.97, -1.05), c(-90.77, 0.57)))
lcp <- least_cost(rtd, origin = pts[1,], dest = pts[2,])
plot(rasts)
plot(pts, col = "gold", add = T)
plot(lcp, col = "gold", add = T, lwd = 2)

