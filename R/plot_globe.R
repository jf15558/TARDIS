#' plot_globe
#'
#' Plot a geoglist on a 3D interactive globe. Geoglists of any spatial scale can
#' be plotted, but is primarily intended for those which cover the total surface
#' of the Earth or at least a decent portion of one of its hemispheres.
#'
#' @param x `geoglist`. The output of `rast_to_geoglist()`.
#' @param y `numeric`. The layer in the geoglist to be plotted, along with
#' its links. Defaults to 1 (the first layer).
#' @param range `vector`. If not `NULL`, the minimum and maximum values to be
#' plotted.
#' @param pal `vector`. A vector of colours to be used for plotting layer values,
#' such as those returned by an R colour palette. The length of this colour palette
#' will also correspond to the number of breaks used when plotting a geoglist
#' with continuous values.
#' @param links `logical`. Should mask links be plotted, if available?
#' @param lcol `character` or `integer`. The colour to be used for plotting links.
#' @param lwd `integer`. The line width to be used for plotting links.
#' @param lty `integer`. The line type to be used for plotting links.
#' @param bg `character` or `integer`. The colour to use for the globe background.
#' @param graticule `logical`. Should a lon-lat graticular be added to the sphere?
#' Defaults to `TRUE`.
#' @param grat.col `character` or `integer`. The colour to be used for the graticule lines.
#' @param add `logical`. Should the geoglist data be added to an existing rgl window?
#' #' @param ... Other arguments passed to `rgl` primitive plotting functions. Note
#' that some of these are set internally and so using this argument may cause errors.
#' Defaults to `FALSE`, which will intiate plotting on a blank sphere in a new
#' rgl window.
#' @return None.
#' @import sf terra rgl h3jsr
#' @importFrom  h3r cellToBoundary
#' @importFrom  h3r cellToLatLng
#' @export
#'
#' @details
#' Longitude-latitude coordinates are converted to spherical coordinates internally
#' as follows, using an authalic Earth radius of 6371.007 km
#'
#' `radius <- 6371.007`
#' `x <- radius * cos(lat) * cos(long)`
#' `y <- radius * cos(lat) * sin(long)`
#' `z <- radius * sin(lat)`
#'
#' The background sphere is produced using a slightly smaller radius to prevent
#' occulusion and artefacts in features plotted onto its surface

#' @examples
#' \dontrun{
#' gal <- cretaceous()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' rasts <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 3)
#' rasts <- link_islands(rasts, klink = NULL)
#' plot_globe(rasts, 1)
#' }

plot_globe <- function(x, y = 1, range = NULL, pal = sf.colors(10), links = T, lcol = "steelblue",
                       lwd = 1, lty = 1, bg = "aliceblue", graticule = T, grat.col = "grey95", add = F, ...) {

  # x = rst
  # y = 1
  # range = NULL
  # pal = sf.colors(10)
  # links = T
  # lcol = "steelblue"
  # lwd = 1
  # lty = 1
  # hex.border = NA
  # bg = "aliceblue"

  if(!inherits(x, "geoglist")) {
    stop("Please supply x as a geoglist object")
  }
  if(!is.lonlat(x$layers[[1]])) {
    stop("geoglist must be in lon-lat projection")
  }

  map2color <- function(x, pal, limits = range) {
    if(is.null(limits)) limits = range(x)
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
  }

  if(inherits(x$layers, "SpatRaster")) {

    grd <- as.polygons(x$layers[[y]], aggregate = F)
    grd <- crds(grd)
    grd <- grd[-(seq(5, nrow(grd), 5)),]

    #Defaults to the R2 radius of Earth (6371.007km).
    xVar <- cos(grd[,2] / 180 * pi) * cos(grd[,1] / 180 * pi)
    yVar <- cos(grd[,2] / 180 * pi) * sin(grd[,1] / 180 * pi)
    zVar <- sin(grd[,2] / 180 * pi)
    quads <- cbind(xVar, yVar, zVar) * 6371.007

    cls <- x$layers[[y]][][,1]
    cls <- cls[is.finite(cls)]
    cls <- map2color(cls, pal)

  } else {

    grd <- get_grid(x$gdat[1:4], x$gdat[7])
    grd <- grd[which(!is.na(x$layers[[y]][,1][,1]))]
    grd1 <- get_disk(grd)
    tris <- lapply(grd, function(x) {
      crds <- cellToBoundary(x)[[1]]
      crds <- rbind(cellToLatLng(x), crds)
      crds <- crds[cbind(sapply(2:(nrow(crds) - 1), function(y) {c(1, y, y + 1)}), c(1, nrow(crds), 2)),2:1]

      xVar <- cos(crds[,2] / 180 * pi) * cos(crds[,1] / 180 * pi)
      yVar <- cos(crds[,2] / 180 * pi) * sin(crds[,1] / 180 * pi)
      zVar <- sin(crds[,2] / 180 * pi)
      cbind(xVar, yVar, zVar) * 6371.007
    })

    cls <- x$layers[[y]][[1]][,1]
    cls <- cls[is.finite(cls)]
    cls <- map2color(cls, pal)
    cls <- rep(cls, times = sapply(tris, nrow))
    tris <- do.call(rbind, tris)
  }


  # base sphere
  if(!add) {
    res = 100
    lat <- matrix(seq(90, -90, len = res) * pi / 180, res, res, byrow = TRUE)
    long <- matrix(seq(-180, 180, len = res) * pi / 180, res, res)

    # slightly smaller than geographic radius to prevent occlusion and artefacts during plotting
    radius <- 6370
    xVar <- radius * cos(lat) * cos(long)
    yVar <- radius * cos(lat) * sin(long)
    zVar <- radius * sin(lat)
    persp3d(xVar, yVar, zVar, col = "lightblue", axes = F, box = F,
            shininess = 128, specular = "grey10", ambient = "steelblue",
            xlab = "", ylab = "", zlab = "", ...)
  }

  # add geoglist layer
  if(inherits(x$layers, "SpatRaster")) {
    quads3d(quads[,1], quads[,2], quads[,3], col = rep(cls, each = 4),
            shininess = 1, specular = "grey20", ...)

  } else {
    triangles3d(tris[,1], tris[,2], tris[,3], col = cls,
                shininess = 1, specular = "grey20", ...)
  }

  if(links) {
    if(!is.null(x$links)) {
      lns <- x$links[which(x$links$layer == y)]
      if(length(lns) != 0) {
        lns <- crds(lns)
        lns <- lapply(seq(1, nrow(lns), 2), function(x) {
          crd <- rbind(lns[x,], gcIntermediate(lns[x,], lns[x + 1,]), lns[x + 1,])

          xVar <- cos(crd[,2] / 180 * pi) * cos(crd[,1] / 180 * pi)
          yVar <- cos(crd[,2] / 180 * pi) * sin(crd[,1] / 180 * pi)
          zVar <- sin(crd[,2] / 180 * pi)
          crd <- cbind(xVar, yVar, zVar) * 6371.007

          crd[c(1, rep(2:(nrow(crd) - 1), each = 2), nrow(crd)),]
        })
        lapply(lns, segments3d, col = "steelblue", lwd = 2, ...)
      }
    }
  }

  if(graticule) {

    crd <- st_coordinates(st_graticule())
    xVar <- cos(crd[, 2]/180 * pi) * cos(crd[, 1]/180 * pi)
    yVar <- cos(crd[, 2]/180 * pi) * sin(crd[, 1]/180 * pi)
    zVar <- sin(crd[, 2]/180 * pi)
    crds <- cbind(cbind(xVar, yVar, zVar) * 6375, crd[,3])
    crds <- lapply(1:max(crds[,4]), function(x) {
      y <- crds[which(crds[,4] == x),1:3]
      y[c(1, rep(2:(nrow(y) - 1), each = 2), nrow(y)),]
    })

    lapply(crds, rgl::segments3d, col = grat.col)
  }
}
