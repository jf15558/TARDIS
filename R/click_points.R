#' click_points
#'
#' Click locations on a geoglist layer and return the coordinates and cell IDs
#' of those locations. All locations are resolved to the nearest available
#' cell, including any which fall outside of unmasked areas of the layer.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param layer `numeric`. The layer in the geoglist to be plotted, along with
#' its links.
#' @param points `numeric`. The number of links you wish to generate.
#' Simply rerun the function if you need to add more lines.
#' @param print.only `logical`. Defaults to `TRUE` and the cell indices and values
#' are only printed to the console. If `FALSE`, then those points are returned
#' as a `SpatVector` object.
#' @param ... Additional arguments passed to `plot.geoglist()`
#' @return Either `NULL` (default) or a `SpatVector` of the clicked points.
#' @import terra sf h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' # load data
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # make geoglist
#' rasts <- rast_to_geoglist(gal, gal_m)
#'
#' # click to get two points
#' click_points(rasts, points = 2)
#' }

click_points <- function(geog, layer = 1, points = 1, print.only = T, ...) {
  #
  # geog = rasts
  # layer = 1
  # points = 2
  #
  if(!exists("geog")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!is.atomic(layer) | length(layer) != 1) {
    stop("layer should be a single integer")
  }
  if(!is.numeric(layer)) {
    stop("layer should be a single integer")
  }
  if(layer %% 1 != 0) {
    stop("layer should be a single integer")
  }
  if(layer > length(geog$layers)) {
    stop("The value of layer exceeds of the number of layers in geog")
  }

  if(!is.atomic(points) | length(points) != 1) {
    stop("points should be a single integer")
  }
  if(!is.numeric(points)) {
    stop("points should be a single integer")
  }
  if(points %% 1 != 0) {
    stop("points should be a single integer")
  }

  if(inherits(geog$layers[[1]], "SpatRaster")) {

    islands <- patches(geog$layers[[layer]], directions = 8, allowGaps = F)
    bounds <- as.points(mask(islands, classify(boundaries(islands), cbind(0, NA))))

  } else {

    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    bounds <- geog$layers[[layer]]
    bounds <- na.omit(bounds, field = names(bounds)[1])
    bar <- relate(bounds, bounds, "intersects", pairs = T)
    bounds <- centroids(bounds[which(table(bar[,1]) != 7)])
  }

  plot(geog, layer, ...)

  lnk <- list()
  for(i in 1:points) {

    crd <- click(n = 1)

    vl <- extract(geog$layers[[1]], crd)
    if(is.na(vl[2])) {
      crd <- crds(nearest(vect(crd, crs = crs(bounds)), bounds)[,5:6])
      vl <- extract(geog$layers[[1]], crd)
    }

    if(inherits(geog$layers, "SpatRaster")) {
      cl <- cellFromXY(geog$layers, crd)
    } else {
      cl <- suppressMessages(match(point_to_cell(crd, geog$gdat[7]), grid))
    }
    pt <- vect(crd)
    pt$cell <- cl
    pt$value <- vl[2]
    lnk[[i]] <- pt
  }
  lnk <- do.call(rbind, lnk)

  if(print.only) {
    print(as.data.frame(lnk))
  } else {
    return(lnk)
  }
}
