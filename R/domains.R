#' domains
#'
#' Calculate Voronoi polygons around islands in a geoglist. These 'domains' are
#' the basis for how `link_islands` determines connections across masked areas
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @return A `SpatVectorCollection` with each element containing the Voronoi
#' polygons for the corresponding layer in the input geoglist.
#' @import terra
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- galapagos()
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # create a geoglist with hexagonal resampling and mask the sea
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#'
#' dms <- domains(hexes)
#' plot(hexes, 1)
#' plot(dms[[1]], add = T)
#' }

domains <- function(geog) {

  #geog = rasts

  if(!exists("geog")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(inherits(geog$layers, "SpatRaster")) {

    bounds <- lapply(geog$layers, function(x) {
      islands <- patches(x, directions = 8, allowGaps = F)
      pt <- as.points(mask(islands, classify(boundaries(islands), cbind(0, NA))))
      aggregate(pt, by = "patches")
    })

  } else {

    bounds <- lapply(geog$layers, function(z) {
      z <- na.omit(z, field = names(z)[1])
      bar <- relate(z, z, "intersects", pairs = T)
      z$patches <- components(graph_from_edgelist(bar))$membership
      z2 <- centroids(z[which(table(bar[,1]) != 7)])
      aggregate(z2, by = "patches")
    })
  }

  vors <- lapply(bounds, function(x) {
    vv <- voronoi(x, bnd = x)
    aggregate(vv, by = "patches")
  })
  return(svc(vors))
}
