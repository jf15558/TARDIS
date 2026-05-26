#' project_geoglist
#'
#' Transform a geoglist and its links to a CRS other than lon-lat. This may be
#' desirable for plotting purposes, but will render the resulting geoglist
#' incompatible with most TARDIS functions.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param crs `character`. A valid crs character string to which the geoglist
#' will be transformed
#' @return `geoglist`. The input geoglist transformed to the target crs.
#' @import sf terra
#' @export
#'
#' @examples
#' gal <- cretaceous()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts, klink = 1)
#'
#' regs <- project_geoglist(geog = rasts, crs = "+proj=eqearth")
#' plot.geolist(regs)

project_geoglist <- function(geog, crs) {

  # geog <- rasts
  # crs ="+proj=eqearth"

  if(!exists("geog")) {
    stop("Supply geog as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }

  if(!exists("crs")) {
    stop("Please supply a crs")
  }
  if(!is.atomic(crs) | length(crs) != 1) {
    stop("Please supply crs as a valid crs string")
  }
  if(!is.character(crs)) {
    stop("Please supply crs as a valid crs string")
  }

  if(inherits(geog$layers[[1]], "SpatRaster")) {
    geog$layers <- project(geog$layers, crs)
  } else {
    geog$layers <- lapply(geog$layers, st_transform, crs = crs("+proj=eqearth"))
  }
  if(!is.null(geog$links)) {
    geog$links <- st_transform(geog$links, crs)
  }

  return(geog)
}
