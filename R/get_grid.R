#' get_grid
#'
#' Return all H3 cell IDs at a specified resolution within a rectangular bounding
#' box.
#'
#' @param bbox `numeric`. A vector of four numbers denoting the xmin, xmax, ymin
#' and ymax coordinates of the bounding box within which cells are to be retrieved.
#' @param hex `integer`. An integer denoting the H3 resolution at which to retrieve cells.
#' This will be in the range of 1 to 15. Note that retrieving high resolution
#' grids (`hex > 4`) for large (i.e., approaching global) extents will take a
#' long time and may become prohibitive.
#' @return `character`. A vector of H3 cell IDs.
#' @import h3jsr
#' @export
#'
#' @details
#' The function can be viewed as the functional equivalent to generating a raster
#' grid of a specified resolution within that same bounding box. The main
#' difference is that the number of hexagonal cells in given bounding box cannot
#' be easily calculated, unlike for a raster grid. The function additionally
#' controls for misbehaviour of bounding boxes approaching global extents,
#' avoiding subsequent complications for H3 cell retrieval.
#'
#'
#' @examples
#' \dontrun{
#' library(rTARDIS)
#
#' get_grid(bbox = c(-10, 13.2, -0.5, 30), hex = 3)
#' }

get_grid <- function(bbox, hex) {

  # bbox = geog$gdat[1:4]
  # hex = geog$gdat[7]

  # control bounding box for extreme latitudes to avoid cell retrieval errors
  if(bbox[1] < -179.99) {bbox[1] <- -179.99}
  if(bbox[2] > 179.99) {bbox[2] <- 179.99}
  if(bbox[3] < -89.99) {bbox[3] <- -89.99}
  if(bbox[4] > 89.99) {bbox[4] <- 89.99}

  # get all cells if raster is essentially global (faster than the polygon method)
  if(bbox[1] <= -175 & bbox[2] >= 175 & bbox[3] <= -86 & bbox[4] >= 86) {
    return(unlist(get_children(get_res0(), as.vector(hex))))

    # otherwise retrieve all cells in bounding box (still fast for local rasters)
  } else {

    # construct bounding polygon
    len1 <- diff(bbox[1:2]) / 0.01
    len2 <- diff(bbox[3:4]) / 0.01
    bounds <- cbind(
      c(seq(bbox[1], bbox[2], length.out = len1), rep(bbox[2], len2), seq(bbox[2], bbox[1], length.out = len1),  rep(bbox[1], len2)),
      c(rep(bbox[4], len1), seq(bbox[4], bbox[3], length.out = len2), rep(bbox[3], len1), seq(bbox[3], bbox[4], length.out = len2))
    )
    bounds <- st_sfc(st_polygon(list(bounds)), crs = "EPSG:4326")
    return(polygon_to_cells(bounds, hex)[[1]])
  }
}
