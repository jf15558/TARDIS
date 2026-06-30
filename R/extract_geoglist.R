#' extract_geoglist
#'
#' Extract values from the layers in a geoglist object, using an sf geometry
#' collection. This geometry is expected to be from another TARDIS function
#' (e.g, `least_cost()`, `isochrone()`, ect), but could also be a user-designed
#' geometry.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param geom `sf simple features collection`. A set of sf geometries which
#' will be used to extract values from geog.
#' @param layer `numeric`. If not NULL, then an integer specifying from which
#' layer in `geog` values are to be extracted. This argument is intended for
#' use with user-designed `geom` objects which do not already contain layer
#' assigments, unlike returns from other TARDIS functions.
#' @return `sf simple features collection`. A simple features collection of points
#' corresponding to the centroids of all cells in geog intersected by an input
#' geometry (denoted by `$feature`) in its specified layer (`$layer`), and the
#' value present at that point (`$value`).
#' @import sf terra h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' rasts <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#' rlink <- link_islands(rasts)
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2))
#' dst <- rbind(c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' rpts <- point_check(rtd, rbind(org, dst))
#' rlcp <- least_cost(rtd, origin = rpts[1,], dest = rpts[3,])
#' vals <- extract_geoglist(rasts, rlcp)
#' }

extract_geoglist <- function(geog, geom, layer = NULL) {

  # geog = rasts
  # geom = ln
  # layer = NULL

  if(!exists("geog")) {
    stop("Supply geog as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!exists("geom")) {
    stop("Supply geom as an sf object")
  }
  if(!is.null(layer)) {
    if(!is.atomic(layer) | length(layer) != 1) {
      stop("If not NULL, layer should be a single integer")
    }
    if(!is.numeric(layer)) {
      stop("If not NULL, layer should be a single integer")
    }
    if(!layer %% 1 != 0) {
      stop("If not NULL, layer should be a single integer")
    }
    geom$layer <- layer
  }
  if(is.null(geom$feature)) {
    geom$feature <- 1:nrow(geom$feature)
  }

  if(!inherits(geog$layers[[1]], "SpatRaster")) {
    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
  }

  vals <- lapply(1:max(geom$feature), function(x) {
    pth <- geom[which(geom$feature == x),]
    vals2 <- lapply(min(pth$layer):max(pth$layer), function(y) {
      prt <- pth[which(pth$layer == y),]
      lyr <- geog$layers[[y]]
      if(inherits(lyr, "SpatRaster")) {
        vl <- extract(lyr, vect(prt), cells = T)
        cls <- cbind.data.frame(x, y, vl[,2])
        colnames(cls) <- c("feature", "layer", "value")
        cls$geometry <- st_as_sf(vect(xyFromCell(lyr, vl$cell)))$geometry

      } else {
        vl <- unlist(st_intersects(prt, lyr))
        cls <- st_drop_geometry(lyr[vl,1])
        cls <- cbind.data.frame(x, y, st_drop_geometry(lyr[vl,1]))
        colnames(cls) <- c("feature", "layer", "value")
        cls$geometry <- cell_to_point(grid[as.numeric(rownames(lyr)[vl])], geog$gdat[7])$geometry
      }
      cls <- st_as_sf(cls)
      st_crs(cls) <- "+proj=lonlat"
      cls
    })
    do.call(rbind, vals2)
  })
  return(vect(do.call(rbind, vals)))
}
