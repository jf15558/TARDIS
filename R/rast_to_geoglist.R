#' rast_to_geoglist
#'
#' Convert a set of rasters to a `geoglist` compatible with downstream TARDIS
#' functions. Typically, these rasters will record topography and/or bathymetry
#' measured in metres, but could record other geographic properties instead.
#' Some downstream functions will also assume that the set represents the same
#' geographic area ordered forwards in time (i.e., the first raster in the set
#' is the oldest).
#'
#' @param geog `SpatRaster`. A set of geographic rasters. These must be in
#' longitude-latitude projection. Missing values are not permitted and should
#' be replaced (e.g. with the average of the surrounding cells, or a dummy value).
#' @param mask `SpatRaster` or `NULL`. If not `NULL`, this will be used to
#' designate non-accessible areas in `geog`. It must be fully contiguous with
#' `geog` (i.e., share the same resolution, extent and number of layers) and
#' contain only `1` (unmasked, accessible) or `NA` (masked, non-accessible) values.
#' @param as.hex `logical`. Should `geog` be resampled to the hexagonal grid
#' system defined by Uber's H3 library? Defaults to `FALSE`.
#' @param hex `"auto"` or `integer`. The desired H3 resolution to be used for
#' resampling rasters. Defaults to `"auto"`, resulting in the function selecting
#' the H3 resolution closest to the resolution of `geog`. Otherwise an integer in
#' the range 1 - 15.
#' @param method `character`. The function to be used for resampling the raster grid.
#' This must be compatible with `exactextractr::exact_extract()`.
#' Defaults to `"mean"`.
#' @param verbose `logical`. Should function progress should be reported to the user?
#' @param ... Additional arguments passed internally to `exactextractr::exact_extract()`
#' for resampling of raster grids.
#' @return A `geoglist` with list elements `gdat` and `layers`. The former records
#' spatial properties of the input rasters used throughout downstream TARDIS functions.
#' The latter is a set of geographic layers, either as a `SpatRaster`, or an `sf data.frame` of
#' polygons if resampling was implemented.
#' @import terra exactextractr h3jsr
#' @export
#'
#' @details
#' Masking is a key feature of TARDIS. Often landscapes will contain
#' areas which we want to exclude from traversal. This is desirable for two reasons.
#' The first is that operations on these landscapes can be spatially constrained
#' in a realistic manner (e.g., restricting terrestrial organisms to the land surface,
#' or marine organisms to the oceans. The other is that masking can dramatically
#' reduce the number of cells and so memory required for landscape representation,
#' improving computational efficiency.
#'
#' Resampling to a hexagonal grid may be desirable when working with landscapes
#' approaching global extents in reduce the number of cells required to represent
#' polar latitudes, again helping to improve computational efficiency. The trade
#' off is that landscape features may be altered or lost depending on the grid
#' resolution used, although this an inevitable risk of any resampling procedure.
#' In addition, a large `SpatRaster` may be noticeably faster to plot than an
#' `sf data.frame` of similar extent and resolution, which can affect downstream
#' functions.
#'
#' Resampling is weighted by the fraction of each raster cell covered by
#' a given hexagonal cell. This is implemented using `exactextractr::exact_extract()`
#' as this is currently much faster compared to the terra's native weighted resampling
#' function. Weighting is not meaningful for all functions, e.g., `"max"`.
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # create a geoglist with hexagonal resampling
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#'
#' # create a geoglist using rasters directly
#' rasts <- rast_to_geoglist(gal, gal_m)
#' }

rast_to_geoglist <- function(geog, mask = NULL, as.hex = FALSE, hex = "auto", method = "mean", verbose = T, ...) {

  # geog  = gal
  # mask = gal_m
  # as.hex = T
  # hex = 6
  # method = "mean"
  # verbose = T

  # check geography
  if(!exists("geog")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!inherits(geog, "SpatRaster")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!is.lonlat(geog)) {
    stop("geog should be in geographic (long-lat) projection")
  }
  crs(geog) <- "EPSG:4326"
  if(any(is.na(geog[]))) {
    stop("NA values present in geog")
  }

  if(!is.null(mask)) {

    if (!inherits(mask, "SpatRaster")) {
      stop("If not NULL, supply mask as a SpatRaster")
    }
    if (any(!unique(mask[]) %in% c(1, NA))) {
      stop("Mask layers can only contain 1 or NA values")
    }
    if (!is.lonlat(mask)) {
      stop("mask should be in geographic (lon-lat) projection")
    }
    crs(mask) <- "EPSG:4326"
    if (!all(dim(geog) == dim(mask)) | ext(geog) != ext(mask)) {
      stop("geog and mask must all have the same extent, resolution, and number of layers")
    }
    geog <- mask(geog, mask)
  }
  if(!is.atomic(as.hex) | length(as.hex) != 1) {
    stop("as.hex should be a logical indicating if to resample the raster grid to a hexagonal grid")
  }

  if(as.hex) {

    res <- max(cellSize(geog)[])
    res <- which.min(abs(res - h3jsr::h3_info_table$avg_area_sqm)) - 1
    if(!is.atomic(hex) | length(hex) != 1) {
      stop("hex should either be set to 'auto', or an integer specifying the desired H3 resolution for resampling [1-15]")
    }
    if(hex != "auto") {
      if(!hex %in% 0:15) {
        stop("If not 'auto', then hex should be a number specifying the desired H3 resolution for resampling [1-15]")
      }
      if(hex > res) {
        warning("The chosen H3 resolution is much smaller than the raster grid resolution. Resampling could take a while")
      }
    } else {
      cat(paste0("Autoselecting H3 resolution ", res, " "))
      hex <- res
    }

    if(!is.null(mask)) {
      grid <- lapply(mask, function(x) {xyFromCell(x, which(x[] == 1))})
      pol <- lapply(1:nlyr(mask), function(x) {st_as_sf(as.polygons(mask[[x]]))})
    } else {
      pol <- grid <- lapply(1:length(geog), function(x) {xyFromCell(geog[[1]], 1:ncell(geog))})
    }

    clist <- get_grid(as.vector(ext(geog)), hex)
    hex_list <- list()
    for(i in 1:nlyr(geog)) {

      if (verbose) {
        cat(paste0("Resampling layer [", i, "/", nlyr(geog), "]\r"))
        if (i == nlyr(geog)) {cat("\n")}
      }

      # both polygon and point to ensure that no cells fail to locate to a hex
      cls <- unique(na.omit(c(unlist(suppressMessages(polygon_to_cells(pol[[i]], hex))),
                              unlist(suppressMessages(point_to_cell(grid[[i]], hex))))))

      # very rarely, some cells are recovered which do not lie in bounds - these are dropped
      cls <- intersect(cls, clist)
      clsp <- cell_to_polygon(cls)
      clsp <- st_make_valid(clsp)
      clsp <- st_wrap_dateline(clsp, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
      id <- match(cls, clist)
      vrs <- exact_extract(geog[[i]], clsp, fun = method, weights = "area", ...)
      dat <- data.frame(vrs)
      colnames(dat)[1] <- names(geog[[i]])
      dat$id <- id
      st_geometry(dat) <- clsp
      hex_list[[i]] <- vect(dat[st_is_valid(dat),])
    }
    out <- list(gdat = c(as.vector(ext(geog)), ncell = length(clist), ncol = NA, hex = hex), layers = svc(hex_list))

  } else {
    out <- list(gdat = c(as.vector(ext(geog)), ncell = ncell(geog), ncol = ncol(geog), hex = NA), layers = geog)
  }
  class(out) <- "geoglist"
  return(out)
}
