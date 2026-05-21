#' click_to_snip
#'
#' Interactively add links to a geoglist by clicking their start and end points
#' in the plotting window. Clicked locations are automatically resolved to the
#' nearest edge cell in the geoglist.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param layer `numeric`. The layer in the geoglist to be plotted, along with
#' its links.
#' @param nlinks `numeric`. The number of links you wish to generate.
#' Simply rerun the function if you need to add more lines.
#' @return The input `geoglist` with added links.
#' @import terra sf h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' # load data
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # make geoglist
#' rasts <- rast_to_geoglist(gal, gal_m)
#'
#' # click the start and end points for two links
#' click_to_link(rasts, nlinks = 2)
#' }

click_to_link <- function(geog, layer = 1, nlinks = 1) {
  #
  # geog = rasts
  # layer = 1
  # nlinks = 2
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

  if(!is.atomic(nlinks) | length(nlinks) != 1) {
    stop("nlinks should be a single integer")
  }
  if(!is.numeric(nlinks)) {
    stop("nlinks should be a single integer")
  }
  if(nlinks %% 1 != 0) {
    stop("nlinks should be a single integer")
  }

  if(inherits(geog$layers[[1]], "SpatRaster")) {
    plot(geog$layers[[layer]])

    islands <- patches(geog$layers[[layer]], directions = 8, allowGaps = F)
    bounds <- as.points(mask(islands, classify(boundaries(islands), cbind(0, NA))))

  } else {
    plot(geog$layers[[layer]]$geometry, border = NA)
    plot(geog$layers[[layer]][,1], add = T)

    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    bounds <- geog$layers[[layer]]
    bar <- st_touches(bounds)
    bounds <- centroids(vect(bounds[sapply(bar, length) != 6,]))
  }
  if(!is.null(geog$links)) {
    plot(geog$links[which(geog$links$layer == layer),"geometry"], add = T)
  }

  lnk <- list()
  for(i in 1:nlinks) {

    crd <- click(n = 2)
    ln <- as.lines(nearest(vect(crd, crs = crs(bounds)), bounds)[,5:6])

    if(inherits(geog$layers[[1]], "SpatRaster")) {
      vals <- extract(geog$layers[[layer]], ln)
      isvalid <- na.omit(vals[,2])
    } else {
      isvalid <- unlist(st_intersects(st_as_sf(ln), geog$layers[[layer]]))
    }

    if(length(isvalid) == 2) {

      lnk[[i]] <- ln
      plot(ln, add = T, col = 2, lwd = 2)

    } else {
      message("Link intersected unmasked regions other than at its start and end points and was not added")
    }
  }

  if(all(sapply(lnk, is.null))) {
    message("No valid links were generated and geog will be returned unmodified")

  } else {
    lnk <- do.call(rbind, lnk[!sapply(lnk, is.null)])
    lnk <- st_as_sf(lnk, crs = "+proj=lonlat")

    if(inherits(geog$layers, "SpatRaster")) {

      cl <- cellFromXY(geog$layers[[1]], as.matrix(st_coordinates(lnk)[,1:2]))
      cls <- as.data.frame(matrix(cl, ncol = 2, byrow = T))

    } else {

      cl <- suppressMessages(point_to_cell(st_coordinates(lnk)[,1:2], res = geog$gdat[7]))
      cl <- match(cl, grid)
      cls <- as.data.frame(matrix(cl, ncol = 2, byrow = T))
    }

    colnames(cls) <- c("srt", "end")
    cls$layer <- layer
    cls$distance <- perim(vect(lnk))
    st_geometry(cls) <- lnk$geometry

    geog$links <- unique(rbind(geog$links, cls))
  }
  return(geog)
}
