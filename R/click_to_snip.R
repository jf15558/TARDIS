#' click_to_snip
#'
#' Interactively remove links from a geoglist by drawing lines across them.
#' Lines are drawn by clicking their start and end points on the geog layer
#' plotted in the plotting window.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param layer `numeric`. The layer in the geoglist to be plotted, along with
#' its links.
#' @param nsnips `numeric`. The number of snipping lines you wish to generate.
#' Simply rerun the function if you need to remove more lines, or add unwanted
#' remaining lines such that they do not intersect any links.
#' @return The input `geoglist` with the snipped links removed.
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
#' # make geoglist and add links
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#'
#' # click to either side of a link to remove it
#' click_to_snip(rasts)
#' }

click_to_snip <- function(geog, layer = 1, nsnips = 1) {

  # geog = rasts
  # layer = 1
  # nsnips = 1

  if(!exists("geog")) {
    stop("Supply geog as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(is.null(geog$links)) {
    message("geog contains no links and will be returned unmodified")

  } else {

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

    if(!is.atomic(nsnips) | length(nsnips) != 1) {
      stop("nsnips should be a single integer")
    }
    if(!is.numeric(nsnips)) {
      stop("nsnips should be a single integer")
    }
    if(nsnips %% 1 != 0) {
      stop("nsnips should be a single integer")
    }

    if(inherits(geog$layers[[1]], "SpatRaster")) {
      plot(geog$layers[[layer]])
    } else {
      plot(geog$layers[[layer]]$geometry, border = NA)
      plot(geog$layers[[layer]][,1], add = T)
    }
    plot(geog$links[which(geog$links$layer == layer),"geometry"], add = T)

    snips <- list()
    for(i in 1:nsnips) {
      ln <- st_sfc(st_linestring(click(n = 2)), crs = "+proj=lonlat")
      plot(ln, add = T, col = 2, lwd = 2)
      lnks <- which(geog$links$layer == layer)
      snips[[i]] <- lnks[unlist(st_intersects(ln, geog$links[lnks,]))]
    }
    snips <- unlist(snips)
    if(length(snips) == 0) {
      message("Snipping lines did not intersect with any links")
    } else {
      geog$links <- geog$links[-snips,]
    }
  }
  return(geog)
}
