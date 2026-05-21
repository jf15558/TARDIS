#' click_lcp
#'
#' Interactively calculate and plot `n` least cost paths by clicking their start
#' and end points on the landscape displayed by the function. Clicked points
#' falling in masked regions are automatically resolved to the nearest available
#' cell.
#'
#' @param tardis `tardis`. An object of class 'tardis', produced by create_tardis
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' @param geog `geoglist`. A geoglist
#' @param time `integer`. The tardis time slice to plot and interact with.
#' Defaults to `NULL`, in which case the first slice is used.
#' @param n `integer`. the number of point pairs to run (only for click_lcp)
#' @param col `character`. The colour to use for plotting interactive features.
#' @import terra sf scales
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' hlink <- link_islands(hexes)
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # click two start and end points on the map
#' click_iso(tardis = rtd, geog = hexes, time = 2, n = 2)
#' }

click_lcp <- function(tardis, weights = "gdist", geog, time = NULL, n = 1, col = "gold") {

  # tardis = rtdw
  # geog = rasts
  # time = 2
  # col = "gold"
  # n = 1

  if(is.null(time)) {
    if(!is.null(tardis$tdat)) {
      time <- sum(tardis$tdat[1:2]) / 2
    }
    bin <- 1
  } else {
    if(time > tardis$tdat[1] | time < tardis$tdat[length(tardis$tdat)]) {
      stop("Time falls outside the range of tardis")
    }
    bin <- sum(time < tardis$tdat)
  }

  if(inherits(geog$layers[[1]], "SpatRaster")) {
    plot(geog$layers[[bin]])
  } else {
    plot(geog$layers[[bin]]$geometry, border = NA)
    plot(geog$layers[[bin]][,1], add = T)
  }
  if(!is.null(geog$links)) {
    plot(geog$links[which(geog$links$layer == bin),"geometry"], add = T)
  }

  org <- cbind(click(n = n), rep(time, n))
  dst <- cbind(click(n = n), rep(time, n))

  hpts <- point_check(tardis, rbind(org, dst))

  hlcp <- least_cost(tardis, weights = weights, origin = hpts[1:n,], dest = hpts[(n + 1):(n * 2),])

  plot(hpts$geometry, col = col, pch = 16, add = T)
  plot(st_wrap_dateline(hlcp$geometry, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")), add = T, col = col, lwd = 2)
}
