#' click_optim
#'
#' Interactively calculate and plot an isochrone around a clicked point on the#
#' landscape displayed by the function. Clicked points falling in masked regions are automatically resolved to the nearest available
#' cell.
#'
#' interactive function under development
#' @param tardis `tardis`. An object of class 'tardis', produced by create_tardis
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' @param geog `geoglist`. A geoglist
#' @param time `integer`. The tardis time slice to plot and interact with.
#' Defaults to `NULL`, in which case the first slice is used.
#' @param n `integer`. the number of points to link with an optimal path. This
#' should be greater than 2.
#' @param loop `logical`. Should the optimal route additionally be closed from
#' end point to start point to return a polygon? Defaults to `FALSE`.
#' @param col `character`. The colour to use for plotting interactive features.
#' @param ... Additional arguments passed to `plot.geoglist()`
#' @import terra sf
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' hlink <- link_islands(hexes)
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # click a point on the map
#' click_iso(tardis = rtd, geog = hexes, time = 2, cost = 1e5)
#' }

click_optim <- function(tardis, weights = "gdist", geog, time = NULL, n = 3, loop = FALSE, col = "gold", ...) {

  # tardis = rtd
  # geog = rasts
  # weights = "gdist"
  # time = 2
  # n = 3
  # col = "gold"
  # cost = 1e7

  if(n < 3) {
    stop("n should be greater than 2")
  }

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
  plot.geoglist(geog, bin, ...)

  org <- cbind(click(n = n), rep(time, n))

  hpts <- point_check(tardis, org)

  hlcp <- optim_route(tardis, points = hpts, loop = loop)

  plot(hpts$geometry, col = col, pch = 16, add = T)
  plot(st_wrap_dateline(hlcp$geometry, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")), add = T, col = col, lwd = 2)
}
