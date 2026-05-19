#' click functions
#'
#' interactive function under development
#' @param tardis `tardis`. An object of class 'tardis', produced by create_tardis
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' @param geog `geoglist`. A geoglist
#' @param time `integer`. The tardis time slice to plot and interact with
#' @param col `character`. The colour to use for plotting interactive features.
#' @param cost `numeric`. the maximum cost for calculating a reach isochrone (only for click_iso)
#' @import terra sf scales

click_iso <- function(tardis, weights = "gdist", geog, time = NULL, cost = 1e6, col = "gold") {

  # tardis = rtd
  # geog = rasts
  # weights = "gdist"
  # time = 115
  # n = 1
  # col = "gold"
  # cost = 1e7

  if(is.null(time)) {
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
    plot(geog$links[which(geog$links$bin == bin),"geometry"], add = T)
  }

  org <- cbind(click(n = 1), rep(time, 1))

  hpts <- point_check(tardis, org)

  hlcp <- isochrone(tardis, origin = hpts, cost = cost)

  plot(hpts$geometry, col = col, pch = 16, add = T)
  plot(hlcp$geometry, add = T, border = col, col = scales::alpha(col, 0.2), lwd = 2)
}
