#' click functions
#'
#' interactive function under development
#' @param tardis `tardis`. An object of class 'tardis', produced by create_tardis
#' @param geog `geoglist`. A geoglist
#' @param time `integer`. The tardis time slice to plot and interact with
#' @param n `integer`. the number of point pairs to run (only for click_lcp)
#' @param col `character`. The colour to use for plotting interactive features.
#' @import terra sf scales

click_lcp <- function(tardis, geog, time = NULL, n = 1, col = "gold") {

  # tardis = rtdw
  # geog = rasts
  # time = 2
  # col = "gold"
  # n = 1

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

  org <- cbind(click(n = n), rep(time, n))
  dst <- cbind(click(n = n), rep(time, n))

  hpts <- point_check(tardis, rbind(org, dst))

  hlcp <- least_cost(tardis, origin = hpts[1:n,], dest = hpts[(n + 1):(n * 2),])

  plot(hpts$geometry, col = col, pch = 16, add = T)
  plot(st_wrap_dateline(hlcp$geometry, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")), add = T, col = col, lwd = 2)
}
