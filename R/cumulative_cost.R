#' cumulative_cost
#'
#' From a given starting point, calculate the cost of reaching every other location
#' in the same time slice in a tardis graph and return the cumulative cost surface
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param origin `sf data.frame` A simple features collection produced by `point_check()`,
#' denoting a single point from which to calculate cumulative costs
#' @param verbose `logical` Should function progress be reported to the user?
#' @return A `geoglist` with a single layer recording the costs required to reach
#' each cell in the same time slice as the origin point.
#' them from the origin point.
#' in each landscape layer through time
#' @import terra sf cppRouting h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' rasts <- rast_to_geoglist(gal, gal_m)

#' hlink <- link_islands(hexes)
#' rlink <- link_islands(rasts)

#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))

#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2))
#' dst <- rbind(c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- point_check(htd, rbind(org, dst))
#' rpts <- point_check(rtd, rbind(org, dst))
#'
#' hlcp <- lcp(htd, origin = hpts[1:2,], dest = hpts[3:4,])
#' rlcp <- lcp(rtd, origin = rpts[1:2,], dest = rpts[3:4,])
#' }

cumulative_cost <- function(tardis, weights = "gdist", origin, verbose = TRUE) {
  #
  # tardis = rtd
  # weights = "tdist"
  # origin = org
  # verbose = TRUE

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.na(tardis$gdat[7])) {
    grid <- get_grid(tardis$gdat[1:4], tardis$gdat[7])

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
  }

  if (!class(origin)[1] == c("sf")) {
    stop("Supply origin as the output of stp")
  }
  if(!is.na(tardis$gdat[7])) {

    pcell <- point_to_cell(origin, tardis$gdat[7])
    origin$cell <- match(pcell, grid) + (tardis$gdat[5] * (origin$bin - 1))

  } else {

    origin$cell <- cellFromXY(samprast, st_coordinates(origin)) +
      (tardis$gdat[5] * (origin$bin - 1))
  }

  if (!all(origin$cell %in% tardis$edges[, 1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }
  if (verbose) {
    cat("Initialising graph\n")
  }
  tardis <- instantiate_tardis(tardis = tardis, weights = weights)

  edges <- tardis$edges[which(tardis$edges[,3] == 0),]
  dests <- edges[which(edges[,2] %/% tardis$gdat[[5]] + 1 == origin$bin),2]


  if (verbose) {
    cat("Running cumulative costs\r")
  }
  paths <- get_distance_matrix(tardis$tgraph, from = origin$cell, to = dests)[1,]

  if(!is.na(tardis$gdat[7])) {
    cls <- st_wrap_dateline(cell_to_polygon(grid[as.numeric(names(paths)) %% tardis$gdat[5]]),
                            options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    tmp <- data.frame(paths)
    colnames(tmp) <- weights
    st_geometry(tmp) <- cls

  } else {
    tmp <- samprast
    tmp[] <- NA
    names(tmp) <- weights
    tmp[][as.numeric(names(paths)) %% tardis$gdat[5], 1] <- paths
  }

  out <- list(gdat = tardis$gdat, layers = list(tmp))
  class(out) <- "geoglist"
  return(out)
}
