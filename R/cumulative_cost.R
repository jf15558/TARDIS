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
#' @return A `geoglist` recording the costs required to reach, from the origin
#' point, each cell in the same time slice as that point.
#' @import terra sf cppRouting h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' # load data
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # build a geoglist and add links
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' hexes <- link_islands(hexes)
#'
#' # build a tardis graph
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # resolve a destination point to that graph
#' org <- rbind(c(-89.78873, -1.420627, 2))'
#' hpts <- point_check(htd, org)
#'
#' # get the cumulative cost around that point
#' cumulative_cost(htd, origin = hpts)
#' }

cumulative_cost <- function(tardis, weights = "gdist", origin, verbose = TRUE) {

   #tardis = rtd
   #weights = "gdist"
   #origin = pt[1,]
   #verbose = TRUE

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

  edges <- tardis$edges[which(tardis$edges[,3] != 2),]
  dests <- edges[which((edges[,2] %/% tardis$gdat[[5]]) + 1 == origin$bin),2]


  if (verbose) {
    cat("Running cumulative costs\r")
  }
  paths <- get_distance_matrix(tardis$tgraph, from = origin$cell, to = dests)[1,]

  if(!is.na(tardis$gdat[7])) {
    tmp <- vect(st_wrap_dateline(cell_to_polygon(grid[as.numeric(names(paths)) %% tardis$gdat[5]]),
                            options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")))
    tmp$distance <- paths
    names(tmp)[1] <- weights
    tmp <- svc(tmp)

  } else {
    tmp <- samprast
    tmp[] <- NA
    names(tmp) <- weights
    tmp[][as.numeric(names(paths)) %% tardis$gdat[5], 1] <- paths
  }
  out <- list(gdat = tardis$gdat, layers = tmp)
  class(out) <- "geoglist"
  return(tmp)
}
