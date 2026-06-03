#' isochrone
#'
#' Calculate a cost isochrone around a point, that is all cells in the same
#' time layer which fall within a certain cost of access of that point.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param origin `SpatVector`. The output of `point_check()`, denoting the points
#' around which to calculate isochrones.
#' @param cost The maximum cumulative cost of below which cells will be
#' included within the isochrone in a point's time-specific origin landscape
#' layer, either a single number or a vector with as many elements as points in
#' origin to enable different costs for each point. By default, this value is
#' in metres, corresponding to the TARDIS default weighting scheme.
#' @param verbose `logical` Should function progress be reported to the user?
#' @return A `SpatVector` containing isochrone polygons, recording to which input
#' point in `origin` they correspond to (`$feature`) and their time layer (`$layer`)
#' @import terra sf h3jsr
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
#' # create a geoglist and link islands
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#'
#' org <- rbind(c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- point_check(htd, org)
#' foo <- isochrone(htd, origin = hpts, cost = 100000)
#' }

isochrone <- function(tardis, weights = "gdist", origin, cost = 1e5, verbose = TRUE) {

  # tardis = rtd
  # weights = "gdist"
  # origin = point_check(rtd, rbind(c(-160.559, 82.33975, 115)))
  # cost = 1e7
  # restrict = T
  # verbose = T

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

  if (!class(origin)[1] == c("SpatVector")) {
    stop("Supply origin as the output of point_check")
  }

  if(!is.na(tardis$gdat[7])) {

    pcell <- point_to_cell(crds(origin), tardis$gdat[7])
    origin$cell <- match(pcell, grid) + (tardis$gdat[5] * (origin$layer - 1))

  } else {

    origin$cell <- cellFromXY(samprast, crds(origin)) +
      (tardis$gdat[5] * (origin$layer - 1))

  }

  if (!all(origin$cell %in% tardis$edges[, 1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }

  if(!is.atomic(cost) | !is.numeric(cost)) {
    stop("If not NULL, cost should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }
  if(length(cost) == 1) {
    cost <- rep(cost, nrow(origin))
  }
  if(length(cost) != length(cost)) {
    stop("If not NULL, cost should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }
  if(!all(is.finite(cost) | cost <= 0)) {
    stop("If not NULL, cost should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }

  if (verbose) {
    cat("Initialising graph\n")
  }
  tardis <- instantiate_tardis(tardis, weights)

  ob_list <- t_list <- list()
  for(i in 1:nrow(origin)) {

    if(verbose) {cat(paste0("Running isochrones [", i, "/", nrow(origin), "]\r"))
      if(i == nrow(origin)) {cat("\n")}
    }
    res <- as.numeric(get_isochrone(tardis$tgraph, from = origin$cell[i], lim = cost[i])[[1]])
    rwt <- res %/% tardis$gdat[5] + 1
    rwp <- res %% tardis$gdat[5]
    rwp <- rwp[which(rwt == origin$cell[i] %/% tardis$gdat[5] + 1)]
    rwt <- rwt[which(rwt == origin$cell[i] %/% tardis$gdat[5] + 1)]

    out <- tapply(rwp, rwt, function(x) {
      if(!is.na(tardis$gdat[7])) {

        baz <-  st_wrap_dateline(cell_to_polygon(grid[x]), options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
        iso_xy <- st_sfc(st_union(baz))

      } else {
        tmprast <- samprast
        tmprast[x] <- 1
        iso_xy <- st_as_sf(as.polygons(tmprast))$geometry
      }
      iso_xy
    }, simplify = F)
    if(length(out) > 1) {out <- do.call(rbind, out)} else {out <- out[[1]]}
    ob_list[[i]] <- out
    t_list[[i]] <- unique(rwt)
  }

  ids <- paste0(rep(1:nrow(origin), sapply(ob_list, length)), "_", unlist(t_list), "-", unlist(t_list))
  out <- data.frame(feature = rep(1:nrow(origin), unlist(lapply(ob_list, length))), layer = unlist(t_list))
  st_geometry(out) <- Reduce(c, ob_list)
  rownames(out) <- ids

  # summarise and return
  return(vect(out))
}
