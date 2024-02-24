#' isochrone
#'
#' Calculate a cost isochrone around a point, that is all cells which fall
#' within a certain cost of access of that point, based on the cumulative
#' cost of traversing to those surrounding cells. Isochrones are calculated
#' within the time-specific landscape layer for a given point by default,
#' rather than also extending along temporal links. Full space-time
#' isochrones can be obtained, but it is not necessarily obvious how the cost
#' should ideally change with temporal distance from the origin.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in tardis. All entries must be >= 0
#' and finite (NaN or Inf), or NA. Typically the output of weight_tgraph
#' @param origin A simple features collection produced by stp, denoting
#' the starting points of the random walks.
#' @param cost The maximum cumulative cost of below which cells will be
#' included within the isochrone in a point's time-specific origin landscape
#' layer, either a single number or a vector with as many elements as pindex,
#' to enable different costs for each point. Cost will increase proportionally
#' in successive layers away from the origin as described above.
#' @param restrict A logical indicating whether the isochrone should be
#' constrained to the same time layer as its origin. TRUE by default.
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return If restrict = TRUE (default), then a SpatRaster with each layer
#' recording the frequencies of cell visits. If restrict = FALSE, a list of
#' SpatRasters, with each named layer recording the frequencies of cells visited
#' in each landscape layer through time
#' @import terra pbapply parallel
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#'
#' #gal <- galapagos()
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1),
#' #                                    ncol = 3, byrow = T), right = F)
#' #gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
#'
#' #vars = list(elev = classify(gal, cbind(-Inf, 0, 0)))
#' #gtw <- weight_tardis(test2, vars = vars,
#' #                     mfun = function(origin, dest, lnum, ...) {
#' #                               (origin$hdist^2 + abs(origin$vdist)^2) * 10})
#'
#' #org <- rbind(c(-89, -1.05, 2), c(-89.5, -0.7, 2))
#' #dst <- rbind(c(-91.2, -1, 0), c(-91.6, -0.4, 0))
#' #pts <- stp(test2, rbind(org, dst))
#'
#' #foo <- lcp(tardis = gt, weights = gtw, pts[1:2,], pts[3:4,])
#' #foo2 <- isochrone(tardis = gt, weights = gtw, foo)

isochrone <- function(tardis, weights = NULL, origin, cost = 1e5, restrict = TRUE, verbose = TRUE) {

  # tardis = test
  # weights = test_weight
  # origin = pts[3:4,]
  # cost = 1e4
  # restrict = T
  # verbose = T

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }

  if(!is.null(weights)) {

    if(!is.atomic(weights)) {
      stop("weights must be a vector")
    }
    if(!is.numeric(weights)) {
      stop("weights (the vector or chosen column) must be numeric")
    }
    if(any(is.infinite(weights) | weights < 0 | is.nan(weights))) {
      stop("All entries in weights must be positive, finite numerics")
    }
    if(length(weights) != nrow(tardis$edges)) {
      stop("weights must be the same length as the number of edges in tardis")
    }
  } else {
    weights <- tardis$edges[,5]
  }

  if(!class(origin)[1] == c("sf")) {
    stop("Supply origin as the output of stp")
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

  # get cell id from geographic position and age
  origin$cell <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(origin)) +
    (prod(tardis$gdat[1:2]) * (origin$bin - 1))

  # check point accessibility to ensure they come from the correct tardis object
  if(!all(origin$cell %in% tardis$edges[,1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }

  # initialise graph
  if(any(is.na(weights))) {
    tardis$tgraph$src <- tardis$tgraph$src[which(!is.na(weights))]
    tardis$tgraph$dst <- tardis$tgraph$dst[which(!is.na(weights))]
    tardis$edges <- tardis$edges[which(!is.na(weights)),]
    tardis$tgraph$dict <-  tardis$tgraph$dict[which(tardis$tgraph$dict$ref %in% c(tardis$tgraph$src, tardis$tgraph$dst))]
    tardis$tgraph$id <- 0:(length(tardis$tgraph$dict$ref) - 1)
    tardis$tgraph$nbnode <- length(tardis$tgraph$id)
  }
  tardis$tgraph$attrib$aux <- tardis$edges[!is.na(weights),5]
  tardis$edges[,5] <- weights[!is.na(weights)]
  tardis$tgraph$data <- data.frame(from = tardis$tgraph$dict$id[match(tardis$tgraph$src, tardis$tgraph$dict$ref)],
                                   to   = tardis$tgraph$dict$id[match(tardis$tgraph$dst, tardis$tgraph$dict$ref)],
                                   dist = tardis$edges[,5])
  tardis$tgraph <- tardis$tgraph[1:5]

  ob_list <- t_list <- list()
  for(i in 1:nrow(origin)) {

    if(verbose) {cat(paste0("Running isochrones [", i, "/", nrow(origin), "]\r"))
      if(i == nrow(origin)) {cat("\n")}
    }
    res <- as.numeric(get_isochrone(tardis$tgraph, from = origin$cell[i], lim = cost[i])[[1]])
    rwt <- res %/% prod(tardis$gdat[1:2]) + 1
    rwp <- res %% prod(tardis$gdat[1:2])
    if(restrict) {
      rwp <- rwp[which(rwt == origin$cell[i] %/% prod(tardis$gdat[1:2]) + 1)]
      rwt <- rwt[which(rwt == origin$cell[i] %/% prod(tardis$gdat[1:2]) + 1)]
    }

    tmp <- rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8]))
    out <- tapply(rwp, rwt, function(x) {
      ob <- tmp
      ob[x] <- 1
      ob
      st_as_sf(terra::as.polygons(ob))$geometry
    })
    ob_list[[i]] <- out
    t_list[[i]] <- unique(rwt)
  }

  ids <- paste0(rep(1:nrow(origin), unlist(lapply(ob_list, length))), "_", unlist(t_list), "-", unlist(t_list))
  out <- data.frame(path = rep(1:nrow(origin), unlist(lapply(ob_list, length))), bin = unlist(t_list))
  st_geometry(out) <- st_sfc(unlist(ob_list, recursive = F))
  rownames(out) <- ids

  # summarise and return
  return(out)
}
