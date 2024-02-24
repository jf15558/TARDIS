#' random_walk
#'
#' Conduct random walks within a tardis graph. Full traverses in time and space
#' are possible, but as walks in both these dimensions not are not necessarily
#' easy to interpret, walks are restricted to the time layer in which they
#' originate as a default. Walk lengths can be specified by the number of
#' steps or a target distance to traverse. In the case of the latter, the
#' walk is terminated at the closest distance to the target as the walks will
#' not necessarily be able to precisely hit the target distance.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in tardis. All entries must be >= 0
#' and finite (NaN or Inf), or NA. Typically the output of weight_tgraph
#' @param origin A simple features collection produced by stp, denoting
#' the starting points of the random walks.
#' @param mode The nature of the random walk, either the number of 'steps' to
#' take or the 'cost' in terms of the cumulative weight of the edges traversed.
#' @param rwlen The length of the random walk as selected by @param mode, either
#' a single number or a vector with as many elements as pindex, to enable
#' different walk lengths for each point. When weights = NULL and mode = 'cost',
#' this  can be directly interpreted as metres covered by the random walk, but
#' may also be expressed in relation to a user-supplied weighting scheme.
#' @param restrict A logical indicating whether the random walk should be
#' constrained to the same time layer as its origin. TRUE by default.
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return If restrict = TRUE (default), then a SpatRaster with each stack
#' #' layer recording the frequencies of cell visits. If restrict = FALSE,
#' a list of SpatRasters, with each named layer recording the frequencies of
#' cells visited in each landscape layer through time
#' @import terra
#' @importFrom igraph graph_from_edgelist delete.edges random_walk get.edge.ids components
#' @importFrom methods as
#' @importFrom stats setNames
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
#' #foo <- random_walk(tardis = gt, weights = gtw, pts[3:4,], rwlen = 1e6)

random_walk <- function(tardis, weights = NULL, origin, mode = "steps", rwlen = 1000, restrict = TRUE, verbose = TRUE) {

  # tardis = test
  # weights = NULL
  # origin = pts[3:4,]
  # rwlen = 1000
  # mode = "steps"
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

  if(length(mode) != 1 | !is.atomic(mode) | !is.character(mode)) {
    stop("mode should be a single character string, one of 'steps' or 'cost'")
  }
  if(!mode %in% c("steps", "cost")) {
    stop("mode should be a single character string, one of 'steps' or 'cost'")
  }

  if(!is.atomic(rwlen) | !is.numeric(rwlen)) {
    stop("If not NULL, rwlen should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }
  if(length(rwlen) == 1) {
    rwlen <- rep(rwlen, nrow(origin))
  }
  if(length(rwlen) != length(rwlen)) {
    stop("If not NULL, rwlen should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }
  if(!all(is.finite(rwlen) | rwlen <= 0)) {
    stop("If not NULL, rwlen should be a single positive, finite numeric, or a vector of the same with as many rows as origin")
  }

  # get cell id from geographic position and age
  origin$cell <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(origin)) +
    (prod(tardis$gdat[1:2]) * (origin$bin - 1))

  # check point accessibility to ensure they come from the correct tardis object
  if(!all(origin$cell %in% tardis$edges[,1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }

  # get igraph version of graph (edge weights are conductive, rather than resistive as in cppRouting)
  if(verbose) {cat("Building graph\n")}
  tardis$edges <- tardis$edges[which(!is.na(weights)),]
  ig <- graph_from_edgelist(as.matrix(tardis$edges[,1:2]))
  igraph::E(ig)$weight <- 1 / weights
  true_w <- tardis$edges[which(!is.na(weights)),5]

  det_list <- list()
  for(i in 1:nrow(origin)) {

    if(verbose) {cat(paste0("Running random walks [", i, "/", nrow(origin), "]\r"))
      if(i == nrow(origin)) {cat("\n")}
    }

    if(restrict) {
      rng <- c(origin$cell[i] %/% prod(tardis$gdat[1:2]), (origin$cell[i] %/% prod(tardis$gdat[1:2]) + 1)) * prod(tardis$gdat[1:2])
      grp <- delete.edges(ig, which(tardis$edges[,1] < rng[1] + 1 | tardis$edges[,2] > rng[2]))
      tw <- true_w[-which(tardis$edges[,1] < rng[1] + 1 | tardis$edges[,2] > rng[2])]
    } else {
      grp <- ig
      tw <- true_w
    }
    if(mode == "cost") {
      steps <- ceiling((rwlen[i] / median(tw)) * 2)
    } else {
      steps <- rwlen[i]
    }

    fail <- 1
    while(fail) {
      rw <- as.vector(igraph::random_walk(grp, start = as.character(origin$cell[i]), steps = steps))
      if(mode == "cost") {
        dst <- cumsum(tw[get.edge.ids(grp, c(rw[1], rep(rw[2:(length(rw) - 1)], each = 2), rw[length(rw)]))])
        if(any(dst >= steps)) {fail <- 0}
        rw <- rw[1:which.min(abs(cumsum(tw[igraph::get.edge.ids(grp, c(rw[1], rep(rw[2:(length(rw) - 1)], each = 2), rw[length(rw)]))]) - steps))]
      } else {
        fail <- 0
      }
    }

    # convert to counts and rasterize
    rwt <- rw %/% prod(tardis$gdat[1:2]) + 1
    rwp <- rw %% prod(tardis$gdat[1:2])
    tmp <- rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8]))
    out <- tapply(rwp, rwt, function(x) {
      ob <- tmp
      vl <- table(x)
      ob[as.numeric(names(vl))] <- as.vector(vl)
      ob
    })
    out <- as(setNames(Reduce(c, out), names(out)), "SpatRaster")
    det_list[[i]] <- out
  }

  if(all(unlist(lapply(det_list, nlyr)) == 1)) {
    det_list <- as(setNames(Reduce(c, det_list), 1:nrow(origin)), "SpatRaster")
  }

  # summarise and return
  return(det_list)
}
