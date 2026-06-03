#' rand_walk
#'
#' Conduct random walks within a tardis graph. Full traverses in time and space
#' are possible, but as walks in both these dimensions not are not necessarily
#' easy to interpret, walks are restricted to the time layer in which they
#' originate as a default. Walk lengths can be specified by the number of
#' steps or a target distance to traverse. In the case of the latter, the
#' walk is terminated at the closest distance to the target as the walks will
#' not necessarily be able to precisely hit the target distance.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param origin `SpatVector`. The output of `point_check()`, denoting the
#' starting points for the random walks.
#' @param mode `character`. The nature of the random walk, either the number of 'steps' to
#' take or the 'cost' in terms of the cumulative weight of the edges traversed.
#' @param rwlen The length of the random walk as selected by @param mode, either
#' a single number or a vector with as many elements as origin points, to enable
#' different walk lengths for each point. When weights = NULL and mode = 'cost',
#' this  can be directly interpreted as metres covered by the random walk, but
#' may also be expressed in relation to a user-supplied weighting scheme.
#' @param restrict A logical indicating whether the random walk should be
#' constrained to the same time layer as its origin. TRUE by default.
#' @param trials An integer indicating the number of random walk trials to
#' conduct before quitting, if mode = "cost" and the random walk length proves
#' insufficient to reach to cost value.
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return If restrict = TRUE (default), then a SpatRaster with each stack
#' #' layer recording the frequencies of cell visits. If restrict = FALSE,
#' a list of SpatRasters, with each named layer recording the frequencies of
#' cells visited in each landscape layer through time
#' @import terra h3jsr
#' @importFrom igraph graph_from_edgelist delete_edges random_walk get.edge.ids components
#' @importFrom methods as
#' @importFrom stats setNames
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
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#' org <- rbind(c(-89, -1.05, 2))
#' pts <- point_check(htd, org)
#'
#' foo <- rand_walk(htd, origin = pts, rwlen = 1e4)
#' }

rand_walk <- function(tardis, weights = "gdist", origin, mode = "steps", rwlen = 1000, restrict = TRUE, trials = 100, verbose = TRUE) {

   #tardis = rtd
   #weights = "gdist"
   #origin = pt
   #rwlen = 1000
   #mode = "steps"
   #restrict = T
   #verbose = T
   #trials = 100

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
    stop("Supply origin as the output of point_check()")
  }

  if(!is.na(tardis$gdat[7])) {

    pcell <- suppressMessages(point_to_cell(crds(origin), tardis$gdat[7]))
    origin$cell <- match(pcell, grid) + (tardis$gdat[5] * (origin$layer - 1))

  } else {

    origin$cell <- cellFromXY(samprast, crds(origin)) +
      (tardis$gdat[5] * (origin$layer - 1))

  }

  if (!all(origin$cell %in% tardis$edges[, 1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
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

  if(!is.atomic(trials) | !is.numeric(trials)) {
    stop("trials should be a single positive integer")
  }
  if(length(trials) != 1) {
    stop("trials should be a single positive integer")
  }

  # get igraph version of graph (edge weights are conductive, rather than resistive as in cppRouting)
  if(verbose) {cat("Building graph\n")}

  tardis$edges <- tardis$edges[which(!is.na(tardis$edges[,weights])),]
  ig <- graph_from_edgelist(as.matrix(tardis$edges[,1:2]))
  igraph::E(ig)$weight <- 1 / tardis$edges[,weights]
  true_w <- tardis$edges[,5]

  det_list <- list()
  for(i in 1:nrow(origin)) {

    if(verbose) {cat(paste0("Running random walks [", i, "/", nrow(origin), "]\r"))
      if(i == nrow(origin)) {cat("\n")}
    }

    if(restrict) {
      rng <- c(origin$cell[i] %/% tardis$gdat[5], (origin$cell[i] %/% tardis$gdat[5] + 1)) * tardis$gdat[5]
      grp <- delete_edges(ig, which(tardis$edges[,1] < rng[1] + 1 | tardis$edges[,2] > rng[2]))
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

    fail <- T
    iter <- 1
    while(fail) {
      rw <- as.vector(igraph::random_walk(grp, start = origin$cell[i], steps = steps))
      iter <- iter + 1
      if(mode == "cost") {
        dst <- cumsum(tw[get.edge.ids(grp, c(rw[1], rep(rw[2:(length(rw) - 1)], each = 2), rw[length(rw)]))]) - rwlen[i]
        if(any(dst >= 0)) {fail <- F}
        rw <- rw[1:which(dst > 0)[1]]
      } else {
        fail <- F
      }
      if(iter > trials) {
        stop(paste0("Random walk failed to reach the requested cost after ", trials, " trials"))
      }
    }

    # convert to counts and rasterize
    rw <- table(rw)
    rwt <- as.numeric(names(rw)) %/% tardis$gdat[5] + 1
    rwp <- as.numeric(names(rw)) %% tardis$gdat[5]

    out <- lapply(unique(rwt), function(x) {

      freq <- rw[which(rwt == x)]
      if(!is.na(tardis$gdat[7])) {
        freqs <- vect(st_wrap_dateline(cell_to_polygon(grid[rwp[which(rwt == x)]]), options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")))
        freqs$point <- i
        freqs$bin <- x
        freqs$freq <- as.vector(rw[which(rwt == x)])

      } else {
        freqs <- samprast
        freqs[] <- NA
        freqs[][rwp[which(rwt == x)],] <- rw[which(rwt == x)]
        names(freqs) <- "freq"
        freqs <- as.polygons(freqs, aggregate = F)
        freqs$point <- i
        freqs$bin <- x
        freqs[,c(2, 3, 1)]
      }
      freqs
    })
    det_list[[i]] <- out
  }
  # summarise and return
  return(svc(unlist(det_list, recursive = F)))
}
