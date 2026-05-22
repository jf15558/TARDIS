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
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param origin `sf data.frame` A simple features collection produced by `point_check()`,
#' denoting the starting points for the random walks.
#' @param mode `character`. The nature of the random walk, either the number of 'steps' to
#' take or the 'cost' in terms of the cumulative weight of the edges traversed.
#' @param rwlen The length of the random walk as selected by @param mode, either
#' a single number or a vector with as many elements as origin points, to enable
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
#' gal <- crop(gal, extent(-92, -88, -2, 1))
#' gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1),
#'                                     ncol = 3, byrow = T), right = F)
#' gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
#'
#' vars = list(elev = classify(gal, cbind(-Inf, 0, 0)))
#' gtw <- weight_tardis(test2, vars = vars,
#'                      mfun = function(origin, dest, lnum, ...) {
#'                                (origin$hdist^2 + abs(origin$vdist)^2) * 10})
#'
#' org <- rbind(c(-89, -1.05, 2), c(-89.5, -0.7, 2))
#' dst <- rbind(c(-91.2, -1, 0), c(-91.6, -0.4, 0))
#' pts <- stp(test2, rbind(org, dst))
#'
#' foo <- random_walk(tardis = gt, weights = gtw, pts[3:4,], rwlen = 1e6)
#' }

random_walk <- function(tardis, weights = "gdist", origin, mode = "steps", rwlen = 1000, restrict = TRUE, verbose = TRUE) {

  # tardis = rtd
  # weights = "gdist"
  # origin = rpt[3:4,]
  # rwlen = 100000
  # mode = "cost"
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
    while(fail) {
      rw <- as.vector(igraph::random_walk(grp, start = as.character(origin$cell[i]), steps = steps))
      if(mode == "cost") {
        dst <- cumsum(tw[get.edge.ids(grp, c(rw[1], rep(rw[2:(length(rw) - 1)], each = 2), rw[length(rw)]))]) - rwlen[i]
        if(any(dst >= 0)) {fail <- F}
        rw <- rw[1:which(dst > 0)[1]]
      } else {
        fail <- F
      }
    }

    # convert to counts and rasterize
    rw <- table(rw)
    rwt <- as.numeric(names(rw)) %/% tardis$gdat[5] + 1
    rwp <- as.numeric(names(rw)) %% tardis$gdat[5]

    out <- lapply(unique(rwt), function(x) {

      freq <- rw[which(rwt == x)]
      if(!is.na(tardis$gdat[7])) {
        cls <- st_wrap_dateline(cell_to_polygon(grid[rwp[which(rwt == x)]]), options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
        freqs <- cbind.data.frame(point = i, bin = x, freq = as.vector(rw[which(rwt == x)]))
        st_geometry(freqs) <- cls
        colnames(freqs) <- c("point", "bin", "freq", "geometry")

      } else {
        freqs <- samprast
        freqs[] <- NA
        freqs[][rwp[which(rwt == x)],] <- rw[which(rwt == x)]
        names(freqs) <- "freq"
        freqs <- st_as_sf(as.polygons(freqs, aggregate = F))
        freqs <- st_sf(cbind.data.frame(point = i, freqs))
      }
      freqs
    })
    det_list[[i]] <- out
  }
  # summarise and return
  return(do.call(rbind, unlist(det_list, recursive = F)))
}
