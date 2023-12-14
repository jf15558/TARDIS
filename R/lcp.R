#' lcp
#'
#' Determine least cost paths through a spatiotemporal landscape graph. Least
#' cost paths are identified between pairs of origin and destination coordinates
#' using Djikstra's algorithm. Costs are either geographic distances or a custom
#' weighting scheme supplied by the user.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in tardis. All entries must be >= 0
#' and finite (NaN or Inf), or NA. Typically the output of weight_tgraph
#' @param origin A simple features collection produced by stp, denoting
#' the starting points of the least cost paths
#' @param dest As for origin. The user should be careful to the time ordering
#' of point pairs matches the time linking mode if tardis contains multiple
#' layers. For example, points in younger layers (higher bin numbers) will not
#' be accessible from points in older layers (lower bin numbers) if the linking
#' mode of tardis is backwards in time (tlink = 2)
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return A list with three elements: 1. 'paths', a simple features linestring
#' geometry collection representing the time-discrete portions of each least
#' cost path; 2. 'costs', a data.frame containing the geographic and weighted
#' costs of each least cost path (identical when weights = NULL); 3. 'cells',
#' a list of cell IDs comprising each least cost path.
#' @import terra sf cppRouting
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#'
#' #data("galapagos")
#' #gal <- unwrap(galapagos)
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

lcp <- function(tardis, weights = NULL, origin, dest, verbose = TRUE) {

  # tardis = test2
  # weights = test_weight
  # origin = pts[1:2,]
  # dest = pts[3:4,]
  # verbose = TRUE

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
  if(!class(dest)[1] == c("sf")) {
    stop("Supply dest as the output of stp")
  }
  if(nrow(origin) != nrow(dest)) {
    stop("The number of origin and destination points should be the same (i.e. paired points")
  }

  # get cell id from geographic position and age
  origin$cell <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(origin)) +
    (prod(tardis$gdat[1:2]) * (origin$bin - 1))

  dest$cell <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(dest)) +
    (prod(tardis$gdat[1:2]) * (dest$bin - 1))

  # check point accessibility to ensure they come from the correct tardis object
  if(!all(origin$cell %in% tardis$edges[,1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }
  if(!all(dest$cell %in% tardis$edges[,1])) {
    stop("One or more points in dest do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for dest")
  }
  if(any(dest$bin > origin$bin) & tardis$link.mode[2] == 2) {
    stop("tardis is linked backwards in time, but some paths terminate in younger layers than their starting point")
  }
  if(any(dest$bin < origin$bin) & tardis$link.mode[2] == 1) {
    stop("tardis is linked forwards in time, but some paths terminate in older layers than their starting point")
  }

  # initialise graph
  if(verbose) {cat("Initialising graph\n")}
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

  if(verbose) {cat("Running paths\n")}
  paths <- lapply(suppressMessages(get_path_pair(tardis$tgraph, from = origin$cell, to = dest$cell)), as.numeric)
  wpaths <- lapply(paths, function(x) {c(get_distance_pair(tardis$tgraph, x[-length(x)], x[-1]), 0)})
  tpaths <- lapply(paths, function(x) {c(get_distance_pair(tardis$tgraph, x[-length(x)], x[-1], aggregate_aux = T), 0)})
  pords <- ifelse(origin$bin - dest$bin <= 0, 1, 0)

  if(verbose) {cat("Summarising\n")}
  path_groups <- path_ids <- wvec <- tvec <- orgs <- dsts <- list()
  for(i in 1:length(paths)) {

    # cells to lines
    path_xy <- xyFromCell(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), paths[[i]] %% prod(tardis$gdat[1:2]))
    pbin <- ceiling(paths[[i]] / prod(tardis$gdat[1:2]))
    pseq <- rep(1:length(rle(pbin)$length), rle(pbin)$length)
    wcst <- rev(rev(c(rbind(tapply(wpaths[[i]], INDEX = pseq, sum), 0)))[-1])
    dcst <- rev(rev(c(rbind(tapply(tpaths[[i]], INDEX = pseq, sum), 0)))[-1])

    pgrp <- list()
    within <- tapply(1:length(pseq), INDEX = pseq, function(y) {
      if(length(y) == 1) {st_linestring(path_xy[c(y, y),])} else {st_linestring(path_xy[y,])}
    })
    for(j in 1:length(seq(1, length(wcst), 2))) {pgrp[[seq(1, length(wcst), 2)[j]]] <- within[[j]]}
    pids <- paste0(i, "_", pbin[1], "-", pbin[1])

    if(length(pbin) > 1) {
      between <- lapply(which(diff(pseq) != 0) + 1, function(y) {st_linestring(as.matrix(path_xy[c(y - 1, y), ,drop = FALSE]))})
      for(j in 1:length(seq(2, length(wcst), 2))) {pgrp[[seq(2, length(wcst), 2)[j]]] <- between[[j]]}
      pids <- pbin[c(1, which(diff(pbin) != 0) + 1)]
      pids <- rev(rev(c(t(cbind(paste0(i, "_", pids, "-", pids), c(paste0(i, "_", pids[-length(pids)], "-", pids[-1]), NA)))))[-1])
    }

    pord <- if(pords[i] == 1) {1:length(pgrp)} else {rev(1:length(pgrp))}
    path_groups[[i]] <- pgrp[pord]
    path_ids[[i]] <- pids[pord]
    wvec[[i]] <- wcst[pord]
    tvec[[i]] <- dcst[pord]
  }

  # restructure and return
  path_ids <- unlist(path_ids)
  ob <- cbind.data.frame(path = as.numeric(unlist(lapply(strsplit(path_ids, "_"), function(y) {y[[1]]}))),
                         bin = as.numeric(unlist(lapply(strsplit(path_ids, "-"), function(y) {y[[2]]}))),
                         order = unlist(lapply(path_groups, function(x) {1:length(x)})),
                         cost = unlist(wvec), distance = unlist(tvec))
  row.names(ob) <- path_ids
  st_geometry(ob) <- st_sfc(unlist(path_groups, recursive = F), crs = "+proj=longlat")
  return(ob)
}
