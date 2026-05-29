#' min_span
#'
#' Determine the least cost minimum spanning arborescence between a set of points
#' in a spatiotemporal landscape graph. Least cost paths are identified between
#' pairs of points using Djikstra's algorithm, then the minimum spanning
#' arborescence selected from these sets of pairwise distances.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param points `sf data.frame` A simple features collection produced by `point_check()`,
#' denoting the points to be linked with a minimum spanning arborescence.
#' @param verbose `logical` Should function progress be reported to the user?
#' @return An `sf data.frame` containing time-discrete linestrings that comprising the minimum
#' spanning arboresence for the set of input points.
#' @import terra sf cppRouting rlemon h3jsr
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph as_edgelist
#' @importFrom igraph E
#' @importFrom igraph delete_edges
#' @importFrom rlemon MinCostArborescence
#' @export
#'
#' @details
#' If geographic distances are used (default), then there will be a single least
#' cost path between any two points and the minimum spanning arborescence will be
#' equivalent to a minimum spanning tree. However, as TARDIS graphs are
#' directed, it is possible that the shortest path from A to B may not be the
#' same as the shortest path from B to A depending on the weighting scheme,
#' hence why minimum spanning arborescence is used as the more general case.
#' Note that the time directionality in the TARDIS graph will also affect the
#' minimum spanning arborescence structure.
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
#' hlcp <- min_span(htd, origin = hpts[1:2,], dest = hpts[3:4,])
#' rlcp <- min_span(rtd, origin = rpts[1:2,], dest = rpts[3:4,])
#' }

min_span <- function(tardis, weights = "gdist", points, verbose = TRUE) {

  # tardis = rtd
  # weights = "gdist"
  # points = pt
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

  if (!class(points)[1] == c("sf")) {
    stop("Supply origin as the output of stp")
  }

  if(!is.na(tardis$gdat[7])) {

    points$cell <- match(point_to_cell(points, tardis$gdat[7]), grid) +
      (tardis$gdat[5] * (points$bin - 1))

  } else {

    points$cell <- cellFromXY(samprast, st_coordinates(points)) +
      (tardis$gdat[5] * (points$bin - 1))

  }

  if (!all(points$cell %in% tardis$edges[, 1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }
  if (verbose) {
    cat("Initialising graph\n")
  }
  tardis <- instantiate_tardis(tardis = tardis, weights = weights)

  if (verbose) {
    cat("Running paths\r")
  }

  dists <- get_distance_matrix(tardis$tgraph, points$cell, points$cell)
  # set impossible paths to infinite
  dists[is.na(dists)] <- Inf

  gr <- graph_from_adjacency_matrix(dists, weighted = T)
  # drop inaccessible cells
  gr <- delete_edges(gr, which(E(gr)$weight == Inf))
  el <- as_edgelist(gr)
  el <- matrix(match(el, as.character(points$cell)), ncol = 2)
  ws <- E(gr)$weight

  mca <- MinCostArborescence(el[,1], el[,2], ws, el[1,1], nrow(points))

  srt <- points$cell[mca$sources]
  end <- points$cell[mca$targets]

  srt_t <- srt %/% tardis$gdat[5] + 1
  end_t <- end %/% tardis$gdat[5] + 1

  paths <- lapply(suppressMessages(get_path_pair(tardis$tgraph, from = srt, to = end)),
                  function(x) {as.numeric(x)})

  costs <- lapply(1:length(paths), function(x) {
    if (verbose) {
      cat(paste0("Running paths [", x, "/", length(paths),
                 "]\r"))
    }
    pth <- paths[[x]]
    wp <- c(get_distance_pair(tardis$tgraph, pth[-length(pth)],
                              pth[-1], algorithm = "Djikstra"), 0)
    if (weights != "tdist") {
      tp <- c(get_distance_pair(tardis$tgraph, pth[-length(pth)],
                                pth[-1], algorithm = "Djikstra", aggregate_aux = T),
              0)
    }
    else {
      tp <- wp
    }
    if (verbose) {
      if (x == length(paths)) {
        cat("\n")
      }
    }
    list(wp, tp)
  })
  wpaths <- lapply(costs, function(x) {
    x[[1]]
  })
  tpaths <- lapply(costs, function(x) {
    x[[2]]
  })

  pords <- ifelse(srt_t - end_t <= 0, 1, 0)
  if (verbose) {
    cat("Summarising\n")
  }

  path_groups <- path_ids <- wvec <- tvec <- orgs <- dsts <- list()
  for (i in 1:length(paths)) {

    if(!is.na(tardis$gdat[7])) {
      path_xy <- st_coordinates(cell_to_point(grid[paths[[i]] %% tardis$gdat[5]]))

    } else {
      path_xy <- xyFromCell(samprast, paths[[i]] %% tardis$gdat[5])
    }

    pbin <- ceiling(paths[[i]] / tardis$gdat[5])
    pseq <- rep(1:length(rle(pbin)$length), rle(pbin)$length)
    wcst <- rev(rev(c(rbind(tapply(wpaths[[i]], INDEX = pseq,
                                   sum), 0)))[-1])
    dcst <- rev(rev(c(rbind(tapply(tpaths[[i]], INDEX = pseq,
                                   sum), 0)))[-1])
    pgrp <- list()
    within <- tapply(1:length(pseq), INDEX = pseq, function(y) {
      if (length(y) == 1) {
        st_linestring(path_xy[c(y, y), ])
      }
      else {
        st_linestring(path_xy[y, ])
      }
    })
    for (j in 1:length(seq(1, length(wcst), 2))) {
      pgrp[[seq(1, length(wcst), 2)[j]]] <- within[[j]]
    }
    pids <- paste0(i, "_", pbin[1], "-", pbin[1])
    if (length(unique(pbin)) > 1) {
      between <- lapply(which(diff(pseq) != 0) + 1, function(y) {
        st_linestring(as.matrix(path_xy[c(y - 1, y),
                                        , drop = FALSE]))
      })
      for (j in 1:length(seq(2, length(wcst), 2))) {
        pgrp[[seq(2, length(wcst), 2)[j]]] <- between[[j]]
      }
      pids <- pbin[c(1, which(diff(pbin) != 0) + 1)]
      pids <- rev(rev(c(t(cbind(paste0(i, "_", pids, "-",
                                       pids), c(paste0(i, "_", pids[-length(pids)],
                                                       "-", pids[-1]), NA)))))[-1])
    }
    pord <- if (pords[i] == 1) {
      1:length(pgrp)
    } else {
      rev(1:length(pgrp))
    }
    path_groups[[i]] <- pgrp[pord]
    path_ids[[i]] <- pids[pord]
    wvec[[i]] <- wcst[pord]
    tvec[[i]] <- dcst[pord]
  }
  path_ids <- unlist(path_ids)
  ob <- cbind.data.frame(feature = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                     "_"), function(y) {
                                                                       y[[1]]
                                                                     }))), srt_bin = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                                                                       "_|-"), function(y) {
                                                                                                                         y[[2]]
                                                                                                                       }))), end_bin = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                                                                                                                         "_|-"), function(y) {
                                                                                                                                                                           y[[3]]
                                                                                                                                                                         }))), distance = unlist(tvec), cost = unlist(wvec))
  ob$layer <- (ob$srt_bin + ob$end_bin) / 2
  ob <- ob[,c("feature", "layer")]
  st_geometry(ob) <- st_sfc(unlist(path_groups, recursive = F),
                            crs = "+proj=longlat")
  return(ob)
}
