#' least_cost
#'
#' Determine least cost paths between pairs of origin and destination coordinates
#' using Djikstra's algorithm. Costs are either geographic distances or a custom
#' weighting scheme supplied by the user.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @param origin `sf data.frame`. A simple features collection produced by `point_check()`, denoting
#' the starting points of the least cost paths.
#' @param dest `sf data.frame`. As for origin. The user should be careful to the time ordering
#' of point pairs matches the time linking mode if tardis contains multiple
#' layers (i.e., points in older layers cannot be accessed from points in younger
#' layers if the linking mode is forwards in time (`tlink = 1`).
#' @param verbose `logical`. Should function progress be reported to the user?
#' @return An `sf data.frame` of time-discrete linestrings representing the least
#' cost paths between each point pair, along with the path costs and geographic
#' distances (identical if `weights = gdist`).
#' @import terra sf cppRouting h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
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
#' hlcp <- least_cost(htd, origin = hpts[1:2,], dest = hpts[3:4,])
#' rlcp <- least_cost(rtd, origin = rpts[1:2,], dest = rpts[3:4,])
#' }

least_cost <- function(tardis, weights = "gdist", origin, dest, verbose = TRUE) {

  # tardis = rtd
  # weights = "gdist"
  # origin = pt[1,]
  # dest = pt[2,]
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
  if (!class(dest)[1] == c("sf")) {
    stop("Supply dest as the output of stp")
  }
  if (nrow(origin) != nrow(dest)) {
    stop("The number of origin and destination points should be the same (i.e. paired points")
  }

  if(!is.na(tardis$gdat[7])) {

    pcell <- point_to_cell(origin, tardis$gdat[7])
    origin$cell <- match(pcell, grid) + (tardis$gdat[5] * (origin$bin - 1))
    pcell <- point_to_cell(dest, tardis$gdat[7])
    dest$cell <- match(pcell, grid) + (tardis$gdat[5] * (dest$bin - 1))

  } else {

    origin$cell <- cellFromXY(samprast, st_coordinates(origin)) +
      (tardis$gdat[5] * (origin$bin - 1))
    dest$cell <- cellFromXY(samprast, st_coordinates(dest)) +
      (tardis$gdat[5] * (dest$bin - 1))

  }

  if (!all(origin$cell %in% tardis$edges[, 1])) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }
  if (!all(dest$cell %in% tardis$edges[, 1])) {
    stop("One or more points in dest do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for dest")
  }
  if (any(dest$bin > origin$bin) & tardis$tlink == 2) {
    stop("tardis is linked backwards in time, but some paths terminate in younger layers than their starting point")
  }
  if (any(dest$bin < origin$bin) & tardis$tlink == 1) {
    stop("tardis is linked forwards in time, but some paths terminate in older layers than their starting point")
  }
  if (verbose) {
    cat("Initialising graph\n")
  }
  tardis <- instantiate_tardis(tardis, weights)

  if (verbose) {
    cat("Running paths\r")
  }
  paths <- lapply(suppressMessages(get_path_pair(tardis$tgraph,
                                                 from = origin$cell, to = dest$cell)), as.numeric)
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
  pords <- ifelse(origin$bin - dest$bin <= 0, 1, 0)
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

    pbin <- (paths[[i]] %/% tardis$gdat[5]) + 1
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
