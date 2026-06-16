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
#' @param origin `SpatVector`. The output of `point_check()`, denoting the
#' starting points of the least cost paths.
#' @param dest `SpatVector`. The output of `point_check()`, denoting the
#' destination points of the least cost paths. The user should be careful to ensure
#' that the time ordering of point pairs matches the time linking mode if
#' tardis contains multiple layers (i.e., points in older layers cannot be
#' accessed from points in younger layers if the linking mode is forwards in time
#' (`tlink = 1`).
#' @param verbose `logical`. Should function progress be reported to the user?
#' @return An `SpatVector` of time-discrete lines representing the least
#' cost paths between each point pair, recording which overall path they belong
#' to (`$feature`), the costs along each line (`$cost`) and their geographic
#' distances (`$distance`, identical if `weights = gdist`).
#' @import terra sf cppRouting h3jsr
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
#' rasts <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#' rlink <- link_islands(rasts)
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#' pts <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-88.70836, -0.2627832, 2))
#'
#' rpts <- point_check(rtd, pts)
#' rlcp <- least_cost(rtd, origin = pts[1,], dest = pts[2,])
#' }

least_cost <- function(tardis, weights = "gdist", origin, dest, verbose = TRUE) {

   # tardis = rtd
   # weights = "gdist"
   # origin = pts[phy$edge[410,1],]
   # dest = pts[phy$edge[410,2],]
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

  if (!class(origin)[1] == c("SpatVector")) {
    stop("Supply origin as the output of point_check")
  }
  if (!class(dest)[1] == c("SpatVector")) {
    stop("Supply dest as the output of point_check")
  }
  if (nrow(origin) != nrow(dest)) {
    stop("The number of origin and destination points should be the same (i.e. paired points")
  }

  if(!is.na(tardis$gdat[7])) {

    pcell <- suppressMessages(point_to_cell(crds(origin), tardis$gdat[7]))
    origin$cell <- match(pcell, grid) + (tardis$gdat[5] * (origin$layer - 1))
    pcell <- suppressMessages(point_to_cell(crds(dest), tardis$gdat[7]))
    dest$cell <- match(pcell, grid) + (tardis$gdat[5] * (dest$layer - 1))

  } else {

    origin$cell <- cellFromXY(samprast, crds(origin)) +
      (tardis$gdat[5] * (origin$layer - 1))
    dest$cell <- cellFromXY(samprast, crds(dest)) +
      (tardis$gdat[5] * (dest$layer - 1))

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
  pords <- ifelse(origin$layer - dest$layer <= 0, 1, 0)
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
  ob <- ob[,c("feature", "layer", "cost", "distance")]
  st_geometry(ob) <- st_sfc(unlist(path_groups, recursive = F),
                            crs = "+proj=longlat")
  ob <- st_wrap_dateline(ob, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  return(vect(ob))
}
