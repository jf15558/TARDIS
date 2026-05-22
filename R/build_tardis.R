#' build_tardis
#'
#' Generate a TARDIS graph from a geoglist, representing landscape connectivity
#' across space and through time. Weights in the graph represent the geographic
#' distances between cells. As such, the function assumes that the input
#' `geoglist` contains topographic and/or bathymetric data measured in metres.
#' Connections through time can be spatially constant, or variable.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param times `numeric` or `NULL`. If there is only one layer in geog, then
#' times is not used. Otherwise, a  vector with `nlayers(geog) + 1` positive
#' elements, expressing the temporal boundaries of each layer as time in the past.
#' The vector need not end in the present (i.e. `0`), but time must flow from
#' oldest to youngest.
#' @param tlink `integer`. The linking mode between layers, either `1` (forwards-in-time),
#' `2` (backwards-in-time) or `3` (bidirectional). The forwards-in-time case is
#' the default.
#' @param island.check `logical`. If `geog` already contains links, then this argument
#' will be ignored. Otherwise should the layers of `geog` be checked for
#' islands and connected using `link_islands()`?
#' @param klink `integer`. The number of island connections to generate, as
#' called by `link_islands()`.
#' @param rotations `list` or `NULL`. By default `NULL`, indictating that temporal
#' links are spatially constant. Otherwise, a list with `nlayers(geog) - 1` elements
#' recording the shift in cell locations between layers (see @details).
#' @param verbose `logical`. Should function progress be to the user? This may
#' be useful when dealing with large rasters (high resolution and/or many layers).
#' @return A `tardis` spatiotemporal object with elements `edges`, `gdat`, `tdat`
#' and `tlink`. `edges` is the main graph object (see @details), while the other elements record
#' the spatial, temporal and linkage properties of the graph for internal use by
#' downstream TARDIS functions.
#' @import terra h3jsr sf
#' @importFrom geosphere distGeo bearing
#' @importFrom stats complete.cases
#' @export
#'
#' @details
#' The resulting graph can be thought of as a 3D lattice. Graph weights for
#' horizontal edges within layers record the great circle distances between
#' adjacent cells in metres, adjusted' for differences in cell elevation using
#' Pythagoras's theorem. For raster inputs, 8-degree (Queen's case) adjacency
#' is used, while for hexagonally resampled inputs, 6-degree adjacency is used.
#' Single layer cases (i.e. no time element) are allowed. Otherwise vertical
#' edges between time layers are assigned weights of zero so that they do not
#' affect downstream distance calculations.
#'
#' In many cases, the positions of cells will remain constant within the extent
#' of each layer in `geog`. For global landscapes over geological timescales,
#' however, the positions of landmasses chance noticeably due to continental
#' drift. In these cases, the edges can be altered so that they connect pairs of
#' geographically homologous cells through time. This is implemented using the
#' `rotations` argument.
#'
#' If rotations is not `NULL`, then each element in the list is a two-column
#' numeric matrix containing the IDs of geographically homologous cells between
#' successive pairs of layers in x. No list element can be blank as there must
#' be at least one edge between each pair of layers. Otherwise, any number of
#' edges between layers can be specified, although it  makes sense for there to
#' be maximally as many edges as cells within a layer. As this is quite a specific
#' format, the function `get_rotations` can be used to expedite its construction.
#' Examples of the rotations object are available via its documentation
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # create a geoglist with hexagonal resampling and mask the sea
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' # build a tardis with hexagonal cells
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # create a geoglist in raster format and mask the sea
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#'
#' # build a tardis from raster cells
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#' }

build_tardis <- function(geog, times = NULL, tlink = 1, island.check = TRUE, klink = NULL, rotations = NULL, verbose = TRUE) {

  # geog = rasts
  # #times = c(seq(2.25, 0, -0.5), 0)
  # times = c(117, 114, 112)
  # tlink = 1
  # island.check = F
  # klink = 2
  # rotations = NULL
  # verbose = TRUE

  nlayers <- length(geog$layers)
  if(is.na(geog$gdat[7])) {nlayers <- nlyr(geog$layers)}

  if(nlayers > 1) {
    if (!exists("times")) {
      stop("If there are multiple layers in geog, then times must be specified")
    }
    if (!is.numeric(times) | length(times) != nlayers + 1) {
      stop("Please supply times as a vector of time bin boundaries with n elements in geog$layers + 1")
    }
    if (any(diff(times) > 0)) {
      stop("All elements of times should be positive (i.e. before present) and in descending age order")
    }
    if (!is.numeric(tlink) | length(tlink) != 1) {
      stop("tlink should be one of 1 (forward-in-time), 2 (backward-in-time) or 3 (bidirectional")
    }
    if (!tlink %in% c(1:3)) {
      stop("tlink should be one of 1 (forward-in-time), 2 (backward-in-time) or 3 (bidirectional")
    }
    if (is.null(rotations)) {
      rotations <- lapply(1:(nlayers - 1), function(x) {matrix(rep(1:geog$gdat[5], 2), ncol = 2)})
    } else {
      if (!is.list(rotations)) {
        stop("Rotations should be a list")
      }
      if (!all(unlist(lapply(rotations, function(x) {class(x)[1]})) == "matrix")) {
        stop("All elements of rotations should be two-column matrices")
      }
      if (!all(unlist(lapply(rotations, function(x) {mode(x)[1]})) == "numeric")) {
        stop("All elements of rotations should be two-column matrices")
      }
      if (!all(unlist(lapply(rotations, ncol)) == 2)) {
        stop("All elements of rotations should be two-column matrices")
      }
      if (length(rotations) != nlayers - 1) {
        stop("Rotations should be one element shorter than the number of layers in geog")
      }
      if (any(unlist(lapply(rotations, max)) > geog$gdat[5])) {
        stop("One of more of the rotation cell IDs in rotations exceed the maximum permissible grid cell ID geog")
      }
    }
  }
  if(!is.null(klink)) {
    if (length(klink) != 1 | !inherits(klink, "numeric")) {
      stop("If not NULL, klink should be an integer")
    }
    if (!klink%%1 == 0) {
      stop("If not NULL, klink should be an integer")
    }
  }
  if(!is.logical(island.check) | length(island.check) != 1) {
    stop("island.check should be a single logical")
  }

  if(is.null(geog$links)) {

    if(island.check) {
      geog <- link_islands(geog, klink = klink)
      add_links <- lapply(1:length(geog$layers), function(x) {
        if(x %in% geog$links$layer) {
          lnk <- geog$links[which(geog$links$layer == x),]
          lnk
        } else {
          NULL
        }
      })

    } else {
      warning("layers were not checked for islands. TARDIS paths fail unexpectedly. Consider running with island.check = TRUE or adding links to the geoglist with link_islands()")
      add_links <- lapply(1:length(geog$layers), function(x) {
        NULL
      })
    }

  } else {

    if (!inherits(geog$links, "sf")) {
      stop("geog$links should be an sf data.frame of linestrings")
    }
    if (!inherits(geog$links, "data.frame")) {
      stop("geog$links should be an sf data.frame of linestrings")
    }
    if (ncol(geog$links) != 5) {
      stop("geog$links should match the format of the output of link_islands()")
    }
    if (!all(colnames(geog$links) == c("srt", "end", "layer", "distance","geometry"))) {
      stop("geog$links should match the format of the output of link_islands()")
    }
    if (any(is.na(geog$links))) {
      stop("One or more columns in geog$links contains NA values")
    }
    if (any(geog$links$layer < 1) | any(geog$links$layer%%1 != 0)) {
      stop("Only positive integers are permitted in geog$links$layer")
    }
    if (nlayers == 1) {
      if (length(unique(geog$links$layer)) > nlayers) {
        stop("geog$links contains links for more layers than those present in geog")
      }
    } else {
      if (any(geog$links$layer > nlayers)) {
        stop("geog$links contains values exceeding the number of layers present in geog")
      }
    }
    if (!all(as.vector(st_geometry_type(geog$links)) %in% c("LINESTRING"))) {
      stop("All geometries in geog$links should be lines")
    }
    if (!all(table(st_coordinates(geog$links)[, 3]) == 2)) {
      stop("Each line in geog$links can only contain 2 coordinates (start and end)")
    }
    tests <- sapply(1:nlayers, function(x) {
      ext1 <- geog$gdat[1:4]
      ext2 <- st_bbox(geog$links[which(geog$links$layer == x),])[c(1, 3, 2, 4)]
      if(all(is.na(as.vector(ext2)))) {ext2 <- ext1}
      all(c(ext1[1] <= ext2[1], ext1[2] >= ext2[2], ext1[3] <= ext2[3], ext1[4] >= ext2[4]))
    })
    if(!all(tests)) {
      stop("The extent of geog$links does not fall fully within the extent of geog")
    }
    add_links <- lapply(1:nlayers, function(x) {

      if(x %in% geog$links$layer) {
        lnk <- geog$links[which(geog$links$layer == x),]

        if(inherits(geog$layers, "SpatRaster")) {
          vals <- extract(boundaries(geog$layers[[x]], directions = 8), vect(lnk))
          vals <- vals[complete.cases(vals),]
          if (!all(vals[,2] == 1)) {
            stop(paste0("In layer ", x, ", one or more line start/end points do not fall on cells at the edges of islands"))
          }
          # if (any(table(vals[complete.cases(vals), 1]) > 2)) {
          #   stop(paste0("In layer ", x, ", one or more lines intersect non-masked areas other than at their start and end points"))
          # }

        } else {
          neigh <- sapply(st_touches(geog$layers[[x]]), length)
          crds <- unique(unlist(st_drop_geometry(lnk[,c("srt", "end")])))
          if(any(neigh[match(crds, as.numeric(rownames(geog$layers[[x]])))] == 6)) {
            stop(paste0("In layer ", x, ", one or more line start/end points do not fall on cells at the edges of islands"))
          }
          #ints <- st_intersects(lnk, geog$layers[[x]])
          #if(any(sapply(ints, length) == 1)) {
          #  stop(paste0("In layer ", x, ", one or more lines start and terminate in the same cell"))
          #}
          # if(any(sapply(ints, length) > 2)) {
          #   stop(paste0("In layer ", x, ", one or more lines intersect non-masked areas other than at their start and end points"))
          # }
        }
        lnk

      } else {
        NULL
      }
    })
  }

  if(inherits(geog$layers, "SpatRaster")) {
    ed <- adjacent(geog$layers[[1]], cells = 1:ncell(geog$layers), directions = 8, pairs = T)
    ed <- ed[ed[, 1] < ed[, 2], ]
    ed <- matrix(c(t(cbind(ed, ed[, 2:1]))), ncol = 2, byrow = T)
    h_dists <- distGeo(xyFromCell(geog$layers[[1]], ed[,1]), xyFromCell(geog$layers[[1]], ed[,2]))
    h_ang <- distGeo(xyFromCell(geog$layers[[1]], ed[,1]), xyFromCell(geog$layers[[1]], ed[,2]))

  } else {

    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    pts <- cell_to_point(grid)
    ed <- st_touches(cell_to_polygon(grid))
    ed <- cbind(rep(grid, sapply(ed, length)), grid[unlist(ed)])
    ed <- matrix(match(ed, grid), ncol = 2)
    h_dists <- distGeo(st_coordinates(pts[ed[,1]]), st_coordinates(pts[ed[,2]]))
    h_ang <- bearing(st_coordinates(pts[ed[,1]]), st_coordinates(pts[ed[,2]]))
  }

  glinked <- list()
  for (i in 1:nlayers) {

    if (verbose) {
      cat(paste0("Linking geog [", i, "/", nlayers, "]\r"))
      if (i == nlayers) {cat("\n")}
    }

    if(inherits(geog$layers, "SpatRaster")) {
      matched <- which(ed[,1] %in% which(!is.na(geog$layers[[i]][])) & ed[,2] %in% which(!is.na(geog$layers[[i]][])))
    } else {
      matched <- which(ed[,1] %in% as.numeric(rownames(geog$layers[[i]])) & ed[,2] %in% as.numeric(rownames(geog$layers[[i]])))
    }

    ed2 <- ed[matched,]
    h_dists2 <- h_dists[matched]
    h_ang2 <- h_ang[matched]
    type <- rep(0, length(h_dists2))
    if(!is.null(add_links[[i]])) {
      lnk <- as.matrix(st_drop_geometry(add_links[[i]][,1:2]))
      ed2 <- rbind(ed2, lnk)
      h_dists2 <- c(h_dists2, add_links[[i]]$distance)
      angs <- bearing(st_coordinates(add_links[[i]])[seq(1, nrow(add_links[[i]]) * 2, 2),1:2],
                      st_coordinates(add_links[[i]])[seq(2, nrow(add_links[[i]]) * 2, 2),1:2])
      h_ang2 <- c(h_ang2, angs)
      type <- c(type, rep(1, length(add_links[[i]]$distance)))
    }
    rownames(ed2) <- NULL

    if(inherits(geog$layers, "SpatRaster")) {
      v_dists <- geog$layers[[i]][][ed2[, 2]] - geog$layers[[i]][][ed2[, 1]]
    } else {
      v_dists <- geog$layers[[i]][[1]][match(ed2[,1], as.numeric(rownames(geog$layers[[i]])))] - geog$layers[[i]][[1]][match(ed2[,2], as.numeric(rownames(geog$layers[[i]])))]
    }

    t_dists <- sqrt(h_dists2^2 + (abs(v_dists)^2))
    t_dists[which(type == 1)] <- h_dists2[which(type == 1)] + v_dists[which(type == 1)]
    edge <- cbind(from = ed2[, 1], to = ed2[, 2], type = type, bearing = h_ang2, hdist = h_dists2, vdist = v_dists, gdist = t_dists)
    glinked[[i]] <- edge[complete.cases(edge), ]
  }

  if (!is.null(times)) {
    for (i in 1:length(rotations)) {
      glinked[[i + 1]][, 1:2] <- glinked[[i + 1]][, 1:2] +
        (i * geog$gdat[5])
      ob <- rotations[[i]]

      angs <- rep(0, nrow(ob))
      to_do <- which(ob[,1] != ob[,2])
      if(length(to_do) > 0) {
        if(is.na(geog$gdat[7])) {
          angs[to_do] <- bearing(xyFromCell(geog$layers[[1]], ob[to_do,1]), xyFromCell(geog$layers[[1]], ob[to_do,2]))
        } else {
          angs[to_do] <- bearing(st_coordinates(cell_to_point(grid[ob[to_do,1]])), st_coordinates(cell_to_point(grid[ob[to_do,2]])))
        }
      }
      ob[, 1] <- ob[, 1] + ((i - 1) * geog$gdat[5])
      ob[, 2] <- ob[, 2] + (i * geog$gdat[5])
      ob[!ob[, 1] %in% glinked[[i]][, 1], 1] <- NA
      ob[!ob[, 2] %in% glinked[[i + 1]], 2] <- NA
      angs <- angs[complete.cases(ob)]
      ob <- ob[complete.cases(ob), , drop = F]
      angs <- angs[!duplicated(ob)]
      ob <- ob[!duplicated(ob), , drop = F]
      if (nrow(ob) == 0) {
        stop(paste0("No links are available from time layer ", i, ". Check rotations and layer masks"))
      }
      if(tlink == 2) {ob <- ob[, 2:1]}
      if(tlink == 3) {
        ob <- matrix(c(t(cbind(ob, ob[, 2:1]))), ncol = 2, byrow = T)
      }
      glinked[[i]] <- rbind(glinked[[i]], cbind(from = ob[, 1], to = ob[, 2], type = 2, bearing = angs, hdist = 0, vdist = 0, gdist = 0))
    }
  }
  if (verbose) {
    cat("Building graph\n")
  }
  glinked <- do.call(rbind, glinked)
  forbidden <- c(glinked[which((glinked[, 1] %% geog$gdat[5]) ==
                                 0), 1], glinked[which((glinked[, 2] %% geog$gdat[5]) ==
                                                         0), 2])
  if (length(forbidden) > 0) {
    glinked <- glinked[!glinked[, 1] %in% forbidden, ]
    glinked <- glinked[!glinked[, 2] %in% forbidden, ]
  }
  src <- as.character(glinked[, 1])
  dst <- as.character(glinked[, 2])
  nodes <- unique(c(src, dst))
  id <- 0:(length(nodes) - 1)

  out <- list(edges = glinked, gdat = geog$gdat, tdat = times, tlink = tlink,
              tgraph = list(data = NULL, dict = data.frame(ref = nodes, id = id),
                            coords = NULL, nbnode = length(nodes),
                            attrib = list(aux = NULL, cap = NULL, alpha = NULL, beta = NULL),
                            src = src, dst = dst))
  class(out) <- "tardis"
  return(out)


}
