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
#' gal <- galapagos()
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

build_tardis <- function(geog, times = NULL, tlink = 1, rotations = NULL, verbose = TRUE) {

   # geog = hexes
   # times = c(seq(2.25, 0, -0.5), 0)
   # times = c(117, 114, 112)
   # tlink = 1
   # rotations = rots
   # verbose = TRUE

  if(inherits(geog$layers, "SpatRaster")) {
    nlayers <- nlyr(geog$layers)
  } else {
    nlayers <- length(geog$layers)
  }
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

  if(!is.null(geog$links)) {

    if (!inherits(geog$links, "SpatVector")) {
      stop("geog$links should be an SpatVector of lines")
    }
    if (dim(geog$links)[2] != 4) {
      stop("geog$links should contain four value fields: srt, end, layer, distance")
    }
    if (!all(names(geog$links) == c("srt", "end", "layer", "distance"))) {
      stop("geog$links should contain four value fields: srt, end, layer, distance")
    }
    if (any(is.na(geog$links))) {
      stop("One or more values in geog$links is NA")
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
    if (!geomtype(geog$links) == "lines") {
      stop("All geometries in geog$links should be lines")
    }
    #if (!all(table(st_coordinates(geog$links)[, 3]) == 2)) {
    #  stop("Each line in geog$links can only contain 2 coordinates (start and end)")
    #}
    tests <- sapply(1:nlayers, function(x) {
      ext1 <- geog$gdat[1:4]
      ext2 <- ext(geog$links[which(geog$links$layer == x),])
      if(all(is.na(as.vector(ext2)))) {ext2 <- ext1}
      all(c(ext1[1] <= ext2[1], ext1[2] >= ext2[2], ext1[3] <= ext2[3], ext1[4] >= ext2[4]))
    })
    if(!all(tests)) {
      stop("The extent of geog$links does not fall fully within the extent of geog")
    }
    add_links <- lapply(1:nlayers, function(x) {

      ## ideally want some checks here for invalid links, but these were always faulty
      if(x %in% geog$links$layer) {
        geog$links[which(geog$links$layer == x)]
      } else {
        NULL
      }
    })
  } else {
    add_links <- lapply(1:nlayers, function(x) {NULL})
  }

  if(inherits(geog$layers, "SpatRaster")) {

    ed <- adjacent(geog$layers[[1]], cells = 1:ncell(geog$layers), directions = 8, pairs = T)
    ed <- ed[ed[, 1] < ed[, 2], ]
    ed <- matrix(c(t(cbind(ed, ed[, 2:1]))), ncol = 2, byrow = T)
    h_dists <- distGeo(xyFromCell(geog$layers[[1]], ed[,1]), xyFromCell(geog$layers[[1]], ed[,2]))
    h_ang <- distGeo(xyFromCell(geog$layers[[1]], ed[,1]), xyFromCell(geog$layers[[1]], ed[,2]))

  } else {

    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    pts <- vect(cell_to_point(grid))
    cls <- vect(cell_to_polygon(grid))
    ed <- relate(cls, cls, "intersects", pairs = T)
    ed <- ed[which(ed[,1] != ed[,2]),]
    h_dists <- distGeo(crds(pts[ed[,1]]), crds(pts[ed[,2]]))
    h_ang <- bearing(crds(pts[ed[,1]]), crds(pts[ed[,2]]))
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
      matched <- which(ed[,1] %in% geog$layers[[i]]$id & ed[,2] %in% geog$layers[[i]]$id)
    }

    ed2 <- ed[matched,]
    h_dists2 <- h_dists[matched]
    h_ang2 <- h_ang[matched]
    type <- rep(0, length(h_dists2))
    if(!is.null(add_links[[i]])) {
      lnk <- as.matrix(values(add_links[[i]])[,1:2])
      ed2 <- rbind(ed2, lnk)
      h_dists2 <- c(h_dists2, add_links[[i]]$distance)
      if(inherits(geog$layers, "SpatRaster")) {
        crd1 <- xyFromCell(geog$layers[[1]], lnk[,1])
        crd2 <- xyFromCell(geog$layers[[1]], lnk[,1])
      } else {
        crd1 <- st_coordinates(cell_to_point(grid[lnk[,1]], geog$gdat[7]))
        crd2 <- st_coordinates(cell_to_point(grid[lnk[,2]], geog$gdat[7]))
      }
      angs <- bearing(crd1, crd2)
      h_ang2 <- c(h_ang2, angs)
      type <- c(type, rep(1, length(add_links[[i]]$distance)))
    }
    rownames(ed2) <- NULL

    if(inherits(geog$layers, "SpatRaster")) {
      v_dists <- geog$layers[[i]][][ed2[, 2]] - geog$layers[[i]][][ed2[, 1]]
    } else {
      v_dists <- geog$layers[[i]][[1]][match(ed2[,1], geog$layers[[i]]$id),] - geog$layers[[i]][[1]][match(ed2[,2], geog$layers[[i]]$id),]
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
