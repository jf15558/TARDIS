#' link_mask
#'
#' Check a set of geographic mask rasters to detect islands and calculate
#' optimal bridging points between them. The function has two modes. Mode
#' 'cells' returns the pairs of cell IDs forming each bridge. Cell IDs are
#' given as R's standard cell number order for raster objects. Mode 'cells'
#' is called internally by link_geog(), but it may be desirable to inspect
#' where the bridges have been placed. Mode 'lines' returns a set of lines
#' which can be plotted with the masks to visualize the bridging solutions.
#'
#' @param mask A SpatRaster of geographic masks containing 1 (unmasked) or NA
#' (masked) values
#' @param mode One of 'cells' or 'lines'. If you are using this function
#' directly, then the default 'lines' is probably the desired result.
#' @param alg The linkage algorithm to be used. By default, this is "v", the
#' voronoi linkage method which is superior and generally quite robust. This
#' could still fail though for extremely topologically complex masks, so
#' alternatively, "k" can be used to implement k-nearest neighbour links.
#' This produces similar results to the other method, and is guaranteed to
#' return complete linkage.
#' @param verbose A logical to determine whether function progress should be
#' reported. Useful when dealing with large rasters (high resolution and/or many
#' layers).
#' @return If mode = 'cells', a list of matrices with as many elements as
#' layers in x. Each matrix consists of two columns of cell IDs forming the
#' bridge pairs for its corresponding layer in x. If mode = 'lines', then a
#' list of sf multilinestrings, each of which can be plotted with its
#' corresponding layer in x to visualise the bridges.
#' @import terra sf
#' @importFrom nngeo st_connect
#' @export
#'
#' @details The "v"' bridging algorithm generates Voronoi cells around each
#' detected island. An island is then connected to all other islands whose
#' Voronoi cells are adjacent to its own Voronoi cell. The links themselves
#' are made between the closest cells in each island, producing the
#' minimum planar graph representation for the landscape (Fall et al. 2007:
#' Spatial Graphs: Principles and Applications for Habitat Connectivity). This
#' linkage occurs simultaneously and is driven primarily by functions from
#' the sf package
#'
#' The "k" bridging algorithm works by detecting all islands of cells, counting
#' the number of islands, and assigning each cell a unique island-wise ID.
#' While the number of unique island IDs is greater than one: for each island,
#' its four nearest neighbors are identified and the pairs of cells for those
#' neighbor connections are identified. The cells in all the islands comprising
#' the neighborhood are then assigned the same numeric island ID, progressively
#' reducing the number of island IDs. The pairs of cells recorded at each stage
#' are then returned as the bridging solution. This strategy generally provides
#' more realistic connections between closely positioned sets of islands
#' compared to simply joining to a single nearest neighbour, which would be
#' equivalent to creating a minimum spanning tree for that set. k=4 connections
#' produces superficially similar results to the "v" algorithm, but connections
#' are identified using the st_connect function from nngeo, which assumes
#' Euclidean geometry, while the "v" algorithm uses great circle distances
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#' #data("galapagos")
#' #gal <- unwrap(galapagos)
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' #v <- link_mask(gal_m)
#' #k <- link_mask(gal_m, alg = "k")
#' #plot(gal_m[[1]])
#' #plot(v[[1]], add = T)
#' #plot(k[[1]], add = T, col = 2)

link_mask <- function(mask, mode = "lines", alg = "v", verbose = TRUE) {

  # mask = gal_m
  # mode = "lines"
  # alg = "v"
  # verbose = TRUE

  # check x is correctly supplied
  if(!exists("mask")) {
    stop("Supply mask as SpatRaster")
  }
  if(!inherits(mask, "SpatRaster")) {
    stop("Supply mask as a SpatRaster")
  }
  if(!is.lonlat(mask)) {
    stop("mask should be in geographic (long-lat) projection")
  }
  if(any(!unique(mask[]) %in% c(1, NA))) {
    stop("Mask layers can only contain 1 or NA values")
  }
  if(length(mode) != 1 | inherits(mode, "character")) {
    stop("mode should be one of 'lines' or 'cells'")
  }
  if(!mode %in% c("lines", "cells")) {
    stop("mode should be one of 'lines' or 'cells'")
  }
  if(length(alg) != 1 | inherits(alg, "character")) {
    stop("alg should be one of 'v' or 'k'")
  }
  if(!alg %in% c("v", "k")) {
    stop("alg should be one of 'v' or 'k'")
  }

  bar <- rast(lapply(mask, patches, allowGaps = F))
  res_list <- list()
  linest <- NA
  cls <- cbind(NA, NA)
  for(i in 1:nlyr(bar)) {

    if(verbose) {
      cat(paste0("Resolving mask [", i, "/", nlyr(bar), "]\r"))
      if(i == nlyr(bar)) {cat("\n")}
    }

    if(mode == "lines") {
      res_list[[i]] <- linest
    } else {
      res_list[[i]] <- cls
    }

    if(minmax(bar[[i]])[2] > 1) {

      if(alg == "v") {
        # set up mask and island objects
        foo2 <- boundaries(bar[[i]])
        poly <- st_as_sf(as.polygons(bar[[i]]))
        poly2 <- st_transform(poly, "+proj=cea")

        coords <- tapply(which(foo2[] == 1), extract(bar[[i]], which(foo2[] == 1)), function(x) {xyFromCell(bar[[i]], x)})
        coords <- st_as_sfc(lapply(coords, st_multipoint))

        # generate voronoi polygons for all vertices
        vv <- st_voronoi(st_combine(poly2))
        vv <- st_collection_extract(vv, 'POLYGON')
        vv <- st_crop(vv, st_bbox(poly2))

        # determine which voronoi polygons intersect with input polygons
        ii <- st_intersects(poly2, vv)

        # union/dissolve voronoi polygons that belong to the same inputs
        vt <- lapply(ii, function(x) {st_union(st_combine(vv[x]))})
        vt <- suppressWarnings(suppressMessages(lapply(vt, st_collection_extract)))
        vt <- lapply(vt, st_cast, "POLYGON")
        vt <- lapply(vt, st_union)
        vt <- st_sfc(do.call(rbind, vt))
        vt <- st_intersects(vt, vt)

        # get edge cells in each cluster
        coords <- tapply(which(foo2[] == 1), extract(bar[[i]], which(foo2[] == 1)), function(x) {xyFromCell(bar[[i]], x)})
        coords <- st_as_sfc(lapply(coords, st_multipoint))

        # get nearest neighbours for adjacent voronoi polygons
        res <- suppressMessages(lapply(1:length(vt), function(x) {st_nearest_points(coords[x], coords[vt[[x]]])}))
        res <- Reduce(c, res)
        st_crs(res) <- "+proj=longlat"
        touches <- suppressMessages(st_intersects(res, poly$geometry, sparse = F))
        lin <- res[which(apply(touches, 1, sum) == 2)]

        if(nrow(unique(extract(bar[[i]], st_coordinates(lin)[,-3]))) != minmax(bar[[i]])[2]) {
          stop(paste0("v algorithm failed to link all islands on layer ", i, ". Retry with k"))
        }

      } else {

        crds <- list()
        storeval <- 1
        while(minmax(bar[[1]])[2] > 1) {

          # get edge cells in each cluster
          foo <- boundaries(bar[[i]])
          foo2 <- extract(bar[[i]], which(foo[] == 1))
          coords <- tapply(which(foo[] == 1), foo2, function(x) {xyFromCell(bar[[i]], x)})
          coords <- st_as_sfc(lapply(coords, st_multipoint))

          # get nearest k neighbours between clusters as line links
          if(length(coords) > 2) {
            lin <- suppressMessages(st_connect(coords, coords, k = ifelse(length(coords) > 3, 4, 3), progress = FALSE))
          } else {
            lin <- suppressMessages(st_connect(coords[-1], coords[-2], k = 1, progress = FALSE))
          }

          # coerce to source and destination points, dropping self links
          lin <- matrix(c(t(st_coordinates(lin)[,1:2])), ncol = 4, byrow = T)
          lin <- lin[apply(lin, 1, function(x) {x[1] != x[3] | x[2] != x[4]}), ,drop = F]
          crds[[storeval]] <- lin
          storeval <- storeval + 1

          # convert links to graph and find the new linked clumps
          linked <- cbind(extract(bar[[i]], cellFromXY(bar[[i]], lin[,1:2]))[,1],
                          extract(bar[[i]], cellFromXY(bar[[i]], lin[,3:4]))[,1])
          foo3 <- cbind(1:minmax(bar[[1]])[2], components(graph_from_edgelist(linked))$membership)

          # reclassify raster into the new, aggregated clumps
          bar[[i]] <- classify(bar[[i]], foo3)
        }
        crds <- do.call(rbind, crds)

        # drop any lines which intersect islands other than their start and end ones
        lin <- st_multilinestring(lapply(1:nrow(crds), function(x) {st_linestring(matrix(crds[x,], ncol = 2, byrow = T))}))
        touches <- unlist(do.call(rbind, lapply(lin, function(x) {extract(bar[[i]], vect(x), ID = F, fun = function(x) {length(na.omit(x))})})))
        lin <- st_multilinestring(lin[which(touches == 2)])
      }

      if(mode == "lines") {
        res_list[[i]] <- lin
      } else {

        cls <- matrix(c(t(st_coordinates(lin)[,-3])), ncol = 4, byrow = T)
        cls <- cbind(cellFromXY(bar[[1]], cls[,1:2]), cellFromXY(bar[[1]], cls[,3:4]))
        cls <- t(apply(cls, 1, function(x) {c(min(x), max(x))}))
        cls <- cls[!duplicated(cls),, drop = F]
        cls <- cls[order(cls[,1]),, drop = F]
        res_list[[i]] <- matrix(c(t(cbind(cls, cls[,2:1,drop = F]))), ncol = 2, byrow = T)
      }
    }
  }
  return(res_list)
}
