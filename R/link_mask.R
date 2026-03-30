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
#' @param glink The linkage mode for island cell clusters - 4 (Rook's case,
#' only orthogonal neighbours considered as touching) or 8 (Queen's case,
#' orthogonal and diagonal neighbours considered as touching).
#' @param klink The k-nearest neighbours to be linked for each island. If NULL
#' (default), the Voronoi neighbourhood of each island will be linked, i.e,
#' each neighbour which can be accessed in a straight line without passing
#' through any others first. Otherwise a numeric for the desired number of
#' links per island. For klink = 1, the returned links will form a minimum
#' spanning tree. klink can be set as high as the user likes, but any links not
#' part of the Voronoi linkage for a given polygon will be discarded.
#' @param verbose A logical to determine whether function progress should be
#' reported. Useful when dealing with large rasters (high resolution and/or many
#' layers).
#' @return If mode = 'cells', a list of matrices with as many elements as
#' layers in x. Each matrix consists of two columns of cell IDs forming the
#' bridge pairs for its corresponding layer in x. If mode = 'lines', then a
#' list of sf multilinestrings, each of which can be plotted with its
#' corresponding layer in x to visualise the bridges.
#' @import terra sf
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#' #gal <- galapagos()
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' #v <- link_mask(gal_m)
#' #k <- link_mask(gal_m, alg = "k")
#' #plot(gal_m[[1]])
#' #plot(v[[1]], add = T)
#' #plot(k[[1]], add = T, col = 2)

link_mask <- function (mask, glink = 8, klink = NULL, verbose = TRUE) {

  mask = msk
  glink = 8
  klink = 1
  verbose = TRUE

  if (!exists("mask")) {
    stop("Supply mask as SpatRaster")
  }
  if (!inherits(mask, "SpatRaster")) {
    stop("Supply mask as a SpatRaster")
  }
  if (!is.lonlat(mask)) {
    stop("mask should be in geographic (lon-lat) projection")
  }
  if (any(!unique(mask[]) %in% c(1, NA))) {
    stop("Mask layers can only contain 1 or NA values")
  }
  if (!is.numeric(glink) | length(glink) != 1) {
    stop("glink should be one of 4 (Rook's case) or 8 (Queen's case)")
  }
  if (!glink %in% c(4, 8)) {
    stop("glink should be one of 4 or 8")
  }
  if (!is.null(klink)) {
    if (length(klink) != 1 | !inherits(klink, "numeric")) {
      stop("If not NULL, klink should be an integer")
    }
    if (!klink%%1 == 0) {
      stop("If not NULL, klink should be an integer")
    }
  }
  bar <- rast(lapply(mask, patches, directions = glink, allowGaps = F))
  res_list <- list()
  for (i in 1:nlyr(bar)) {
    if (verbose) {
      cat(paste0("Resolving mask [", i, "/", nlyr(bar),
                 "]\r"))
      if (i == nlyr(bar)) {
        cat("\n")
      }
    }
    if (minmax(bar[[i]])[2] > 1) {
      crds <- list()
      iter <- 1
      while (as.logical(iter)) {

        # get the voronoi neighbourhood of each island
        poly <- as.polygons(bar[[i]])
        vv <- voronoi(poly)
        vv <- crop(vv, bar[[i]])
        vv <- aggregate(vv, by = "patches")

        # get the distances between islands in adjacent patches
        ii <- adjacent(vv)
        nr <- do.call(rbind, apply(ii, 1, function(x) {
          nearest(poly[x[1]], poly[x[2]], centroids = F, lines = T)
        }))
        vals <- extract(mask[[i]], nr)
        dists <- perim(nr)
        nr <- geom(nr)
        nr <- sapply(seq(2, nrow(nr), 2), function(x) {nr[c(x-1, x),3:4]}, simplify = F)
        links <- cbind.data.frame(ii, dists)
        nr <- nr[order(links$from, links$dists)]
        links <- links[order(links$from, links$dists),]

        # adjust coordinates of lines from cell corners to cell centres
        bnd <- mask(bar[[i]], classify(boundaries(bar[[i]]), cbind(0, NA)))
        bnd[which(!is.na(bnd[]))] <- which(!is.na(bnd[]))
        to_cell <- as.points(as.polygons(bnd))
        to_cell <- cbind(to_cell$patches, geom(to_cell)[,3:4])
        nr <- lapply(nr, function(x) {
          st_linestring(xyFromCell(bar[[i]], to_cell[c(which(to_cell[,2] == x[1,1] & to_cell[,3] == x[1,2])[1],
                                                       which(to_cell[,2] == x[2,1] & to_cell[,3] == x[2,2])[1]),1]))
        })

        # nr <- lapply(nr, function(x) {
        #   cl <- xyFromCell(bar[[i]], to_cell[c(which(to_cell[,2] == x[1,1] & to_cell[,3] == x[1,2])[1],
        #                                        which(to_cell[,2] == x[2,1] & to_cell[,3] == x[2,2])[1]),1])
        #   rbind(cl[1,], gcIntermediate(cl[1,], cl[2,]), cl[2,])
        # })

        # reject links which cut through other islands
        checks <- extract(mask[[i]], vect(nr))
        checks <- which(table(checks[complete.cases(checks),1]) > 2)
        if(length(checks) != 0) {
          nr <- nr[-checks]
          links <- links[-checks,]
        }

        # take up to klink links for each island patch
        newlink <- tapply(links$dists, links$from, function(x) {
          1:ifelse(klink < length(x), klink, length(x))
        })
        newlink <- unlist(mapply(x = newlink, y = as.list((which(!duplicated(links$from)) - 1)), function(x, y) {x + y}))
        crds[[iter]] <- nr[newlink]
        newlink <- as.matrix(links[newlink,1:2])

        # reclassify island membership by links
        newid <- cbind(1:minmax(bar[[i]])[2], components(graph_from_edgelist(newlink))$membership)
        bar[[i]] <- classify(bar[[i]], newid)

        # stop loop if all islands have the same membership
        iter <- iter + 1
        if(minmax(bar[[i]])[2] == 1) {
          iter <- FALSE
        }
      }
      lin <- do.call(c, crds)
      cls <- matrix(cellFromXY(bar[[i]], st_coordinates(lin)[,1:2]), ncol = 2, byrow = 2)
      lin <- st_sf(data.frame(srt = cls[,1], end = cls[,2], bin = rep(i, length(lin))), distance = as.vector(st_length(lin)), geometry = lin)
      res_list[[i]] <- lin
    }
  }
  res_list <- do.call(rbind.data.frame, res_list)
  return(res_list)
}
