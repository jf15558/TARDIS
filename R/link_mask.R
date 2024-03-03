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
#' @importFrom nngeo st_connect
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

link_mask <- function(mask, glink = 8, klink = NULL, verbose = TRUE) {

  # mask = masks
  # glink = 8
  # klink = 1
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
  if(!is.numeric(glink) | length(glink) != 1) {
    stop("glink should be one of 4 (Rook's case) or 8 (Queen's case)")
  }
  if(!glink %in% c(4, 8)) {
    stop("glink should be one of 4 or 8")
  }
  if(!is.null(klink)) {
    if(length(klink) != 1 | !inherits(klink, "numeric")) {
      stop("If not NULL, klink should be an integer")
    }
    if(!klink %% 1 == 0) {
      stop("If not NULL, klink should be an integer")
    }
  }

  bar <- rast(lapply(mask, patches, directions = glink, allowGaps = F))
  res_list <- list()
  for(i in 1:nlyr(bar)) {

    if(verbose) {
      cat(paste0("Resolving mask [", i, "/", nlyr(bar), "]\r"))
      if(i == nlyr(bar)) {cat("\n")}
    }

    if(minmax(bar[[i]])[2] > 1) {

      crds <- list()
      iter <- 1
      while(as.logical(iter)) {

        # set up mask and island objects
        bnd <- boundaries(bar[[i]], directions = glink)
        poly <- st_as_sf(as.polygons(bar[[i]]))
        poly2 <- st_transform(poly, "+proj=cea")


        # pts <- lapply(poly$geometry, function(x) {st_coordinates(x)[,1:2]})
        # pts <- cbind(do.call(rbind, pts), unlist(lapply(1:length(pts), function(x) {rep(x, nrow(pts[[x]]))})))
        # vert <- pts[!duplicated(pts),]
        # xlat = as.double(vert[,2]) * pi/180
        # xlon = as.double(vert[,1]) * pi/180
        # x = cos(xlat) * cos(xlon)
        # y = cos(xlat) * sin(xlon)
        # z = sin(xlat)
        # outvec = cbind(x, y, z)
        #
        # vor <- VoronoiOnSphere(outvec)
        # vor <- lapply(vor, function(x) {
        #   xyz <- t(x$cell)
        #   lat <- 180 * asin(xyz[,3]) / pi
        #   lon <- 180 * atan2(xyz[,2], xyz[,1]) / pi
        #   lonlat <- cbind(lon, lat)
        #   st_polygon(list(lonlat[c(1:nrow(lonlat), 1),]))
        # })
        # vor <- (st_sfc(vor, crs = "+proj=lonlat"))


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
        coords <- tapply(which(bnd[] == 1), terra::extract(bar[[i]], which(bnd[] == 1)), function(x) {xyFromCell(bar[[i]], x)})
        coords <- st_as_sfc(lapply(coords, st_multipoint))

        # get nearest neighbour lines for adjacent voronoi polygons
        res <- suppressMessages(lapply(1:length(vt), function(x) {st_nearest_points(coords[x], coords[vt[[x]]])}))
        res <- lapply(res, function(x) {st_crs(x) <- "+proj=lonlat"; x})

        # get the lengths of each line, drop zero-length connections, return desired connection number (or the max if kcon exceeds the available connections)
        len <- lapply(res, function(x) {as.vector(st_length(x))})
        res <- mapply(x = res, y = len, function(x, y) {
          x <- x[order(y)]
          y <- y[order(y)]
          x <- x[which(y > 0)]
          y <- y[which(y > 0)]
          if(is.null(klink)) {kn <- length(y)} else {
            if(klink > length(y)) {kn <- length(y)} else {kn <- klink}
          }
          x[1:kn]
        }, SIMPLIFY = F)
        res <- Reduce(c, res)

        # remove any lines that intersect the map polygons
        touches <- suppressMessages(st_intersects(res, poly$geometry, sparse = F))
        lin <- res[which(apply(touches, 1, sum) == 2)]

        # coerce to source and destination points, dropping self links
        lin <- matrix(c(t(st_coordinates(lin)[,1:2])), ncol = 4, byrow = T)
        lin <- lin[apply(lin, 1, function(x) {x[1] != x[3] | x[2] != x[4]}), ,drop = F]

        # ensuring that duplicate connections (or reciprocal) are removed
        cl <- cbind(cellFromXY(mask, lin[,1:2]), cellFromXY(mask, lin[,3:4]))
        cl <- t(apply(cl, 1, function(x) {c(min(x), max(x))}))
        crds[[iter]] <- lin[!duplicated(cl),]
        iter <- iter + 1

        # convert links to graph and find the new linked clumps
        linked <- cbind(terra::extract(bar[[i]], cellFromXY(bar[[i]], lin[,1:2, drop = F]))[,1],
                        terra::extract(bar[[i]], cellFromXY(bar[[i]], lin[,3:4, drop = F]))[,1])
        foo3 <- cbind(1:minmax(bar[[i]])[2], components(graph_from_edgelist(linked))$membership)

        # reclassify raster into the new, aggregated clumps
        bar[[i]] <- classify(bar[[i]], foo3)
        if(minmax(bar[[i]])[2] == 1) {iter <- FALSE}
      }
      crds <- do.call(rbind, crds)

      # drop any lines which intersect islands other than their start and end ones
      #lin <- st_multilinestring(lapply(1:nrow(crds), function(x) {st_linestring(matrix(crds[x,], ncol = 2, byrow = T))}))
      #touches <- unlist(do.call(rbind, lapply(lin, function(x) {terra::extract(bar[[i]], vect(x), ID = F, fun = function(x) {length(na.omit(x))})})))
      #lin <- st_multilinestring(lin[which(touches == 2)])
      #lin <- st_sfc(lin, crs = "+proj=lonlat")

      # create lines and store
      lin <- st_sfc(lapply(1:nrow(crds), function(x) {st_linestring(matrix(crds[x,], ncol = 2, byrow = T))}), crs = "+proj=lonlat")
      lin <- st_sf(data.frame(srt = cellFromXY(bar[[1]], crds[,1:2]),
                              end = cellFromXY(bar[[1]], crds[,3:4]),
                              bin = rep(i, length(lin))),
                              distance = as.vector(st_length(lin)), geometry = lin)
      res_list[[i]] <- lin
    }
  }
  res_list <- do.call(rbind.data.frame, res_list)
  return(res_list)
}
