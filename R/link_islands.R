#' link_islands
#'
#' Detect islands in the layers in a `geoglist` and calculate closest bridging
#' points to maintain connectivity between these isolated regions.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param klink `integer` or `NULL`. The k-nearest neighbours to be linked for
#' each island. Defaults to `NULL`, in which case the Voronoi neighbourhood of
#' each island will determine its k-nearest neighbours automatically (see @details).
#' @param verbose `logical`. Should function progress be to the user? This may
#' be useful when dealing with large rasters (high resolution and/or many layers).
#' @return A `geoglist` as supplied in `geog` with the additional element, `links`.
#' This contains an `sf data.frame` of linestrings representing the island links,
#' recording their start and end cell IDs, the landscape layer to which they belong
#' and their lengths in metres. If links were already generated using another
#' function from TARDIS, they are appended to the existing element `links` and
#' duplicate entries removed. If no landscape layers contained islands, then
#' `geog` is returned without this additional slot.
#' @import terra sf h3jsr
#' @importFrom igraph components
#' @importFrom igraph graph_from_edgelist
#' @export
#'
#' @details
#' Like masking, island linkage is another key feature of TARDIS.
#' Landscapes may contain geographically isolated regions, but practically we
#' may still want to permit movement between them without considering the
#' structure of the intervening space. This function assumes that movement across
#' masked space should occur along the geographically shortest routes between islands.
#'
#' Linkage between islands is based on their Voronoi neighbourhoods. This can
#' be thought of as detecting each island within the line-of-site of every
#' other island, then ranking those line-of-site neighbours by the shortest
#' distances required to reach them. Links which intersect other islands or their
#' Voronoi neighbourhoods are considered invalid and discarded internally. The
#' number of links between neighbours is then determined by `klink`.
#'
#' For klink = 1, the returned links will form a minimum spanning tree, ensuring
#' that every island can be reached by some route from any other island. This is
#' a conservative method of linkage and islands may instead have many neighbours
#' of similar distances. Setting `klink` to higher values will therefore create
#' an increasingly densely connected neighbourhood of islands. `klink` can be set
#' as high as the user likes, but the Voronoid neighbourhood of an island will set
#' the maximum number of links which the function will return.
#'
#' Currently, Voronoi neighbourhoods are calculated using `terra::voronoi()`. While
#' terra generally uses spherical geometry, this function appears to be a rare exception
#' where spherical geometry is not yet supported. As such, Voronoi neighbourhoods
#' on global scales are approximate, although they still appear to return reasonable
#' results.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
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
#' # create a geoglist in raster format and mask the sea
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#'
#' # plot the first layer of the geoglist and add the island links
#' layer_ind = 1
#' plot(rasts$layers[[layer_ind]])
#' plot(rasts$links[[layer_ind]]]$geometry[which(rasts$links$bin == layer_ind)], add = T)
#'
#' plot(hexes$layers[[layer_ind]])
#' plot(hexes$links[[layer_ind]]]$geometry[which(hexes$links$bin == layer_ind)], add = T)
#' }

link_islands <- function(geog, klink = NULL, verbose = T) {
  #
  # geog = rasts
  # klink = 1
  # verbose = T
  #
  if(!exists("geog")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!is.null(klink)) {
    if (length(klink) != 1 | !inherits(klink, "numeric")) {
      stop("If not NULL, klink should be an integer")
    }
    if (!klink%%1 == 0) {
      stop("If not NULL, klink should be an integer")
    }
  }
  if(!is.logical(verbose) | length(verbose) != 1) {
    stop("verbose should be logical")
  }

  if(inherits(geog$layers, "SpatRaster")) {

    islands <- rast(lapply(geog$layers, patches, directions = 8, allowGaps = F))
    bounds <- lapply(1:nlyr(islands), function(x) {
      pt <- as.points(mask(islands[[x]], classify(boundaries(islands[[x]]), cbind(0, NA))))
      aggregate(pt, by = "patches")
    })

  } else {

    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    dat <- lapply(geog$layers, function(z) {
      bar <- st_touches(z)
      bar2 <- do.call(rbind, sapply(1:length(bar), function(x) {rbind(c(x, x), cbind(rep(x, length(bar[[x]])), bar[[x]]))}))
      z$patches <- components(graph_from_edgelist(bar2))$membership
      z <- vect(z[,"patches"])
      z2 <- centroids(z[sapply(bar, length) != 6,])
      list(z, aggregate(z2, by = "patches"))
    })
    islands <- lapply(dat, `[[`, 1)
    bounds <- lapply(dat, `[[`, 2)
  }

  res_list <- list()
  for (i in 1:length(bounds)) {

    if (verbose) {
      cat(paste0("Resolving layers [", i, "/", length(bounds), "]\r"))
      if (i == length(bounds)) {cat("\n")}
    }

    if(!all(bounds[[i]]$patches == 1)) {

      crds <- list()
      iter <- 1

      while(!all(bounds[[i]]$patches == 1)) {

        poly <- bounds[[i]]
        vv <- voronoi(poly, bnd = poly)

        if(FALSE) {
          ######## SPHERICAL METHOD (NOT IMPLEMENTED) ########
          #
          # library(sphereTessellation) # removed from CRAN, manual installation needed, possibly unreliable
          # pts <- as.points(poly)
          #
          # #pt <- icosa::PolToCar(geom(pts)[,3:4]) # convert from lon-lat to Cartesian
          # xlat = geom(pts)[,4] * pi/180
          # xlon = geom(pts)[,3] * pi/180
          # x = cos(xlat) * cos(xlon)
          # y = cos(xlat) * sin(xlon)
          # z = sin(xlat)
          # outvec = c(x, y, z)
          # pt <- matrix(outvec/sqrt(sum(outvec^2)), ncol = 3)
          #
          # #vor <- VoronoiOnSphere(pt, radius = 6371.007)
          # vor <- VoronoiOnSphere(unique(pt)) # spherical voronoi, then get cell coordinates
          # vor <- lapply(vor, `[[`, "cell")
          #
          # #vor <- lapply(vor, function(x) {icosa::CarToPol(t(x))[,1:2]}) # convert from Cartesian to lon-lat
          # vor <- lapply(vor, function(x) {
          #   xyz <- t(x)
          #   x = xyz[,1]
          #   y = xyz[,2]
          #   z = xyz[,3]
          #   latitude = 180 * asin(z)/pi
          #   longitude = 180 * atan2(y, x)/pi
          #   cbind(longitude, latitude)
          # })
          #
          # vor <- lapply(vor, function(x) {x[c(1:nrow(x), 1),]})
          # vor <- lapply(vor, function(x) {st_polygon(list(x))})
          # vor <- st_sfc(vor, crs = "+proj=lonlat")
          # vor2 <- data.frame(patches = pts$patches)
          # st_geometry(vor2) <- vor
          # vor2 <- st_wrap_dateline(vor2, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
          # vor3 <- vect(vor2)
          # vor3 <- makeValid(vor3)
          # vor3 <- erase(vor3, poly)
          # vv <- vor3[which(tapply(geom(vor3)[,2], geom(vor3)[,1], function(x) {length(unique(x))}) == 1)]
          #
          ##################################
        }

        vv <- aggregate(vv, by = "patches")
        ii <- adjacent(vv)

        # dropping duplicates helps resolve cases where start-end cell for a link pair are not the same
        ii <- unique(t(apply(ii, 1, sort)))
        dl <- do.call(rbind, apply(ii, 1, function(x) {
          ln <- nearest(poly[x[1]], poly[x[2]], centroids = F, lines = T)
          ln[which.min(ln$distance)]
        }))

        if(inherits(geog$layers, "SpatRaster")) {
          nr <- extract(islands[[i]], dl)
          nr <- subset(nr, complete.cases(nr))
          to_keep <- which(tapply(nr$patches, nr$ID, length) == 2)
          #to_keep <- which(tapply(nr$patches, nr$ID, function(x) {length(unique(x))}) == 2)
          dl <- dl[to_keep]
          ii <- ii[to_keep,]
          cl <- cellFromXY(geog$layers[[1]], geom(dl)[,3:4])
          cls <- as.data.frame(matrix(cl, ncol = 2, byrow = T))

        } else {

          nr <- apply(relate(dl, islands[[i]], "intersects"), 1, which, simplify = F)
          to_keep <- which(sapply(nr, length) == 2)
          dl <- dl[to_keep]
          ii <- ii[to_keep,]
          cl <- suppressMessages(point_to_cell(geom(dl)[,3:4], res = geog$gdat[7]))
          cl <- match(cl, grid)
          cls <- as.data.frame(matrix(cl, ncol = 2, byrow = T))
        }

        dists <- perim(dl)
        links <- cbind(ii, dists)

        # duplicate
        links <- rbind.data.frame(links, links[,c(2, 1, 3)])
        dl <- rbind(dl, dl)
        cls <- rbind.data.frame(cls, cls[,2:1])
        colnames(links) <- c("from", "to", "dists")

        dl <- dl[order(links$from, links$dists)]
        cls <- cls[order(links$from, links$dists),]
        links <- links[order(links$from, links$dists),]

        newlink <- tapply(links$dists, links$from, function(x) {
          getlink <- ifelse(is.null(klink), length(x), klink)
          1:ifelse(getlink < length(x), getlink, length(x))
        })
        newlink <- unlist(mapply(x = newlink, y = which(!duplicated(links$from)) - 1, function(x, y) {x + y}, SIMPLIFY = F))

        colnames(cls) <- c("srt", "end")
        cls$layer <- i
        cls$distance <- links$dists
        st_geometry(cls) <- st_as_sf(dl)$geometry
        finals <- cls[newlink,]
        crds[[iter]] <- finals

        newlink <- rbind(as.matrix(links[newlink, 1:2]), cbind(1:max(bounds[[i]]$patches),
                                                               1:max(bounds[[i]]$patches)))
        newid <-  components(graph_from_edgelist(newlink))$membership
        bounds[[i]]$patches <- newid
        bounds[[i]] <- aggregate(bounds[[i]], by = "patches")

        iter <- iter + 1
        if (all(bounds[[i]]$patches == 1)) {
          iter <- FALSE
        }
      }
      res_list[[i]] <- do.call(rbind, crds)
    }
  }
  if(all(sapply(res_list, is.null))) {
    message("No islands found in any layers, no links will be returned")
    return(NULL)
  } else {
    geog$links <- unique(rbind(geog$links, do.call(rbind, res_list)))
    geog$links <- st_wrap_dateline(lnks, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    return(geog)
  }
}
