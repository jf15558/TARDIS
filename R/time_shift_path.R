#' time_shift_path
#'
#' For a given input cell and destination layer, get all its intervening positions
#' from the inter-layer connections in a TARDIS graph. If the cell location is
#' lost in an intervening layer, the path will terminate at its last available
#' location, rather than in the target layer. This function works for any TARDIS
#' graph, but only produces meaningful rotation paths when locations are not
#' homologous between layers (i.e., where locations shift due to plate tectonic
#' rotation or some other process).
#'
#' @param tardis `tardis` A tardis graph as produced by build_tardis().
#' @param points `numeric` or `SpatVector`. Either a numeric vector of cell IDs
#' a `SpatVector` of points with a cell ID variable named 'cell'.
#' @param layer `numeric`. The layer number to shift points to. If a single
#' number, then all points (regardless of their starting layer) will be shifted
#' to that layer. Alternatively a vector of numbers with as many elements as
#' points to permit different destination layers for each point.
#' @param time `numeric`. Instead of a layer number, numeric ages can be
#' supplied. These are resolved to layer numbers internally, then treated like
#' the `layer` argument.
#' @param as.lines `logical`. Should the point-to-point shifts in each layer be
#' returned as line objects? These may be convenient for plotting purposes.
#' `FALSE` by default.
#' @return A SpatVector of points or lines, recording the input point to which
#' they belong (`$feature`), the time layer to which each point or line segment
#' belongs (`$layer`), and the starting and ending cell IDs for each row-wise
#' rotation (`$from`, `$to`).
#' @import sf terra h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' #' library(terra)
#' library(TARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' hexes <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2),
#'              c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- point_check(htd, org)
#'
#' hpts2 <- time_shift_path(htd, hpts, time = 0)
#' }

time_shift_path <- function(tardis, points, layer, time = NULL, as.lines = F) {

  if(!exists("tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if(length(tardis$tdat < 3)) {
    stop("tardis only contains a single layer. No shifts available")
  }

  if(!exists("points")) {
    stop("Please supply points as a vector of TARDIS cell IDs")
  }
  if(inherits(points, "SpatVector")) {
    if(geomtype(points) != "points") {
      stop("points must be a SpatVector of points geometries only")
    }
    if(!"cell" %in% names(points)) {
      stop("If supplying points as a SpatVector, this must contain a column of numeric cell IDs named 'cell'")
    }
    if(!is.numeric(points$cell)) {
      stop("points$cell must be numeric")
    }
    if(any(is.na(points$cell))) {
      stop("points$cell cannot contain missing values")
    }
    pts <- points$cell
  } else {
    pts <- points
  }

  if(!is.atomic(pts) | !is.character(pts)) {
    stop("Please supply points as a vector of TARDIS cell IDs")
  }
  if(any(pts %% tardis$gdat[5]) > tardis$gdat[5]) {
    stop("Some points exceed the maximum permitted cell layer number. Check if these points correspond to the supplied tardis object")
  }
  if(any(pts %/% tardis$gdat[5] > (length(tardis$tdat) - 1))) {
    stop("Some points exceed the number of layers in tardis. Check if these points correspond to the supplied tardis object")
  }
  tests <- as.character(pts) %in% tardis$tgraph$dict$ref
  if(!all(tests)) {
    stop(paste0(length(tests) - sum(tests), " points do not correspond to tardis cell IDs. Use point_check() to validate starting cells first"))
  }

  if(!exists("layer") & is.null(time)) {
    stop("One of layer or time must be supplied")
  }
  if(!is.null(time)) {
    if(!is.atomic(time)) {
      stop("time should be a singe positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(!is.numeric(time)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(any(time < 0)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(length(time) != 1 | length(time) != length(pts)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(any(time < min(tardis$tdat) | time > max(tardis$tdat))) {
      stop("One or more elements in time does not fall within the temporal range of tardis")
    }
    layer <- sum(which(time > tardis$tdat))
  }

  if(!is.atomic(layer)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(!is.numeric(layer)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(any(layer %% 1 != 0 | layer < 0)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(any(layer > (length(tardis$tdat) - 1))) {
    stop("One or more elements in layer exceeds the number of layers in tardis")
  }

  pts <- cbind(1:length(pts), pts, (pts %/% tardis$gdat[5]) - 1, layer)
  rot <- tardis$edges[which(tardis$edges[,3] == 2),1:2]
  if(!is.na(tardis$gdat[7])) {
    grd <- get_grid(tardis$gdat[1:4], tardis$gdat[7])
  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
  }

  pth <- list()
  for(i in 1:nrow(pts)) {

    pth1 <- pt <- pts[i]
    bns <- ((pts[i] %/% tardis$gdat[5]) - 1):layer[i]
    for(k in 1:(length(bns) - 1)) {

      # if not the same bin
      if(bns[1] != bns[length(bns)]) {

        # forwards in time
        if(bns[1] < bns[2]) {
          rt <- rot[which(((rot[,1] %/% tardis$gdat[5]) + 1) == k),]

          # backwards in time
        } else {
          rt <- rot[which(((rot[,2] %/% tardis$gdat[5]) + 1) == k),2:1]
        }
        pt <- rt[match(pt, rot[,2]),1]
        if(!is.na(pt)) {
          pth1 <- c(pth1, pt)
        }
      }
    }
    if(length(pth1) == 1) {
      pth1 <- c(pth1, pth1)
    }
    if(!is.na(tardis$gdat[7])) {
      crd <- st_coordinates(cell_to_point(grd[pth1]))[,3:4, drop = F]
    } else {
      crd <- xyFromCell(samprast, pth1)
    }
    if(as.lines) {
      geoms <- vect(cbind(crd[-nrow(crd),,drop = F], crd[-1,,drop=F]))
    } else {
      geoms <- vect(crd[-1])
    }
    geoms$feature <- i
    geoms$layer <- bns
    geoms$from <- pth1[-length(pth1)]
    geoms$to <- pth1[-1]
    ptl[[i]] <- geoms
  }
  pth <- do.call(rbind, pth)
  return(pth)
}
