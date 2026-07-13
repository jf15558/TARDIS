#' point_check
#'
#' Check a set of coordinates to determine whether they fall within an accessible
#' cell within a TARDIS object. Discrepant points are adjusted to their nearest
#' accessible cell based on great circle distance.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param points `matrix` or `data.frame`. A two or three column matrix or
#' data.frame of coordinates. Column ordering is assumed to be longitude
#' (decimal degrees), latitude (decimal degrees) and time (positive, time before
#' present). The time column is only required if `tardis` contains multiple layers.
#' @param max.dist `numeric`. The maximum distance in metres permitted for adjusting
#' points to accessible cells. Points with adjustments above this threshold will
#' be discarded. `NULL` by default, meaning that no points are discarded.
#' @param verbose `logical` Should function progress be reported to the user?
#' @return A `SpatVector` of points, recording which input point they correspond
#' to in case points were discarded (`$feature`), the id of the tardis cell they
#' fall into (`$cell`), the layer to which they belong (`$layer`) and their
#' adjusted distance (`$adj`).
#' @import geosphere terra sf h3jsr
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
#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0), mlink = hlink)
#'
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2),
#'              c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- point_check(htd, org)
#' }

point_check <- function(tardis, points, max.dist = NULL, verbose = TRUE) {


  #org <- hpts[1,]
  #dst <- hpts[2,]
  #tardis = rtd
  #points = chel$biogeography[1:13,3:5]
  #verbose = T
  #max.dist = NULL

  if (!exists("tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if (!inherits(points, "data.frame") & !inherits(points, "matrix")) {
    stop("Supply origin as a data.frame or matrix with two or three columns")
  }
  if (ncol(points) < 2 | ncol(points) > 3) {
    stop("Supply points as a data.frame or matrix with two or three columns")
  }
  if(!is.null(max.dist)) {
    if(!is.atomic(max.dist) | length(max.dist) != 1) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
    if(!is.numeric(max.dist)) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
    if(max.dist <= 0) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
  }

  if (is.null(tardis$tdat)) {
    points <- cbind(points, rep(0.5, nrow(points)))
    tardis$tdat <- c(1, 0)
  }
  if (!is.null(tardis$tdat) & ncol(points) < 3) {
    stop("If tardis contains multiple time layers, then points should contain longitude, latitude and time columns")
  }
  if(is.null(tardis$tdat) & ncol(points) > 2) {
    stop("If tardis consists of one layer, then points should contain longitude and latitude columns only")
  }
  if(!is.na(tardis$gdat[7])) {

    grid <- get_grid(tardis$gdat[1:4], tardis$gdat[7])
    pcell <- suppressMessages(point_to_cell(points[, 1:2, drop = F], tardis$gdat[7]))
    pcell <- match(pcell, grid)
    pt <- unlist(lapply(points[, 3], function(y) {
      sum(y < tardis$tdat[-1]) * tardis$gdat[5]
    }))
    pt[which(points[, 3] > tardis$tdat[1] | points[, 3] < tardis$tdat[length(tardis$tdat)])] <- NA
    ptcell <- pcell + pt
    pmod <- rep(0, length(ptcell))

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
    pcell <- cellFromXY(samprast, points[, 1:2, drop = F])
    pt <- unlist(lapply(points[, 3], function(y) {
      sum(y < tardis$tdat[-1]) * tardis$gdat[5]
    }))
    pt[which(points[, 3] > tardis$tdat[1] | points[, 3] < tardis$tdat[length(tardis$tdat)])] <- NA
    ptcell <- pcell + pt
    pmod <- rep(0, length(ptcell))
  }

  for (i in 1:length(ptcell)) {
    if (verbose) {
      cat(paste0("Checking point [", i, "/", length(ptcell), "]\r"))
      if (i == length(ptcell)) {
        cat("\n")
      }
    }
    if(!is.na(ptcell[i])) {
      if(!ptcell[i] %in% tardis$edges[, 1]) {

        cls <- (pt[i] + 1):(pt[i] + tardis$gdat[5])
        cls <- cls[cls %in% tardis$edges[, 1]] - pt[i]

        if(!is.na(tardis$gdat[7])) {
          cls <- grid[cls]
          cpt <- st_coordinates(cell_to_point(cls))
          dists <- distGeo(points[i,1:2], cpt)
          pmod[i] <- min(dists)
          pcell[i] <- match(suppressMessages(point_to_cell(cpt[which.min(dists),], tardis$gdat[7])), grid)
          ptcell[i] <- pcell[i] + pt[i]

        } else {

          tmprast <- samprast
          tmprast[cls] <- 1
          cpt <- xyFromCell(tmprast, which(tmprast[] == 1))
          dists <- distGeo(points[i,1:2], cpt)
          pmod[i] <- min(dists)
          pcell[i] <- cellFromXY(tmprast, cpt[which.min(dists),,drop = F])
          ptcell[i] <- pcell[i] + pt[i]
        }
      }
    }
  }

  if (all(is.na(ptcell))) {
    stop("No points could be matched to accessible cells in tardis")
  }
  if (any(as.logical(pmod))) {
    warning("Some points were translated to accessible cells in tardis")
  }
  if (any(is.na(pt))) {
    warning("Some points fell outside the temporal bounds of tardis and were discarded")
  }
  out <- cbind.data.frame(point = which(!is.na(ptcell)),
                          bin = pt[!is.na(ptcell)] %/% tardis$gdat[5] + 1,
                          mod = pmod[!is.na(ptcell)])

  if(!is.na(tardis$gdat[7])) {
    geom <- vect(cell_to_point(grid[pcell], tardis$gdat[7]))
    geom$cell <- match(geom$h3_address, grid)

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6], ext = ext(tardis$gdat[1:4]))
    geom <- vect(xyFromCell(samprast, pcell))
    geom$cell <- ptcell
  }
  geom$feature <- out$point
  geom$layer <- out$bin
  geom$adj <- out$mod

  if(!is.null(max.dist)) {
    if(any(geom$adj > max.dist)) {
      warning("Some point adjustments exceeded max.dist and were discarded")
    }
    geom <- geom[which(geom$adj <= max.dist)]
  }
  return(geom[,c(3,2,4,5)])
}
