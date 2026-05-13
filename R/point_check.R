#' point_check
#'
#' Check a set of space-time coordinates to determine whether they fall within
#' an accessible cell within a TARDIS object and adjust discrepant points to
#' their nearest accessible cell using great circle distances.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param points `matrix`. A two or three column matrix giving the spatiotemporal coordinates.
#' Column ordering is assumed to be longitude (decimal degrees), latitude (decimal
#' degrees) and time (positive, time before present).
#' @param verbose `logical` Should function progress be reported to the user?
#' @return An `sf data.frame` containing valid points for the input points,
#' recording the layer to which they belong and their adjusted distance (if
#' shifted to accessible cells).
#' @import geosphere terra sf h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' gal <- TARDIS::galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' hexes <- rast_to_geoglist(gal, gal_m, as.hex = T, hex = 7)
#' rasts <- rast_to_geoglist(gal, gal_m)

#' hlink <- link_islands(hexes)
#' rlink <- link_islands(rasts)

#' htd <- build_tardis(hexes, times = c(seq(2.25, 0, -0.5), 0), mlink = hlink)
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0), mlink = rlink)

#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2))
#' dst <- rbind(c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- stp(htd, rbind(org, dst))
#' rpts <- stp(rtd, rbind(org, dst))
#' }

point_check <- function(tardis, points, verbose = TRUE) {


  # org <- hpts[1,]
  # dst <- hpts[2,]
  # tardis = rtd
  # points = rbind(org, dst)
  # verbose = T

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
    pcell <- point_to_cell(points[, 1:2, drop = F], tardis$gdat[7])
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
          pcell[i] <- match(point_to_cell(cpt[which.min(dists),], tardis$gdat[7]), grid)
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
    geom <- cell_to_point(grid[pcell], tardis$gdat[7])

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6], ext = ext(tardis$gdat[1:4]))
    geom <- xyFromCell(samprast, pcell)
    geom <- st_sfc(apply(geom, 1, st_point, simplify = F), crs = "EPSG:4326")
  }
  st_geometry(out) <- st_geometry(geom)
  return(out)
}
