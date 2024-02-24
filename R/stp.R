#' stp
#'
#' Check a set of space-time coordinates to determine whether they fall within
#' an accessible cell within a TARDIS object. The function first checks whether
#' all points fall within the spatial and temporal extents of the object. From
#' there, it checks whether the points fall within an accessible cell in their
#' respective time layers and will adjust discrepant points to their nearest
#' accessible cell using great circle distances. Adjustment is performed using
#' a rip of nearestLand() from seegSDM, up to a maximum of 10000 km, although
#' it is expected that any adjustments will be nowhere near this large.
#'
#' @param tardis A object of class 'tardis' from create_tardis.
#' @param points A two or three column matrix giving the spatiotemporal coordinates.
#' Column ordering is assumed to be longitude (decimal degrees), latitude (decimal
#' degrees) and time (positive, years before present).
#' @param verbose A logical indicating whether function progress should be
#' reported to the user
#' @return A three column data.frame containing the final point longitudes and
#' latitudes, with the third column recording the point adjustment distance in
#' metres. For unadjusted points, this will be zero. For points which could
#' not be successfully adjusted, this will be NA.
#' @importFrom geosphere distGeo
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#'
#' #gal <- galapagos()
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1),
#' #                                    ncol = 3, byrow = T), right = F)
#' #gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
#'
#' #org <- rbind(c(-89, -1.05, 2), c(-89.5, -0.7, 2))
#' #dst <- rbind(c(-91.2, -1, 0), c(-91.6, -0.4, 0))
#' #pts <- stp(gt, rbind(org, dst))

stp <- function(tardis, points, verbose = TRUE) {

  # x = test2
  # points = pairs

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(points, "data.frame") & !inherits(points, "matrix")) {
    stop("Supply origin as a data.frame or matrix with two or three columns")
  }
  if(ncol(points) < 2 | ncol(points) > 3) {
    stop("Supply points as a data.frame or matrix with two or three columns")
  }
  if(is.null(tardis$tdat)) {
    points <- cbind(points, rep(0.5, nrow(points)))
    tardis$tdat <- c(1, 0)
  }
  if(!is.null(tardis$tdat) & ncol(points) < 3) {
    stop("If tardis contains multiple time layers, then points should contain longitude, latitude and time columns")
  }

  samprast <- rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8]))
  pcell <- cellFromXY(samprast, points[,1:2])
  pt <- unlist(lapply(points[,3], function(y) {sum(y < tardis$tdat[-1]) * prod(tardis$gdat[1:2])}))
  pt[which(points[,3] > tardis$tdat[1] | points[,3] < tardis$tdat[length(tardis$tdat)])] <- NA
  ptcell <- pcell + pt
  pmod <- rep(0, length(ptcell))

  # static rip of seegSDM nearestLand()
  nearest_land <- function(points, raster, max_distance) {

    nearest <- function(lis, raster) {

      neighbours <- matrix(lis[[1]], ncol = 2)
      point <- lis[[2]]
      land <- !is.na(neighbours[, 2])
      if (!any(land)) {return(c(NA, NA))} else {
        coords <- xyFromCell(raster, neighbours[land, 1])
        if (nrow(coords) == 1) {
          return(coords[1, ])
        }
        dists <- sqrt((coords[,1] - point[1]) ^ 2 + (coords[,2] - point[2]) ^ 2)
        return(coords[which.min(dists),])
      }
    }
    neighbour_list <- extract(raster, points, buffer = max_distance, cellnumbers = TRUE)
    neighbour_list <- lapply(1:nrow(points), function(y) {
      list(neighbours = neighbour_list[[y]], point = as.numeric(points[y,]))
    })
    return(t(sapply(neighbour_list, nearest, raster)))
  }

  for(i in 1:length(ptcell)) {

    if(verbose) {
      cat(paste0("Checking point [", i, "/", length(ptcell), "]\r"))
      if(i == length(ptcell)) {cat("\n")}
    }

    if(!is.na(ptcell[i])) {
      if(!ptcell[i] %in% tardis$edges[,1]) {

        # get cells in age range
        cls <- (pt[i] + 1):(pt[i] + prod(tardis$gdat[1:2]))
        # drop inaccessible cells
        cls <- cls[cls %in% tardis$edges[,1]] - pt[i]

        # locate nearest accessible cell
        tmprast <- samprast
        tmprast[cls] <- 1
        sub <- nearest_land(xyFromCell(tmprast, pcell[i]), tmprast, max_distance = 1e07)
        pmod[i] <- ptcell[i] <- NA
        if(!is.na(sub[1,1])) {
          pmod[i] <- round(distGeo(xyFromCell(tmprast, pcell[i]), sub), 2)
          pcell[i] <- cellFromXY(samprast, sub)
          ptcell[i] <- pcell[i] + pt[i]
        }
      }
    }
  }

  if(all(is.na(ptcell))) {
    stop("No points could be matched to accessible cells in tardis")
  }
  if(any(as.logical(pmod))) {
    warning("Some points were translated to accessible cells in tardis")
  }
  if(any(is.na(pmod))) {
    warning("Some points could not be translated to accessible cells in tardis and were discarded")
  }

  out <- cbind.data.frame(point = which(!is.na(ptcell)), bin = pt[!is.na(ptcell)] %/% prod(tardis$gdat[1:2]) + 1, mod = pmod[!is.na(ptcell)])
  st_geometry(out) <- st_sfc(lapply(pcell[!is.na(ptcell)], function(x) {st_point(xyFromCell(samprast, x))}), crs = "+proj=longlat")
  return(out)
}
