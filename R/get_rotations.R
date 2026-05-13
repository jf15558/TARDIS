#' get_rotations
#'
#' Generate a rotation list from commonly-used plate rotation models for
#' use with build_tardis(). Ensure that the same plate rotation model on which
#' the input `geoglist` is based is requested for function call.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`. This should contain
#' at least two layers, representing two palaeogeographic time slices where
#' there is appreciable intervening plate roation.
#' @param times `numeric`. A  vector with `nlayers(geog) + 1` positive
#' elements, expressing the temporal boundaries of each layer as millions of
#' years in the past. The vector need not end in the present (i.e. `0`), but time
#' must flow from oldest to youngest.
#' @param model `character`. The desired plate reconstruction model. See `palaeoverse::palaeorotate()`
#' for details.
#' @param method `character`. The reconstruction method to use. See `palaeoverse::palaeorotate()`
#' for details.
#' @param verbose `logical` Should function progress be reported to the user?
#' @param ... Other arguments passed to `palaeoverse::palaeorotate()`.
#' @return A list with as many elements as layers in `geog$layers`. Each element
#' is a two-column matrix containing the IDs of geographically homologous cells
#' between successive landscape layers.
#' @import terra sf h3jsr palaeoverse
#' @export
#'
#' @details
#' Cell reconstruction from layer to layer is performed by `palaeoverse::palaeorotate()`.
#' By default this uses the `"grid"` method, by retrieving reconstructed coordinates
#' from pre-rotated H3 grids at resolution 3 (~119 km). This is almost
#' certainly sufficient for any global scale analyses given that palaeogeographic
#' models of topography and bathymetry rarely if ever exceed 1 degree resolution
#' (~111 km).
#'
#' @examples
#' \dontrun{
#' rast <- rast(nrows = 2, ncols = 2)
#' ages <- c(440, 430, 420)
#' foo <- get_rotations(rast = rast, times = ages, model = "MERDITH2021")
#'}

get_rotations <- function(geog, times, model, method = "grid", verbose = TRUE, ...) {

  # geog <- rasts
  # times <- c(125, 120)
  # model = "PALEOMAP"
  # method = "grid"
  # verbose = TRUE

  if(!exists("geog")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }

  if(!exists("times")) {
    stop("Supply times as a vector of descending reconstruction ages in millions of years before present")
  }
  if(!is.numeric(times) | !is.vector(times)) {
    stop("Supply times as a vector of descending reconstruction ages in millions of years before present")
  }
  if(any(is.na(times)) | any(is.infinite(times)) | any(times < 0)) {
    stop("times cannot contain missing, infinite or negative values")
  }
  if(length(times) < 2) {
    stop("times must contain at least two elements")
  }
  if (!is.numeric(times) | length(times) != length(geog$layers) + 1) {
    stop("Please supply times as a vector of time bin boundaries with n elements in geog$layers + 1")
  }

  if(!exists("model")) {
    stop("Supply model as the character string of the desired plate rotation model (see palaeoverse::palaeorotate for details)")
  }
  if(!inherits(model, "character") | !is.vector(model)) {
    stop("Supply model as the character string of the desired plate rotation model (see palaeoverse::palaeorotate for details)")
  }
  if(length(model) > 1) {
    warning("Only the first element of model will be used")
  }

  # area for point reconstruction
  if(!is.na(geog$gdat[7])) {
    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
    pts <- st_coordinates(cell_to_point(grid))

  } else {
    grid <- rast(nrows = geog$gdat[5] / geog$gdat[6], ncols = geog$gdat[6],
                 ext = ext(geog$gdat[1:4]))
    pts <- xyFromCell(grid, 1:ncell(grid))
  }
  colnames(pts) <- c("lng", "lat")

  # times to reconstruct (midpoints of the time bins)
  mids <- times[-1] - diff(times) / 2

  # reconstruct
  rots <- lapply(1:length(mids), function(x) {
    if(verbose) {
      cat(paste0("Rotating layer [", x,"/", length(mids), "]\r"))
      if(x == length(mids)) {cat("\n")}
    }
    pts2 <- cbind.data.frame(pts, age = mids[x])
    pts2 <- suppressWarnings(palaeorotate(pts2, model = model, method = method))[,c("p_lng", "p_lat")]

    cls <- rep(NA, nrow(pts2))
    to_cell <- which(!is.na(pts2[,1]))
    if(!is.na(geog$gdat[7])) {
      cls[to_cell] <- match(suppressMessages(point_to_cell(pts2[to_cell,], res = geog$gdat[7])), grid)
    } else {
      cls[to_cell] <- cellFromXY(grid, pts2[to_cell,])
    }
    cls
  })
  rots <- cbind(do.call(cbind, rots), 1:geog$gdat[5])

  # format and return
  out <- lapply(1:(ncol(rots) - 1), function(x) {
    vals <- rots[,x:(x + 1)]
    vals[complete.cases(vals),]
  })
  names(out) <- paste0(times[-length(times)], "-", times[-1])
  return(out)
}

