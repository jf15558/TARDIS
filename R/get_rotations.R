#' get_rotations
#'
#' Generate a rotation list using commonly-implemented plate rotation models for
#' use with create_tardis(). The user supplies a SpatRaster with the desired
#' geographic resolution for the resulting rotation list, along with a vector of
#' reconstruction ages in descending order. Conceivably, the SpatRaster and age
#' vector the user intends to supply to create_tardis() could be used here. The
#' function will take the cell centre coordinates for all grid cells in the
#' SpatRaster as starting positions at the first age in the vector, then
#' calculates their geographic positions at the subsequent ages. Calculation
#' relies on the reconstruct() function from the rgplates package and will be
#' very slow if high resolution SpatRasters/lengthy time vectors are supplied.
#' Each SpatRaster layer must have at least some reconstructed grid cell
#' positions so the age vector cannot be made shorter. Instead, it may be
#' advisable to use a SpatRaster with a lower geographic resolution than those
#' the user intends to supply to create_tardis(). This is totally fine as the
#' rotations list does not need to give reconstructions for every grid cell in
#' the landscape raster stack, nor will this even be possible as gplates will be
#' unable to determine some geographic positions due to the limitations of the
#' plate rotation models themselves. NAs are returned where reconstructions fail,
#' although this is again not an issue as these will be automatically filtered
#' by create_tardis(). There is no point performing rotations at a higher
#' geographic resolution than the SpatRaster stack that you will use with
#' create_tardis()
#' #'
#' @param rast A SpatRaster in longitude-latitude projection
#' @param ages A numeric vector of desired reconstruction ages in millions of
#' years before present, in descending age order. These will correspond to the
#' layers in the SpatRaster destined for create_tardis().
#' @param model The desired plate reconstruction model. See
#' rgplates::reconstruct() for details.
#' @param out.res If you are calculating rotations at a reduced resolution to
#' save time, then out.res should be the SpatRaster you intent to use with
#' create_tardis(). This is essential so that the returned cell IDs will conform
#' to that SpatRaster, rather than to that given by 'rast'.
#' @param gpath If NULL, the GPlates API is used. This is very slow for large
#' numbers of points, so the software call be installed locally for faster
#' results. If desired gpath is the file path to the local reconstruction
#' submodule. In case the GPlates executable file is not found at the coded
#' default location, the full path to the executable (gplates-<ver>.exe on
#' Windows) can be entered, e.g. "C:/gplates_2.3.0_win64/gplates.exe".
#' @param verbose Should function progress be reported?
#' @return A list with length(ages) - 1 elements. Each element is a two-column
#' matrix containing the IDs of geographically homologous cells between
#' successive reconstruction ages.
#' @import terra rgplates
#'
#' @examples
#' #rast <- rast(nrows = 2, ncols = 2)
#' #ages <- c(440, 430, 420)
#' #foo <- get_rotations(rast = rast, ages = ages)
#'
get_rotations <- function(rast, ages, model = "PALEOMAP", out.res = NULL, gpath = NULL, verbose = TRUE) {

  # rast <- rast(nrows = 2, ncols = 2)
  # ages <- c(440, 430, 420)
  # model = "PALEOMAP"
  # out.res = NULL
  # gpath = NULL
  # verbose = TRUE

  if(!exists("rast")) {
    stop("Supply rast as a SpatRaster with longitude-latitude projection")
  }
  if(!inherits(rast, "SpatRaster")) {
    stop("Supply rast as a SpatRaster with longitude-latitude projection")
  }
  if(!is.lonlat(rast)) {
    stop("Supply rast as a SpatRaster with longitude-latitude projection")
  }

  if(!exists("ages")) {
    stop("Supply ages as a vector of descending reconstruction ages in millions of years before present")
  }
  if(!inherits(ages, "numeric") | !is.vector(ages)) {
    stop("Supple ages as a vector of descending reconstruction ages in millions of years before present")
  }
  if(any(is.na(ages)) | any(is.infinite(ages)) | any(ages < 0)) {
    stop("ages cannot contain missing, infinite or negative values")
  }
  if(length(ages) < 2) {
    stop("ages must contain at least two elements")
  }

  if(!exists("out.res")) {
    stop("Supply out.res as a SpatRaster with longitude-latitude projection")
  }
  if(is.null(out.res)) {out.res <- rast}
  if(!inherits(out.res, "SpatRaster")) {
    stop("Supply out.res as a SpatRaster with longitude-latitude projection")
  }
  if(!is.lonlat(out.res)) {
    stop("Supply out.res as a SpatRaster with longitude-latitude projection")
  }

  if(!exists("model")) {
    stop("Supply model as the character string of the desired plate rotation model (see rgplates::reconstruct for details)")
  }
  if(!inherits(model, "character") | !is.vector(model)) {
    stop("Supply model as the character string of the desired plate rotation model (see rgplates::reconstruct for details)")
  }
  if(length(model) > 1) {
    warning("Only the first element of model will be used")
  }

  # set up args
  model <- model[1]
  ages <- ages[order(ages, decreasing = T)]
  recon <- list(xyFromCell(rast, 1:ncell(rast)))

  # get present day positions from the first time slice positions
  pres <- reconstruct(recon[[1]], age = ages[1], reverse = T, model = model, path.gplates = gpath, verbose = verbose)

  # get the present day positions at the subsequent reconstruction points
  past <- reconstruct(pres, age = ages[-1], model = model, path.gplates = gpath, verbose = verbose)

  # convert from coordinates to cell IDs
  for(i in 1:length(past)) {recon[i + 1] <- past[i]}
  recon <- lapply(recon, function(x) {cellFromXY(out.res, x)})

  # format and return
  recon <- do.call(cbind, recon)
  out <- lapply(1:(ncol(recon) - 1), function(x) {recon[,x:(x + 1)]})
  names(out) <- paste0(ages[-length(ages)], "-", ages[-1])
  return(out)
}

