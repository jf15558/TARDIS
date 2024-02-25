#' weight_tardis
#'
#' Generate a custom weighting scheme for a TARDIS analysis. Internally,
#' weight_tardis generates two data.frames, origin and dest, which respectively
#' record the landscape properties for each pair of origin and destination
#' cells comprising the edges in a graph layer. The first three columns in both
#' data.frames record the cell ID, the horizontal distance to its partner cell
#' in metres, and the horizontal and vertical distances to its partner cell in
#' metres. The horizontal distance will be the same in both data.frames, but
#' the vertical distance will be positive or negative depending on whether
#' elevation is gained or lost along an edge. All subsequent columns in both
#' data.frames record the cell characteristics as supplied by the user. The
#' columns in these data.frames can then be used to calculate weights with a
#' custom, user-supplied weighting function.
#'
#' Crucially, all returned weights must be finite and greater than zero, or NA,
#' as negative weights are not meaningful for downstream methods and weights
#' of zero are reserved for the interlayer edges. Thus, the mathematics of the
#' function should be carefully considered as it may be very easy to create
#' values of zero based on numeric differences in characteristics between
#' adjacent cells, then produce Inf by zero division. NA values are permitted
#' to allow the designation of impermeable edges (for example to restrict
#' movement above a certain threshold cost). Such values, however, may
#' introduce inaccessible islands into a landscape even after resolution of
#' mask islands by create_tardis, which may lead to certain functions like
#' lcp() failing to return complete traverses between points. A solution to this
#' is to use resistance_surface() to derive raster layers with the positions
#' of cells without accessible edges and include these cells within the
#' initial masking stage. This is a somewhat roundabout solution, but the
#' inclusion of NA weights is preferable for flexible landscape specification.
#'
#' @param tardis A object of class 'tardis' from create_tardis.
#' @param vars A named list of SpatRasters where each element has the same
#' resolution, extent, number of layers and layer order as the original landscape
#' data used to create tardis, and records a single property for that landscape
#' (e.g., temperature). List names are crucial as these are used to make their
#' data accessible during weighting. Any names can be used, aside for 'cell',
#' hdist' and 'vdist' which are reserved.
#' @param wfun A function to calculate the cost of traversal for the edges in
#' each graph layer. This must have the signature:
#' function(origin, dest, lnum, ...) {rules for weighting} - see @details
#' @param mfun A function to calculate the cost of traversal for the edges
#' added to bridge mask islands, should the user desire them to be weighted
#' differently to the result of the later. The function should have the same
#' signature as @param wfun, but its internals can be totally different. NULL
#' by default, so edge weights will come from @param wfun.
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @param ... Additional arguments supplied to wfun() if desired.
#' @return A numeric vector of weights with as many elements as edges in x.
#' @import terra
#' @export
#'
#' @details The core of the weighting function can have as many steps as the
#' user likes, but must consist of vectorised calculations that call on the
#' columns in the data.frames origin and dest. The argument default exemplifies
#' this and will return identical geographic distances to those in within x,
#' using Pythagoras's theorem on the horizontal (origin$hdist) and vertical
#' (origin$vdist) intercell distances. The user function can use one, the other,
#' or any combination of origin and dest columns in any order.
#'
#' Weights are iteratively calculated for each landscape layer in tardis, with
#' the index of the landscape layer internally supplied to the argument lnum to
#' allow the user to design weighting rules that can vary through time. The
#' function can additionally take a dots argument to allow data to be supplied
#' to the weighting function from the global environment, for example an object
#' with elements to be used in conjunction with lnum.
#'
#' The weight_tardis function is heavily inspired by the weighting function
#' used in the gen3sis R package by Oskar Hagen. A major difference between
#' tardis and gen3sis, however, is that the former considers a landscape as a
#' lattice graph where only adjacent cell connections are considered. This
#' produces a sparse distance matrix to which graph algorithms can be applied,
#' rather than the dense distance matrix used by gen3sis, where all intercell
#' connections are given, but without consideration of the intervening space
#' those connections span through.
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
#' #vars = list(elev = classify(gal, cbind(-Inf, 0, 0)))
#' #gtw <- weight_tardis(test2, vars = vars,
#' #                     mfun = function(origin, dest, lnum, ...) {
#' #                               (origin$hdist^2 + abs(origin$vdist)^2) * 10})

weight_tardis <- function(tardis, vars, wfun = function(origin, dest, lnum = NULL, ...) {sqrt(origin$hdist^2 + abs(origin$vdist)^2)}, mfun = NULL, verbose = TRUE, ...) {

  # tardis = ob
  # vars = list(elev = dem)
  # wfun = function(origin, dest, lnum, ...) {sqrt(origin$hdist^2 + abs(origin$vdist)^2)}
  # wfun = dst
  # mfun = NULL
  # verbose = T
  # tardis = ob
  # vars = clim
  # wfun = pcaweight

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.list(vars)) {
    stop("Supply vars as a named list of SpatRasters")
  }

  if(is.null(names(vars))) {
    stop("Supply vars as a named list of SpatRasters")
  }
  if(any(names(vars) %in% c("cell", "hdist", "vdist"))) {
    stop("The names cell, h_dist and v_dist are reserved for the geographic distances in x. Please revise names of vars")
  }
  if(!all(unlist(lapply(vars, inherits, "SpatRaster")))) {
    stop("One or more elements of vars is not a SpatRaster")
  }
  if(!(all(unlist(lapply(vars, is.lonlat))))) {
    stop("One or more elements of vars is not in geographic (long-lat) projection")
  }

  if(!is.null(tardis$tdat)) {
    layers <- length(tardis$tdat) - 1
  } else {
    layers <- 1
  }
  if(!all(unlist(lapply(vars, function(x) {dim(x) == c(tardis$gdat[1:2], layers)})))) {
    stop("Each SpatRaster in vars must match the extent, resolution and number of layers as the SpatRaster used to create tardis")
  }
  if(!all(unlist(lapply(vars, function(x) {as.vector(ext(x)) == tardis$gdat[5:8]})))) {
    stop("Each SpatRaster in vars must match the extent, resolution and number of layers as the SpatRaster used to create tardis")
  }
  if(any(unlist(lapply(vars, function(x) {any(is.na(x[]))})))) {
    stop("NA values are not permitted in vars. Please use dummy values if needed")
  }
  if(!is.function(wfun)) {
    stop("wfun should be a user-supplied function. See documentation for required function signature")
  }
  if(!is.null(mfun)) {
    if(!is.function(mfun)) {
      stop("mfun should be a user-supplied function. See documentation for required function signature")
    }
  }

  src <- ceiling(tardis$edges[,1] / prod(tardis$gdat[1:2]))
  dst <- ceiling(tardis$edges[,2] / prod(tardis$gdat[1:2]))
  wts <- rep(0, length(src))

  for(i in 1:layers) {

    if(verbose) {
      cat(paste0("Weighting layer [", i, "/", layers, "]\r"))
      if(i == layers) {cat("\n")}
    }

    links <- tardis$edges[which(src == i & src == dst),]
    links[,1:2] <- links[,1:2] - ((i - 1) * prod(tardis$gdat[1:2]))
    origin <- lapply(vars, function(y) {y[[i]][links[,1]]})
    origin <- cbind.data.frame(links[,1], links[,3:4], do.call(cbind.data.frame, origin))
    dest <- lapply(vars, function(y) {y[[i]][links[,2]]})
    dest <- cbind.data.frame(links[,2], links[,3:4], do.call(cbind.data.frame, dest))
    colnames(origin) <- colnames(dest) <- c("cell", "hdist", "vdist", names(vars))

    weight <- try(wfun(origin = origin, dest = dest))
    if(class(weight)[1] == "try-error") {
      stop(paste0("An error occurred in wfun() for  layer ", i, "/", layers, ". Check that the column names in wfun() match the names of vars, along with 'hdist' and 'vdist'"))
    }
    if(!is.vector(weight) | length(weight) != nrow(links)) {
      stop(paste0("wfun() did not return a vector with as many elements as edges in layer ", i, "/", layers, ". Ensure that the function returns a vector of correct length"))
    }
    if(any(is.nan(weight) | is.infinite(weight))) {
      stop(paste0("wfun() resulted in a non-finite value (NaN, Inf) in  layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers or NA"))
    }
    if(any(weight < 0)) {
      stop(paste0("wfun() resulted in a negative value in layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers"))
    }

    mlink <- which(!tardis$gdat[2] %% abs(links[,1] - links[,2]) %in% c(0, 1, tardis$gdat[2]))
    if(!is.null(mfun) & length(mlink) != 0) {

      mweight <- try(mfun(origin = origin[mlink,], dest = dest[mlink,], lnum = i))
      if(class(mweight)[1] == "try-error") {
        stop(paste0("An error occurred in mfun() for  layer ", i, "/", layers, ". Check that the column names in mfunc() match the names of vars, along with 'hdist' and 'vdist'"))
      }
      if(!is.vector(mweight) | length(mweight) != length(mlink)) {
        stop(paste0("mfun() did not return a vector with as many elements as edges in layer ", i, "/", layers, ". Ensure that the function returns a vector of correct length"))
      }
      if(any(is.nan(weight) | is.infinite(weight))) {
        stop(paste0("mfun() resulted in a non-finite value (NaN, Inf) in  layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers or NA"))
      }
      if(any(mweight < 0)) {
        stop(paste0("mfun() resulted in a negative value in layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers"))
      }
      weight[mlink] <- mweight
    }
    wts[which(src == i & src == dst)] <- weight
  }

  if(any(is.na(wts))) {
    warning("Some weights are NA. This is permissible, but may produce inaccessible islands. Consider checking with resistance_surface()")
  }

  # return weighting vector
  return(wts)
}
