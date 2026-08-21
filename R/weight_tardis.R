#' weight_tardis
#'
#' Generate a custom weighting scheme for a `tardis` graph This may be based
#' on the properties in `tardis` itself or derived from `geoglist` objects recording
#' alternative properties of a landscape.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param name `character`. The name for weighting scheme to be generated. The names 'cell', 'type', 'layer', 'bearing',
#' 'hdist', 'vdist' and gdist' are reserved.
#' @param vars `list`. A named list of `geoglist` objects, each recording an alternative
#' geographic property of the landscape represented in `tardis`. As such, these
#' must bear the same resolution, extent, number of layers and layer order as the
#' `geoglist` used to create `tardis`. The names 'cell', 'type', 'layer', 'bearing',
#' 'hdist', 'vdist' and gdist' are reserved.
#' @param wfun `function(origin, dest)` A function to calculate the
#' cost of traversal for the edges in each graph layer. See @details for the required
#' function signature.
#' @param verbose `logical`. Should function progress be to the user?
#' @return The input `tardis` object with the new weighting scheme added to
#' `tardis$edges` under the column name given in `weights`.
#' @import terra
#' @export
#'
#' @details
#' The `weight_tardis` function is heavily inspired by the weighting function
#' used in the `gen3sis` R package by Oskar Hagen.
#'
#' Internally, `weight_tardis()` generates two `data.frames`, `origin` and `dest`.
#' These minimally record the properties for each pair of origin and destination cells
#' comprising the edges in a graph layer. Each records the cell ID, the type of edge it forms, the horizontal
#' distance and bearing to its partner cell in metres and degrees respectively,
#' the vertical distance to its partner cell in metres, and the time layer in which that cell resides. Horizontal distance will
#' be the same in `origin` and `dest`, but the vertical distance will be positive
#' or negative depending on whether elevation is gained or lost along an edge.
#'
#' If `vars` are supplied, then subsequent columns in both data.frames record the
#' cell characteristics in each `geoglist`. The columns in these data.frames can
#' then be used to calculate weights with a custom, user-supplied weighting function.
#' If the tardis graph contains mask links, these can be weighted differently if
#' desired.
#'
#' The core of the weighting function can have as many steps as the
#' user likes, but must consist of vectorised calculations that call on the
#' columns in the data.frames `origin` and `dest`. The argument default exemplifies
#' this and will return identical geographic distances to those in within x,
#' using Pythagoras's theorem on the horizontal (origin$hdist) and vertical
#' (origin$vdist) intercell distances. The user function can use one, the other,
#' or any combination of origin and dest columns in any order.
#'
#' Weights are iteratively calculated for each landscape layer in tardis, with
#' the index of the landscape layer available in each data.frame to allow the
#' user to design weighting rules that can vary through time.
#'
#' Crucially, all returned weights should be finite and > 0, as negative
#' weights are not permitted for downstream functions and 0 is reserved for temporal links.
#' `NA` values are permitted to allow the designation of impermeable edges (for example to restrict
#' movement above a certain threshold cost). Such values, however, may introduce
#' inaccessible islands into a landscape even if island linkage has already been
#' performed.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(rTARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- galapagos()
#' gal <- crop(gal, ext(-92, -88, -2, 1))
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # create a geoglist with hexagonal resampling and mask the sea
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # create a dummy alternative raster dataset with random values
#' gal_v <- rast(gal_m, vals = rnorm(abs(ncell(gal_m) * nlyr(gal_m))))
#'
#' # convert the dummy data to a geoglist with the same masking
#' gvars <- rast_to_geoglist(gal_v, gal_m)
#'
#' # define a function where geographic distances are multiplied by 10
#' altfunc <- function(origin, dest, tnum, ...) {
#'   (origin$hdist^2 + abs(origin$vdist)^2) * 10
#' }
#'
#' # place the alternative geoglist data into a list
#' vrs <- list(rand = gvars)
#'
#' # generate a new weighting scheme where geographic distances along mask links are multiplied by 10
#' gtw <- weight_tardis(rts, name = "altweight", vars = vrs, mfun = altfunc())
#' }

weight_tardis <- function(tardis, name, vars = NULL, wfun, verbose = TRUE) {


   # tardis = rtd
   # name = "elev"
   # vars = wvars
   # wfun = wfun


  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.null(tardis$tdat)) {
    layers <- length(tardis$tdat) - 1
  } else {
    layers <- 1
  }

  if(!exists("name")) {
    stop("Supply please supply a name for the weighting scheme to be created")
  }
  if(!is.atomic(name) | length(name) != 1) {
    stop("name should only contain one element")
  }
  if(!is.character(name)) {
    stop("name should be a character string")
  }
  if(name %in% c("from", "to", "cell", "type", "bearing", "hdist", "vdist", "gdist")) {
    stop("The names from, to, cell, type, bearing, hdist, vdist and gdist are reserved for tardis internals. Please revise choose a different value for name")
  }

  if(!is.null(vars)) {
    if(!is.list(vars)) {
      stop("Supply vars as a named list of geoglists")
    }
    if(is.null(names(vars))) {
      stop("Supply vars as a named list of geoglists")
    }
    if(any(names(vars) %in% c("from", "to", "cell", "type", "bearing", "hdist", "vdist", "gdist"))) {
      stop("The names from, to, cell, type, bearing, hdist, vdist and gdist are reserved for tardis internals. Please revise choose a different value for name")
    }
    if(!all(unlist(lapply(vars, inherits, "geoglist")))) {
      stop("One or more elements of vars is not a geoglist")
    }
    if(layers > 1) {
      if(!all(sapply(lapply(vars, `[[`, "gdat"), function(x) {all.equal(tardis$gdat, x)}))) {
        stop("Each geoglist in vars must share the extent, resolution and number of layers as tardis")
      }
    }
  }

  if(!exists("wfun")) {
    stop("wfun should be a user-supplied function. See documentation for required function signature")
  }

  if(!is.function(wfun)) {
    stop("wfun should be a user-supplied function. See documentation for required function signature")
  }

  src <- ceiling(tardis$edges[,1] / tardis$gdat[5])
  dst <- ceiling(tardis$edges[,2] / tardis$gdat[5])
  wts <- rep(0, length(src))

  for(i in 1:layers) {

    if(verbose) {
      cat(paste0("Weighting layer [", i, "/", layers, "]\r"))
      if(i == layers) {cat("\n")}
    }

    links <- tardis$edges[which(src == i & src == dst),]
    links[,1:2] <- links[,1:2] %% tardis$gdat[5]
    origin <- as.data.frame(links[,c(1, 3, 4:6)])
    dest <- as.data.frame(links[,c(2, 3, 4:6)])
    origin <- cbind.data.frame(origin, (origin[,1] %/% tardis$gdat[5] + 1))
    dest <- cbind.data.frame(dest, (dest[,1] %/% tardis$gdat[5] + 1))
    origin <- origin[,c(1, 6, 2, 3, 4, 5)]
    dest <- dest[,c(1, 6, 2, 3, 4, 5)]
    if(!is.null(vars)) {
      if(inherits(vars[[1]]$layers, "SpatRaster")) {
        vrs <- lapply(vars, function(y) {y$layers[[i]][][links[,1],]})
        origin <- cbind.data.frame(origin, vrs)
        vrs <- lapply(vars, function(y) {y$layers[[i]][][links[,2],]})
        dest <- cbind.data.frame(dest, vrs)
      } else {
        vrs <- lapply(vars, function(y) {y$layers[[i]][[1]][match(links[,1], y$layers[[i]]$id)]})
        origin <- cbind.data.frame(origin, vrs)
        vrs <- lapply(vars, function(y) {y$layers[[i]][[1]][match(links[,2], y$layers[[i]]$id)]})
        dest <- cbind.data.frame(dest, vrs)
      }
    }
    colnames(origin) <- colnames(dest) <- c("cell", "layer", "type", "bearing", "hdist", "vdist", names(vars))

    weight <- try(wfun(origin = origin, dest = dest))
    if(class(weight)[1] == "try-error") {
      stop(paste0("An error occurred in wfun() for  layer ", i, "/", layers, ". Check that the column names in wfun() match the names of vars, along with 'hdist' and 'vdist'"))
    }
    weight <- as.vector(weight)
    if(!is.vector(weight) | length(weight) != nrow(links)) {
      stop(paste0("wfun() did not return a vector with as many elements as edges in layer ", i, "/", layers, ". Ensure that the function returns a vector of correct length"))
    }
    if(any(is.nan(weight) | is.infinite(weight))) {
      stop(paste0("wfun() resulted in a non-finite value (NaN, Inf) in  layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers or NA"))
    }
    if(all(is.na(weight))) {
      stop(paste0("wfun() resulted in all NA weights in layer ", i, "/", layers))
    }
    if(any(na.omit(weight) <= 0)) {
      stop(paste0("wfun() resulted in a non-positive value in layer ", i, "/", layers, ". Ensure the function and data returns positive real numbers or NA"))
    }
    wts[which(src == i & src == dst)] <- weight
  }

  if(any(is.na(wts))) {
    warning("Some weights are NA. This is permissible, but may produce inaccessible islands. Consider checking with cost_surface()")
  }

  # attach weighting scheme and return
  if(name %in% colnames(tardis$edges)) {
    tardis$edges[,name] <- wts
  } else {
    tardis$edges <- cbind(tardis$edges, wts)
    colnames(tardis$edges)[ncol(tardis$edges)] <- name
  }
  return(tardis)
}
