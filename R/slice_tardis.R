#' slice_tardis
#'
#' Temporally subset a tardis graph to a range of layers.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param times `numeric`. A vector of two positive numbers denoting the desired time range
#' to subset from tardis. Either this or 'layers' must be specified
#' @param layers `numeric`. A vector of two positive numbers denoting the desired range of
#' layers to subset from tardis. Either this or 'times' must be specified
#' @return A `tardis` object comprising the requested layers
#' @export
#'
#' @details
#' This function was developed for instances where successive analyses do not require
#' the entire graph, in which case it is more efficient to weight and analyse subsets,
#' rather than operate on the entire graph or create subsets from scratch. Note that
#' layers are counted in decreasing age order, so the oldest time layer will be 1
#' and so forth.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#'
#' gts <- slice_tardis(rtd, times = c(1.2, 0))
#' gts <- slice_tardis(rtd, layers = c(1, 2))
#' }

slice_tardis <- function(tardis, times = NULL, layers = NULL) {

  # tardis <- rtd
  # times <- NULL
  # layers <- c(2, 2)

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(is.null(tardis$tdat)) {
    stop("Temporal subsetting can only be applied to TARDIS graphs with multiple layers")
  }

  # check subsetting conflict
  if(is.null(times) & is.null(layers)) {
    stop("One of times or layers must be not be NULL")
  }
  if(!is.null(times) & !is.null(layers)) {
    stop("One of times or layers must be left as NULL")
  }

  # check times
  if(!is.null(times)) {

    if(!is.numeric(times) | length(times) != 2) {
      stop("Please supply times as a vector of 2 numbers denoting the desired temporal subset of the graph")
    }
    if(any(is.na(times))) {
      stop("times cannot contain NA values")
    }
    if(any(times < 0)) {
      stop("times must only contain values >= 0")
    }
    times <- times[order(times, decreasing = T)]
    if(times[1] > tardis$tdat[1] | times[2] < tardis$tdat[length(tardis$tdat)]) {
      stop("times must fall within the temporal range of the TARDIS graph")
    }
    layers <- c(sum(times[1] <= tardis$tdat), sum(times[2] < tardis$tdat))
  }

  # check layers
  if(!is.null(layers)) {

    if(!is.numeric(layers) | length(layers) != 2) {
      stop("Please supply layers as a vector of 2 integers denoting the range of desired graph layers")
    }
    if(any(is.na(layers))) {
      stop("layers cannot contain NA values")
    }
    if(any(layers < 0) | any(layers %% 1 != 0)) {
      stop("layers must only contain positive integers")
    }
    layers <- layers[order(layers, decreasing = F)]
    if(layers[1] > length(tardis$tdat) - 1 | layers[2] > length(tardis$tdat) - 1) {
      stop("The values in layers cannot exceed the number of layers in the TARDIS graph")
    }
  }

  # get the cell id range for the requested layer range
  cls <- as.character(((layers[1] * tardis$gdat[5]) - tardis$gdat[5] + 1):(layers[2] * tardis$gdat[5]))

  # subset edges and dict
  valid <- tardis$tgraph$src %in% cls & tardis$tgraph$dst %in% cls
  tardis$tgraph$src <- tardis$tgraph$src[valid]
  tardis$tgraph$dst <- tardis$tgraph$dst[valid]
  tardis$edges <- tardis$edges[valid,]
  tardis$tgraph$dict <- tardis$tgraph$dict[which(tardis$tgraph$dict$ref %in% c(tardis$tgraph$src, tardis$tgraph$dst)),]

  # adjust cell id parameters
  tardis$edges[,1:2] <- tardis$edges[,1:2] - (as.numeric(cls[1]) - 1)
  tardis$tdat <- tardis$tdat[layers[1]:(layers[2] + 1)]
  tardis$tgraph$nbnode <- nrow(tardis$tgraph$dict)
  tardis$tgraph$dict$id <- (1:tardis$tgraph$nbnode) - 1
  tardis$tgraph$dict$ref <- as.character(as.numeric(tardis$tgraph$dict$ref) - (as.numeric(cls[1]) - 1))
  tardis$tgraph$src <- as.character(as.numeric(tardis$tgraph$src) - (as.numeric(cls[1]) - 1))
  tardis$tgraph$dst <- as.character(as.numeric(tardis$tgraph$dst) - (as.numeric(cls[1]) - 1))

  # return
  return(tardis)
}
