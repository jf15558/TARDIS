#' instantiate_tardis
#'
#' Take a tardis object and generate a graph compatible with functions from the
#' `rcppRouting` package. This function is primarily called internally by other
#' TARDIS package functions, but it may be useful for the user to be able to
#' directly interact with this graph themselves, or to supply a pre-instantiated
#' `tardis` object within a loop to avoid repeated calls to this function.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param weights `character`. The name of the weighting scheme column in
#' `tardis$edges` to use. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of a weighting scheme added to the tardis
#' object with `weight_tardis()`.
#' @return. The input `tardis` object with the additional element, `tgraph`.
#' @export
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- galapagos()
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # create a geoglist in raster format and mask the sea
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts)
#'
#' # build a tardis from raster cells
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#'
#' # instantiate
#' rtd <- instantiate_tardis(rtd)
#' }

instantiate_tardis <- function(tardis, weights = "gdist") {

  if(!exists("tardis")) {
    stop("Supply tardis as a tardis object")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as a tardis object")
  }

  if(!is.atomic(weights) | length(weights) != 1) {
    stop("weights should only contain one element")
  }
  if(!is.character(weights)) {
    stop("weights should be a character string")
  }
  if(!weights %in% colnames(tardis$edges)) {
    stop("weights should be a column name in tardis$edges")
  }

  if(!is.numeric(tardis$edges[,weights])) {
    stop("Some weights are not numeric")
  }
  if(any(is.infinite(tardis$edges[,weights]) | tardis$edges[,weights] < 0 | is.nan(tardis$edges[,weights]))) {
    stop("All weights must be positive, finite numerics, NA, or zero")
  }

  weights <- tardis$edges[,weights]
  if(any(is.na(weights))) {
    tardis$tgraph$src <- tardis$tgraph$src[which(!is.na(weights))]
    tardis$tgraph$dst <- tardis$tgraph$dst[which(!is.na(weights))]
    tardis$edges <- tardis$edges[which(!is.na(weights)),]
    tardis$tgraph$dict <- tardis$tgraph$dict[which(tardis$tgraph$dict$ref %in%
                                                     c(tardis$tgraph$src, tardis$tgraph$dst))]
    tardis$tgraph$id <- 0:(length(tardis$tgraph$dict$ref) - 1)
    tardis$tgraph$nbnode <- length(tardis$tgraph$id)
  }
  tardis$tgraph$attrib$aux <- tardis$edges[!is.na(weights), "gdist"]
  tardis$tgraph$data <- data.frame(from = tardis$tgraph$dict$id[match(tardis$tgraph$src,
                                                                      tardis$tgraph$dict$ref)], to = tardis$tgraph$dict$id[match(tardis$tgraph$dst,
                                                                                                                                 tardis$tgraph$dict$ref)], dist = weights[!is.na(weights)])
  tardis$tgraph <- tardis$tgraph[1:5]
  return(tardis)
}
