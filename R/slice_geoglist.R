#' slice_geoglist
#'
#' Extract layers in a geoglist to their own geoglist, along with any links if present.
#' Unlike with `slice_tardis()`, discontinuous ranges of layers can be selected.
#'
#' @param geog `geoglist`. A geoglist object, with or without links.
#' @param layers `numeric`. The layer index to extract or a range.
#' @return `geoglist`. A subsetted geoglist.
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts1 <- slice_geoglist(rasts, layers = 1)
#' rasts2 <- slice_geoglist(rasts, layers = c(1, 3))
#' }

slice_geoglist <- function(geog, layers) {

  if(!inherits(geog, "geoglist")) {
    stop("geog should be `geoglist` object")
  }
  if(!is.atomic(layers)) {
    stop("layers should be a numeric vector of layer(s) to extract")
  }
  if(!is.numeric(layers)) {
    stop("layers should be a numeric vector of layer(s) to extract")
  }
  if(any(layers %% 1 != 0)) {
    stop("All elements of layers should be positive integers")
  }
  if(any(layers < 1)) {
    stop("All elements of layers should be positive integers")
  }
  nly <- ifelse(inherits(geog$layers, "SpatRaster"), nlyr(geog$layers), length(geog$layers))
  if(any(layers > nly)) {
    stop("One or more elements in layers exceeds the number of available layers in geog")
  }
  geog$layers <- geog$layers[layers]
  if(length(layers) == 1) {
    geog$layers <- svc(geog$layers)
  }
  if(!is.null(geog$links)) {
    geog$links <- geog$links[which(geog$links$layer %in% layers)]
  }
  return(geog)
}
