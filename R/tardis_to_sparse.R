#' tardis_to_sparse
#'
#' Create either a sparse weighted adjacency or transition probability matrix
#' from a tardis graph.
#'
#' @param tardis A tardis graph
#' @param weights `character`. A character string denoting the weighting scheme to use. By
#' default these are true geographic distances (gdist). Alternatively, the name
#' of a weighting scheme added to the tardis object with weight_tardis().
#' @param mode `character`. One of "adjacency" or "transition".
#' @return `matrix` A `Matrix::sparseMatrix` adjacency or transition probability matrix.
#' @import cppRouting Matrix
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' hexes <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' cs <- tardis_to_matrix(htd, mode = "cost")
#' aj <- tardis_to_matrix(htd, mode = "adjacency")
#'}

tardis_to_sparse <- function(tardis, weights = "gdist", mode = "adjacency") {

  #tardis <- h1
  #weights = "gdist"
  #mode = "hitting"

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
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
  if(!is.atomic(mode) | length(mode) != 1) {
    stop("mode should only contain one element")
  }
  if(!is.character(mode) | !mode %in% c("adjacency", "transition")) {
    stop("mode should be one of 'adjacency' or 'transition'")
  }

  tardis <- instantiate_tardis(tardis = tardis, weights = weights)
  if(mode == "adjacency") {
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1,
                        x = tardis$tgraph$data$dist)
  } else {
    # matrix as conductance rather than resistance
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1, x = 1 / tardis$tgraph$data$dist)
    # normalise into probability matrix
    mat <- mat / rowSums(mat)
  }
  return(mat)
}
