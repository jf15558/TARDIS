#' tardis_to_matrix
#'
#' Create either a weighted adjacency matrix (sparse), or a cost matrix,
#' transition probability matrix, or a hitting time matrix (dense) from a TARDIS graph.
#'
#' @param tardis A tardis graph
#' @param weights `character`. A character string denoting the weighting scheme to use. By
#' default these are true geographic distances (gdist). Alternatively, the name
#' of a weighting scheme added to the tardis object with weight_tardis().
#' @param mode `character`. One of "adjacency", "cost", "transition" or "hitting".
#' @return `matrix` Either a `Matrix::sparseMatrix` adjacency matrix or a dense
#' base R distance matrix.
#' @import Matrix cppRouting
#' @importFrom markovchain meanFirstPassageTime
#' @importFrom methods new
#' @importClassesFrom markovchain markovchain
#' @export
#'
#' @details
#' For dense matrix methods, this can very quickly create large objects. For
#' example, a TARDIS graph with ~20,000 accessible cells will give a 20k x 20k
#'  matrix approaching a gigabyte in size. You may wish to optimise the resolution
#'  and masking of your landscapes first.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' hexes <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = T, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' dm <- tardis_to_matrix(htd, mode = "cost")
#' aj <- tardis_to_matrix(htd, mode = "adjacency")
#' tr <- tardis_to_matrix(htd, mode = "transition")
#' ht <- tardis_to_matrix(htd, mode = "hitting")
#'}

tardis_to_matrix <- function(tardis, weights = "gdist", mode = "adjacency") {

  # tardis <- slice_tardis(tardis, c(2, 2))
  # weights = "tdist"

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
  if(!is.character(mode) | !mode %in% c("adjacency", "cost", "transition", "hitting")) {
    stop("mode should be one of 'adjacency', 'cost', 'transition' or 'hitting'")
  }

  tardis <- instantiate_tardis(tardis = tardis, weights = weights)
  if(mode == "adjacency") {
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1,
                        x = tardis$tgraph$data$dist)
  }
  if(mode == "cost") {
    mat <- get_distance_matrix(tardis$tgraph, tardis$tgraph$dict$ref, tardis$tgraph$dict$ref)
  }
  if(mode %in% c("transition" ,"hitting")) {
    # matrix as conductance rather than resistance
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1, x = 1 / tardis$tgraph$data$dist)
    # normalise into probability matrix
    mat <- mat / rowSums(mat)

    if(mode == "hitting") {
      # create the Markov Chain object
      mc <- new("markovchain", states = tardis$tgraph$dict$ref, transitionMatrix = mat)

      # calculate the hitting times matrix
      mat <- meanFirstPassageTime(mc)
    }
  }
  return(mat)
}
