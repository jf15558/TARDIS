#' resistance_surface
#'
#' Convenience function to visualise the weighting scheme for the layers in a
#' tardis object. Two options are available, one to visualise the empirical
#' geographic distance weights in the landscapes, the other to visualise a
#' custom weighting scheme.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in x. All entries must be >= 0 and finite (no NaN or
#' Inf), or NA. Typically the output of weight_tardis
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return A SpatRaster with each layer corresponding to the original landscape
#' layers in the graph, and cell values calculated as the mean of its inbound edge weights.
#' @import terra
#' @export

resistance_surface <- function(tardis, weights = NULL, verbose = T) {

  # tardis = test
  # weights = NULL

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(class(tardis) != "tardis") {
    stop("Supply tardis as the output of create_tardis")
  }

  if(!is.null(weights)) {

    if(!class(weights)[1] == "vector") {
      stop("weights must be a vector")
    }
    if(!is.numeric(weights)) {
      stop("weights (the vector or chosen column) must be numeric")
    }
    if(any(is.na(weights))) {
      stop("NA entries are not permitted in weights")
    }
    if(any(is.infinite(weights) | weights < 0 | is.nan(weights))) {
      stop("All entries in weights must be positive, finite numerics or NA")
    }
    if(length(weights) != nrow(tardis$edges)) {
      stop("weights must be the same length (if a vector) or have a number of rows (if a data.frame) as the number of graph edges in x")
    }
  } else {
    weights <- tardis$edges[,5]
  }

  samprast <- rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8]))
  pt <- ceiling(tardis$edges[,2] / prod(tardis$gdat[1:2]))
  pcell <- tardis$edges[,2] - ((pt - 1) * prod(tardis$gdat[1:2]))

  out <- list()
  for(i in 1:max(pt)) {

    if(verbose) {
      cat(paste0("Calculating surface [", i, "/", max(pt), "]\r"))
      if(i == max(pt)) {cat("\n")}
    }

    # calculate mean weights
    tmprast <- samprast
    wts <- tapply(weights[which(pt == i)], INDEX = pcell[which(pt == i)], mean, na.rm = TRUE)
    tmprast[as.numeric(names(wts))] <- as.vector(wts)
    tmprast[!is.finite(tmprast[])] <- NA

    # store
    out[[i]] <- tmprast
  }
  return(rast(out))
}
