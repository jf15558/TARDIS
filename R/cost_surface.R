#' cost_surface
#'
#' Convenience function to visualise the weighting scheme for the layers in a
#' tardis object.
#'
#' @param tardis An object of class 'tardis'
#' @param weights A character string denoting the weighting scheme to visualise.
#' By default these are true geographic distances (gdist). Alternatively, the
#' name of a weighting scheme added to the tardis object with weight_tardis().
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return A geoglist with layers corresponding to the original landscape
#' layers in the graph, and cell values calculated as the mean of its inbound edge weights.
#' @import terra sf h3jsr
#' @export
#'
#' @examples
#' \dontrun{
#' #library(terra)
#' #library(TARDIS)
#'
#' # load data
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # make geoglist and add links
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rlink <- link_islands(rasts)
#'
#' rtd <- build_tardis(rasts, times = c(seq(2.25, 0, -0.5), 0))
#' cost_surface(rtd)
#' }

cost_surface <- function(tardis, weights = "gdist", verbose = T) {

  # tardis = rtdw
  # weights = "clim"

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.na(tardis$gdat[7])) {
    grid <- get_grid(tardis$gdat[1:4], tardis$gdat[7])

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
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

  edges <- tardis$edges[which(tardis$edges[,3] == 0),]
  pt <- edges[,2] %/% tardis$gdat[[5]] + 1
  pcell <- edges[,2] %% tardis$gdat[5]

  out <- list()
  for(i in 1:max(pt)) {

    if(verbose) {
      cat(paste0("Calculating surface [", i, "/", max(pt), "]\r"))
      if(i == max(pt)) {cat("\n")}
    }

    # calculate mean weights
    wts <- tapply(edges[which(pt == i), weights], INDEX = pcell[which(pt == i)], mean, na.rm = TRUE)
    wts[!is.finite(wts)] <- NA

    if(!is.na(tardis$gdat[7])) {
      tmp <- cell_to_polygon(grid[as.numeric(names(wts))])
      wts <- data.frame(wts)
      colnames(wts) <- weights
      st_geometry(wts) <- st_wrap_dateline(tmp, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
      tmp <- wts

    } else {
      tmp <- samprast
      tmp[] <- NA
      names(tmp) <- weights
      tmp[][as.numeric(names(wts)),1] <- wts
    }

    # store
    out[[i]] <- tmp
  }
  if(is.na(tardis$gdat[7])) {
    out <- rast(out)
  }
  out <- list(gdat = tardis$gdat, layers = out)
  class(out) <- "geoglist"
  return(out)
}
