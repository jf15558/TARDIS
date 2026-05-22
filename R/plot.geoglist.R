#' plot.geoglist
#' 
#' Plotting method for geoglist. If the geoglist contains multiple layers, then
#' default behaviour is to plot the first one.
#' 
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param layer `numeric`. The layer in the geoglist to be plotted, along with
#' its links. Defaults to 1 (the first layer).
#' @param pal `vector`. A vector of colours to be used for plotting layer values,
#' such as those returned by an R colour palette.
#' @param links `logical`. Should mask links be plotted, if available?
#' @param lcol `character` or `integer`. The colour to be used for plotting links
#' @param lwd `integer`. The line width to be used for plotting links.
#' @param lty `integer`. The line type to be used for plotting links.
#' @param hex.border `character` or `integer`. The colour to be used for plotting
#' hexagonal grids. By default none.
#' @return None.

plot.geoglist <- function(geog, layer = 1, pal = terrain.colors(100), links = T,
                          lcol = 1, lwd = 1, lty = 1, hex.border = NA) {
  
  if(!exists("geog")) {
    stop("Supply geog as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  
  if(!is.atomic(layer) | length(layer) != 1) {
    stop("layer should be a single integer")
  }
  if(!is.numeric(layer)) {
    stop("layer should be a single integer")
  }
  if(layer %% 1 != 0) {
    stop("layer should be a single integer")
  }
  if(layer > length(geog$layers)) {
    stop("The value of layer exceeds of the number of layers in geog")
  }
  
  if(inherits(geog$layers[[1]], "SpatRaster")) {
    plot(geog$layers[[bin]], col = pal, legend = F)
  } else {
    plot(geog$layers[[bin]]$geometry, border = NA)
    plot(geog$layers[[bin]][,1], add = T, col = pal, border = hex.border)
  }
  if(links) {
    if(!is.null(geog$links)) {
      plot(geog$links[which(geog$links$layer == bin),"geometry"], add = T,
           col = lcol, lwd = lwd, lty = lty)
    }
  }
}