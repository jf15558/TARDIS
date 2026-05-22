#' plot.geoglist
#'
#' Plotting method for a geoglist layer. If the geoglist contains multiple
#' layers, then default behaviour is to plot the first one.
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
#' @param legend `logical`. Should a legend be added to the plot? Defaults to
#' `TRUE`.
#' @return None.
#' @import sf terra
#' @export

plot.geoglist <- function(geog, layer = 1, pal = sf.colors(10), links = T,
                          lcol = 1, lwd = 1, lty = 1, hex.border = NA, legend = T) {

  # geog = foo
  # layer = 1
  # pal = sf.colors(10)
  # links = T
  # lcol = 1
  # lwd = 1
  # lty = 1
  # hex.border = NA
  # legend = T

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

  if(legend) {
    pr <- par()
    newmar <- pr$mar
    newmar[4] <- 5.1
    par(mar = newmar)
  }
  plot(0, 0, type = "n", xlim = geog$gdat[1:2], ylim = geog$gdat[3:4],
       asp = 1, xaxs = "i", yaxs = "i", xlab = "", ylab = "", cex.axis = 0.7)
  if(inherits(geog$layers[[1]], "SpatRaster")) {
    plot(geog$layers[[layer]], col = pal, legend = F, add = T)
  } else {
    plot(geog$layers[[layer]][,1], add = T, pal = pal, border = hex.border)
    box()
  }
  if(links) {
    if(!is.null(geog$links)) {
      plot(geog$links[which(geog$links$layer == layer),"geometry"], add = T,
           col = lcol, lwd = lwd, lty = lty)
    }
  }
  if(legend) {
    if(inherits(geog$layers[[1]], "SpatRaster")) {
      rng <- c(minmax(geog$layers[[layer]]))
    } else {
      rng <- range(st_drop_geometry(foo2$layers[[1]][,1]), na.rm = T)
    }
    legend_cont("right", legend = rng, col = pal)
    suppressWarnings(par(pr))
  }
}
