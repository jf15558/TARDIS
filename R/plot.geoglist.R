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
                          lcol = 1, lwd = 1, lty = 1, hex.border = NA, legend = T,
                          axes = T) {

  # geog = rasts2
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
    pr <- newmar <- par("mar")
    newmar[4] <- 5.1
    par(mar = newmar)
  }

  # construct bounding polygon in lon lat and transform if needed
  bbox <- geog$gdat[1:4]
  len1 <- diff(bbox[1:2]) / 0.01
  len2 <- diff(bbox[3:4]) / 0.01
  bounds <- cbind(
    c(seq(bbox[1], bbox[2], length.out = len1), rep(bbox[2], len2), seq(bbox[2], bbox[1], length.out = len1),  rep(bbox[1], len2)),
    c(rep(bbox[4], len1), seq(bbox[4], bbox[3], length.out = len2), rep(bbox[3], len1), seq(bbox[3], bbox[4], length.out = len2))
  )
  bounds <- st_sfc(st_polygon(list(bounds)), crs = "EPSG:4326")
  bounds <- st_transform(bounds, crs(geog$layers[[1]]))

  #plot(0, 0, type = "n", xlim = geog$gdat[1:2], ylim = geog$gdat[3:4],
  #     asp = 1, xaxs = "i", yaxs = "i", xlab = "", ylab = "", cex.axis = 0.7, axes = F)
  plot(bounds)
  if(inherits(geog$layers[[1]], "SpatRaster")) {
    plot(geog$layers[[layer]], col = pal, legend = F, add = T)
  } else {
    plot(geog$layers[[layer]][,1], add = T, pal = pal, border = hex.border)
  }
  if(axes) {
    axis(1, pos = geog$gdat[3], cex.axis = 0.7, col = "grey80", padj = -1)
    a2 <- axTicks(2)
    axis(2, pos = geog$gdat[1], at = a2[which(a2 >= geog$gdat[3] & a2 <= geog$gdat[4])],
         cex.axis = 0.7, col = "grey80", padj = 0.8)
  }
  rect(geog$gdat[1], geog$gdat[3], geog$gdat[2], geog$gdat[4])

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
      rng <- range(st_drop_geometry(geog$layers[[1]][,1]), na.rm = T)
    }
    legend_cont("right", legend = rng, col = pal)
    suppressWarnings(par(mar = pr))
  }
}
