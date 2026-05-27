#' plot.geoglist
#'
#' Plotting method for a geoglist layer. If the geoglist contains multiple
#' layers, then default behaviour is to plot the first one.
#'
#' @name plot
#' @method plot geoglist
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
#' @param axes `logical`. Should axes be added to the plot? These will look
#' sensible for lon-lat geoglists, but may look odd for other projection systems.
#' @param bg `character` or `integer`. The colour to use for the map background.
#' White by default.
#' @return None.
#' @param xlim `numeric`. If not `NULL`, then a vector of two numbers to set the
#' minimum and maximum x extent of the plot in terms of the projection system in
#' `geoglist`.
#' @param ylim `numeric`. If not `NULL`, then a vector of two numbers to set the
#' minimum and maximum y extent of the plot in terms of the projection system in
#' `geoglist`.
#' @param zlim `numeric`. If not `NULL`, then a vector of two numbers to set the
#' range on the plotting legend.
#' the range of values
#' @import sf terra
#' @importFrom graphics par
#' @importFrom graphics axis
#' @importFrom graphics axTicks
#' @export
#'
#' @examples
#' \dontrun{
#' gal <- cretaceous()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts, klink = 1)
#'
#' plot(regs)
#' }

plot.geoglist <- function(geog, layer = 1, pal = sf.colors(10), links = T,
                          lcol = "grey", lwd = 1, lty = 1, hex.border = NA,
                          legend = T, axes = T, bg = "white",
                          xlim = NULL, ylim = NULL, zlim = NULL) {

   # geog = rasts2
   # layer = 1
   # pal = sf.colors(10)
   # links = T
   # lcol = 1
   # lwd = 1
   # lty = 1
   # hex.border = NA
   # legend = T
   # axes = T
   # xlim = NULL
   # ylim = NULL
   # #xlim = c(-91, -89.5)
   # #ylim = c(-1.5, 0)
   # #xlim = c(-91, 50)
   # #ylim = c(-50, 80)
   # zlim = NULL
   # xlim = c(-1.5e7, 0)
   # ylim = c(-5e6, 1e7)

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
  if(!is.null(xlim)) {
    if(!is.atomic(xlim) | length(xlim) != 2) {
      stop("xlim should be a vector of length 2")
    }
    if(!is.numeric(xlim)) {
      stop("xlim should be numeric")
    }
  }
  if(!is.null(ylim)) {
    if(!is.atomic(ylim) | length(ylim) != 2) {
      stop("ylim should be a vector of length 2")
    }
    if(!is.numeric(ylim)) {
      stop("ylim should be numeric")
    }
  }
  if(!is.null(zlim)) {
    if(!is.atomic(zlim) | length(zlim) != 2) {
      stop("zlim should be a vector of length 2")
    }
    if(!is.numeric(zlim)) {
      stop("zlim should be numeric")
    }
  }

  if(legend) {
    pr <- newmar <- par("mar")
    newmar[4] <- 5.1
    par(mar = newmar)
  }

  # construct bounding polygon (adjust for global extent to avoid errors)
  bbox <- geog$gdat[1:4]
  if(bbox[1] < -179.99) {bbox[1] <- -179.99}
  if(bbox[2] > 179.99) {bbox[2] <- 179.99}
  if(bbox[3] < -89.99) {bbox[3] <- -89.99}
  if(bbox[4] > 89.99) {bbox[4] <- 89.99}
  len1 <- diff(bbox[1:2]) / 0.1
  len2 <- diff(bbox[3:4]) / 0.1
  frame <- cbind(
    c(seq(bbox[1], bbox[2], length.out = len1), rep(bbox[2], len2), seq(bbox[2], bbox[1], length.out = len1),  rep(bbox[1], len2)),
    c(rep(bbox[4], len1), seq(bbox[4], bbox[3], length.out = len2), rep(bbox[3], len1), seq(bbox[3], bbox[4], length.out = len2))
  )
  frame <- st_make_valid(st_sfc(st_polygon(list(frame)), crs = "EPSG:4326"))
  frame <- st_transform(frame, crs(geog$layers[[1]]))
  st_crs(frame) <- NA

  # crop bounding box to plot limits and plot
  if(is.null(xlim)) {xlim <- st_bbox(frame)[c(1, 3)]}
  if(is.null(ylim)) {ylim <- st_bbox(frame)[c(2, 4)]}
  names(xlim) <- c("xmin", "xmax")
  names(ylim) <- c("ymin", "ymax")

  # set plot canvas
  pol <- st_sfc(st_polygon(list(cbind(xlim[c(1, 1, 2, 2, 1)], ylim[c(1, 2, 2, 1, 1)]))))
  bounds <- st_intersection(frame, pol)
  plot(bounds, col = bg, border = NA)

  # add geoglist layers
  if(inherits(geog$layers[[1]], "SpatRaster")) {

    # mask and crop to plotting bounds
    rst <- suppressWarnings(mask(crop(geog$layers[[layer]], vect(bounds)), vect(bounds)))
    # adjust the boundary polygon so that it conforms to the raster grid resolution
    frame <- st_crop(frame, y = as.vector(ext(rst)))
    # plot
    plot(rst, col = pal, legend = F, add = T)

  } else {
    # crop to plotting bounds
    lyr <- geog$layers[[layer]]
    st_crs(lyr) <- NA
    lyr <- suppressWarnings(st_intersection(lyr, bounds))
    # adjust the boundary polygon frame to plot bounds
    frame <- st_crop(frame, st_bbox(bounds))
    # plot
    plot(lyr[,1], add = T, pal = pal, border = hex.border)
  }

  # add axes
  if(axes) {
    a1 <- axTicks(1)
    axis(1, pos = st_bbox(frame)[2], at = a1[which(a1 >= st_bbox(frame)[1] & a1 <= st_bbox(frame)[3])],
         cex.axis = 0.7, col = "grey", padj = -1)
    a2 <- axTicks(2)
    axis(2, pos = st_bbox(frame)[1], at = a2[which(a2 >= st_bbox(frame)[2] & a2 <= st_bbox(frame)[4])],
         cex.axis = 0.7, col = "grey", padj = 0.8)
  }

  # add links
  if(links) {
    if(!is.null(geog$links)) {
      lnk <- st_crop(geog$links[which(geog$links$layer == layer),"geometry"], st_bbox(frame))
      plot(lnk, add = T, col = lcol, lwd = lwd, lty = lty)
    }
  }

  # add legend
  if(legend) {
    if(is.null(zlim)) {
      if(inherits(geog$layers[[1]], "SpatRaster")) {
        zlim <- c(minmax(geog$layers[[layer]]))
      } else {
        zlim <- range(st_drop_geometry(geog$layers[[1]][,1]), na.rm = T)
      }
    }
    legend_cont(x = seq(st_bbox(frame)[3], par("usr")[2], length.out = 4)[2:3],
                y = st_bbox(bounds)[c(2, 4)], legend = zlim, col = pal)
    suppressWarnings(par(mar = pr))
  }

  # add plot boundary
  plot(frame, add = T)
}

