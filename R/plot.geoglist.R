#' plot.geoglist
#'
#' Plotting method for a geoglist layer. If the geoglist contains multiple
#' layers, then default behaviour is to plot the first one.
#'
#' @name plot
#' @method plot geoglist
#' @param x `geoglist`. The output of `rast_to_geoglist()`.
#' @param y `numeric`. The layer in the geoglist to be plotted, along with
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
#' `NA` (no colour) by default.
#' @param add `logical`. Should the geoglist layer be added to an existing plot?
#' `FALSE` by default.
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
#' @param ... Other arguments passed to `terra::plot()`. Use will probably cause
#' errors.
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

plot.geoglist <- function(x, y = 1, pal = sf.colors(10), links = T,
                          lcol = "grey", lwd = 1, lty = 1, hex.border = NA,
                          legend = T, axes = T, bg = NA, add = F,
                          xlim = NULL, ylim = NULL, zlim = NULL, ...) {

   # x = rasts
   # y = 1
   # pal = sf.colors(10)
   # links = T
   # lcol = 1
   # lwd = 1
   # lty = 1
   # hex.border = NA
   # legend = T
   # axes = T
   # add = F
   # bg = NA
   # xlim = NULL
   # ylim = NULL
   # #xlim = c(-91, -89.5)
   # #ylim = c(-1.5, 0)
   # #xlim = c(-91, 50)
   # #ylim = c(-50, 80)
   # zlim = NULL
   # #xlim = c(-1.5e7, 0)
   # #ylim = c(-5e6, 1e7)

  if(!exists("x")) {
    stop("Supply x as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(x, "geoglist")) {
    stop("Supply x as a geoglist from rast_to_geoglist()")
  }

  if(!is.atomic(y) | length(y) != 1) {
    stop("y should be a single integer")
  }
  if(!is.numeric(y)) {
    stop("y should be a single integer")
  }
  if(y %% 1 != 0) {
    stop("y should be a single integer")
  }
  if(y > ifelse(inherits(x$layers, "SpatRaster"), nlyr(x$layers), length(x$layers))) {
    stop("The value of y exceeds of the number of layers in x")
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
  bbox <- x$gdat[1:4]
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
  frame <- as.polygons(vect(frame, crs = "EPSG:4326"))
  frame <- project(frame, crs(x$layers))
  crs(frame) <- NA

  # crop bounding box to plot limits and plot
  if(is.null(xlim)) {xlim <- ext(frame)[1:2]}
  if(is.null(ylim)) {ylim <- ext(frame)[3:4]}
  names(xlim) <- c("xmin", "xmax")
  names(ylim) <- c("ymin", "ymax")

  # set plot canvas
  pol <- as.polygons(vect(cbind(xlim[c(1, 1, 2, 2, 1)], ylim[c(1, 2, 2, 1, 1)])))
  frame <- crop(frame, pol)
  bounds <- intersect(frame, pol)

  # mask and crop to plotting bounds
  rst <- x$layers[[y]]
  crs(rst) <- NA
  rst <- mask(crop(rst, bounds), bounds)

  # if the frame bisects a raster row or column, adjust so it falls in step
  if(inherits(rst, "SpatRaster")) {
    #foo <- ext(frame)
    #foo2 <- ext(rst)
    #frame <- crop(frame, bounds, ext = T)
  }

  plot(frame, col = bg, border = NA, axes = F, add = add)
  if(inherits(rst, "SpatRaster")) {
    plot(rst, col = pal, xlim = ext(frame)[1:2], ylim = ext(frame)[3:4],
         axes = F, add = T, legend = F, range = zlim, ...)
  } else {
    plot(rst, values = rst[[1]][,1], col = pal, xlim = ext(frame)[1:2], ylim = ext(frame)[3:4],
         border = hex.border, axes = F, add = T, legend = F, range = zlim, ...)
  }

  # add axes
  if(axes) {
    a1 <- axTicks(1)
    axis(1, pos = ext(frame)[3], at = a1[which(a1 >= ext(frame)[1] & a1 <= ext(frame)[2])],
         cex.axis = 0.7, col = "grey", padj = -1)
    a2 <- axTicks(2)
    axis(2, pos = ext(frame)[1], at = a2[which(a2 >= ext(frame)[3] & a2 <= ext(frame)[4])],
         cex.axis = 0.7, col = "grey", padj = 0.8)
  }

  # add links
  if(links) {
    if(!is.null(x$links)) {
      lnk <- x$links[which(x$links$layer == y)]
      crs(lnk) <- NA
      lnk <- intersect(lnk, frame)
      plot(lnk, add = T, col = lcol, lwd = lwd, lty = lty)
    }
  }

  # add legend
  if(legend) {
    if(is.null(zlim)) {
      if(inherits(x$layers[[1]], "SpatRaster")) {
        zlim <- c(minmax(x$layers[[y]]))
      } else {
        zlim <- range(rst, na.rm = T)
      }
    }
    legend_cont(x = seq(ext(frame)[2], par("usr")[2], length.out = 4)[2:3],
                y = ext(bounds)[3:4], legend = zlim, col = pal)
    suppressWarnings(par(mar = pr))
  }

  # add plot boundary
  plot(frame, add = T)
}

