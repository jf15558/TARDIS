#' detours
#'
#' Calculate detours around least cost paths from tardis. The function identifies
#' all cells accessible from the least cost path within a certain proportional
#' additional cost relative to the cost of that least cost path. Detour
#' calculation within the tardis spatiotemporal graphs is performed by the
#' get_detour function of the cppRouting package. As defined there, accessible
#' cells (n) from the least cost path are those that meet the condition:
#' SP(o,n) + SP(n,d) < SP(o,d) + t, where SP is the cost of access, o the origin
#' node, d the destination node and t the extra cost. Any identified cells are
#' then converted to coordinates and the concave hull for that set of coordinates
#' calculated and returned.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in tardis. All entries must be >= 0
#' and finite (NaN or Inf), or NA. Typically the output of weight_tgraph
#' @param paths An object returned by the lcp function.
#' @param pindex The indices of the least cost paths in 'paths' to calculate
#' detours for. If NULL, then all paths will be analysed.
#' @param detour The proportional extra cost to detect accessible cells within,
#' a numeric > 0. The default is 1% (i.e. 0.01). Optionally, a vector with as
#' many elements as pindex, to enable different detour amounts for each path.
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return If no detours could be calculated for any of the paths, the function
#' will return NULL. Otherwise a simple features geometry collection, with a
#' multipolygon for each path where a detour could be calculated. As least cost
#' paths can extend across multiple time bins, the polygons in a multipolygon
#' represent bin-specific segments of the overall detour.
#' @import terra sf concaveman
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#'
#' #gal <- galapagos()
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1),
#' #                                    ncol = 3, byrow = T), right = F)
#' #gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
#'
#' #vars = list(elev = classify(gal, cbind(-Inf, 0, 0)))
#' #gtw <- weight_tardis(test2, vars = vars,
#' #                     mfun = function(origin, dest, lnum, ...) {
#' #                               (origin$hdist^2 + abs(origin$vdist)^2) * 10})
#'
#' #org <- rbind(c(-89, -1.05, 2), c(-89.5, -0.7, 2))
#' #dst <- rbind(c(-91.2, -1, 0), c(-91.6, -0.4, 0))
#' #pts <- stp(test2, rbind(org, dst))
#'
#' #foo <- lcp(tardis = gt, weights = gtw, pts[1:2,], pts[3:4,])
#' #foo2 <- detour(tardis = gt, weights = gtw, foo, detour = 0.01)

detour <- function(tardis, weights = NULL, paths, pindex = NULL, detour = 0.01, verbose = TRUE) {

  # tardis = test2
  # weights = test_weight
  # paths = test3
  # detour = 0.01
  # verbose = T
  # pindex = NULL

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }

  if(!is.null(weights)) {

    if(!is.atomic(weights)) {
      stop("weights must be a vector")
    }
    if(!is.numeric(weights)) {
      stop("weights (the vector or chosen column) must be numeric")
    }
    if(any(is.infinite(weights) | weights < 0 | is.nan(weights))) {
      stop("All entries in weights must be positive, finite numerics")
    }
    if(length(weights) != nrow(tardis$edges)) {
      stop("weights must be the same length as the number of edges in tardis")
    }
  } else {
    weights <- tardis$edges[,5]
  }

  if(!class(paths)[1] == "sf") {
    stop("Supply paths as the output of lcp()")
  }
  if(!all(c("path", "bin", "order", "cost", "distance", "geometry") %in% names(paths))) {
    stop("Supply paths as the output of lcp()")
  }

  if(!is.null(pindex)) {

    if(!is.atomic(pindex) | !is.numeric(pindex)) {
      stop("If not NULL, pindex should be a vector of indexes of the paths for which detours will be calculated")
    }
    if(any(pindex < 1 | pindex > max(paths$path))) {
      stop("One or more indices in pindex falls outside the range of paths")
    }
    pindex <- pindex[order(pindex)]

  } else {
    pindex <- 1:max(paths$path)
  }

  if(!is.atomic(detour) | !is.numeric(detour)) {
    stop("If not NULL, detour should be a single positive, finite numeric, or a vector of the same with as many elements as pindex")
  }
  if(length(detour) == 1) {
    detour <- rep(detour, length(pindex))
  }
  if(length(detour) != length(pindex)) {
    stop("If not NULL, detour should be a single positive, finite numeric, or a vector of the same with as many elements as pindex")
  }
  if(!all(is.finite(detour) | detour <= 0)) {
    stop("If not NULL, detour should be a single positive, finite numeric, or a vector of the same with as many elements as pindex")
  }

  # obtain start and end points, and costs
  paths <- paths[order(paths$path, paths$order, method = "radix"),]
  cls <- tapply(paths$geometry, paths$path, function(x) {ob <- st_coordinates(x); ob[c(1, nrow(ob)),1:2]})
  org <- unlist(lapply(cls, function(x) {cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), x[1,,drop = F])}))
  dst <- unlist(lapply(cls, function(x) {cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), x[2,,drop = F])}))
  org <- data.frame(cell = org, bin = tapply(paths$bin, paths$path, function(x) {x[1]}))
  dst <- data.frame(cell = dst, bin = tapply(paths$bin, paths$path, function(x) {rev(x)[1]}))
  org$cell <- org$cell + (prod(tardis$gdat[1:2]) * (org$bin - 1))
  dst$cell <- dst$cell + (prod(tardis$gdat[1:2]) * (dst$bin - 1))

  # check point accessibility to ensure they come from the correct tardis object
  if(!all(c(org$cell %in% tardis$edges[,1], dst$cell %in% tardis$edges[,2]))) {
    stop("One or more points in paths do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for paths")
  }
  if(any(dst$bin > org$bin) & tardis$link.mode[2] == 2) {
    stop("tardis is linked backwards in time, but some paths terminate in younger layers than their starting point")
  }
  if(any(dst$bin < org$bin) & tardis$link.mode[2] == 1) {
    stop("tardis is linked forwards in time, but some paths terminate in older layers than their starting point")
  }

  # initialise graph
  if(any(is.na(weights))) {
    tardis$tgraph$src <- tardis$tgraph$src[which(!is.na(weights))]
    tardis$tgraph$dst <- tardis$tgraph$dst[which(!is.na(weights))]
    tardis$edges <- tardis$edges[which(!is.na(weights)),]
    tardis$tgraph$dict <-  tardis$tgraph$dict[which(tardis$tgraph$dict$ref %in% c(tardis$tgraph$src, tardis$tgraph$dst))]
    tardis$tgraph$id <- 0:(length(tardis$tgraph$dict$ref) - 1)
    tardis$tgraph$nbnode <- length(tardis$tgraph$id)
  }
  tardis$tgraph$attrib$aux <- tardis$edges[!is.na(weights),5]
  tardis$edges[,5] <- weights[!is.na(weights)]
  tardis$tgraph$data <- data.frame(from = tardis$tgraph$dict$id[match(tardis$tgraph$src, tardis$tgraph$dict$ref)],
                                   to   = tardis$tgraph$dict$id[match(tardis$tgraph$dst, tardis$tgraph$dict$ref)],
                                   dist = tardis$edges[,5])
  tardis$tgraph <- tardis$tgraph[1:5]

  samprast <- rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8]))
  pweight <- tapply(paths$distance, INDEX = paths$path, sum)
  if(!is.null(weights)) {
    tardis$tgraph$data$dist <- weights
    pweight <- tapply(paths$cost, INDEX = paths$path, sum)
  }

  det_list <- list()
  foo <- lapply(pindex, function(i) {

    if(verbose) {cat(paste0("Running detours [", i, "/", length(pindex), "]\r"))
      if(i == length(pindex)) {cat("\n")}
    }

    ob <- NULL
    if(pweight[i] > 0) {
      iso_cell <- as.numeric(suppressMessages(get_detour(tardis$tgraph, from = org$cell[i], to = dst$cell[i], extra = pweight[i] * detour[i]))[[1]])

      if(length(iso_cell) > 0) {

        iso_bin <- ceiling(iso_cell / prod(tardis$gdat[1:2]))
        iso_seq <- rep(1:length(rle(iso_bin)$length), rle(iso_bin)$length)
        iso_ids <- paste0(i, "_", rle(iso_bin)$values, "-", rle(iso_bin)$values)

        iso_cell <- iso_cell - ((iso_bin - 1) * prod(tardis$gdat[1:2]))
        ob <- data.frame(path = rep(i, length(rle(iso_bin)$values)), bin = rle(iso_bin)$values, order = 1:length(iso_ids))
        rownames(ob) <- iso_ids

        isochrone <- tapply(iso_cell, INDEX = iso_seq, function(y) {
          tmprast <- samprast
          tmprast[y] <- 1
          tmprast <- patches(tmprast)
          st_multipolygon(tapply(which(!is.na(tmprast[])), INDEX = tmprast[which(!is.na(tmprast[]))], function(z) {
            iso_crd <- xyFromCell(tmprast, z)
            pts <- cbind(c(iso_crd[,1] - res(tmprast)[1], iso_crd[,1] - res(tmprast)[1],
                           iso_crd[,1] + res(tmprast)[1], iso_crd[,1] + res(tmprast)[1]),
                         c(iso_crd[,2] - res(tmprast)[2], iso_crd[,2] + res(tmprast)[2],
                           iso_crd[,2] + res(tmprast)[2], iso_crd[,2] - res(tmprast)[2]))
            st_polygon(list(concaveman(pts)))
          }, simplify = F))
        }, simplify = F)
        st_geometry(ob) <- st_sfc(isochrone)
      }
    }
    ob
  })

  # finish up
  if(all(unlist(lapply(foo, is.null)))) {
    if(verbose) {warning("No detours were present for any of the paths at the supplied threshold(s), returning NULL")}
    foo <- NULL
  } else {
    foo <- do.call(rbind, foo)
  }
  return(foo)
}
