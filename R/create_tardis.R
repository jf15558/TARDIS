#' create_tardis
#'
#' Generate a 3D lattice graph from a time-ordered set of geographic landscape
#' rasters. Weights for horizontal edges within layers record the great circle
#' distances between adjacent raster cells in metres, adjusted' for differences
#' in cell elevation using Pythagoras's theorem. On global scales this adjustment
#' will be a second order effect, but will be much more important at smaller
#' scales. Weights for edges between layers are assigned as zero so
#' that they do not affect downstream distance calculations. In many cases,
#' the positions of cells will remain constant within the extent of each layer
#' in geog and so the edges through time will essentially be vertical. For
#' global landscapes over geological timescales, however, the positions of
#' landmasses chance noticeably due to continental drift. In these cases, the
#' edges can be altered so that they connect pairs of geographically homologous
#' cells through time. Single layer cases (i.e. no time element) are allowed.
#'
#' @param geog A SpatRaster of digital elevation models recording a changing
#' landscape forwards in time (i.e. the first layer is the oldest). The
#' SpatRaster must be in geographic projection, i.e. longitude-latitude, with
#' elevations in metres. Negative values are permissible (e.g. bathymetry), but
#' missing values are not and should be replaced (e.g. with the average of the
#' surrounding cells).
#' @param times If there is only one layer in geog, then times is not used.
#' Otherwise, a vector of time bin boundaries with nlayers(geog) + 1 elements.
#' Each successive pair of entries in the vector will be treated as the time
#' range represented by successive layer in x. Entries should be positive numbers
#' expressing time in 'years ago'. The vector need not end in the present
#' (i.e. 0), but time must flow from oldest to youngest.
#' @param glink The linking mode for neighbouring grid cells within layers,
#' either 4 (Rook's case, orthogonal links only) or 8 (Queen's case, orthogonal
#' and diagonal).
#' @param tlink The linking mode for homologous grid cells between layers,
#' either 1 (forwards-in-time), 2 (backwards-in-time) or 3 (bidirectional).
#' The forwards-in-time case is the default, but there are cases where the
#' other options may be useful.
#' @param mask A SpatRaster of the same resolution, extent and number of layers
#' as geog, which will be used to designate non-traversable cells in the
#' corresponding layers in geog. This must contain only 1 (unmasked) or NA
#' (masked) values.
#' @param mask.check A logical determining whether mask should be checked for
#' islands using link_mask. This is the recommended default as bridges will be
#' added to ensure that the unmasked portion of the landscape is fully traversable.
#' @param klink The number of island connections to generate, as called by link_mask().
#' @param mlink If not NULL, an sf linestrings data.frame matching that produced
#' by link_mask(). This argument mainly exists to allow the user to supply
#' modified outputs of link_mask(). The validity of links are checked
#' internally, so modifications should be undertaken with care. If mlink is
#' supplied, then mask checking will not be performed, even if mask.check = TRUE
#' @param rotations A list with nlayers(geog) - 1 elements. Each element in the
#' list is a two-column numeric matrix containing the IDs of geographically
#' homologous cells between successive pairs of layers in x. Cell IDs are given
#' as R's standard cell number order for raster objects. Consequently, no cell
#' ID should be larger than ncell(geog). No list element can be blank as there
#' must be at least one edge between each pair of layers. Otherwise, any number
#' of edges between layers can be specified, although it  makes sense for there
#' to be maximally as many edges as cells within a layer.
#' @param verbose A logical to determine whether function progress should be
#' reported. Useful when dealing with large rasters (high resolution and/or many
#' layers)
#' @return A spatiotemporal graph object of class 'tardis'.
#' @import terra
#' @importFrom geosphere distGeo
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' #library(terra)
#' #library(TARDIS)
#' #gal <- galapagos()
#' #gal <- crop(gal, extent(-92, -88, -2, 1))
#' #gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' #gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)

create_tardis <- function(geog, times = NULL, glink = 8, tlink = 1,
                          mask = NULL, mask.check = TRUE, klink = NULL, mlink = NULL,
                          rotations = NULL, verbose = TRUE) {

  # geog = dem
  # times = bins
  # glink = 8
  # tlink = 1
  # mask = masks
  # mask.check = TRUE
  # klink = 1
  # mlink = NULL
  # rotations = NULL
  # verbose = TRUE

  # some precedent in Maartensen et al. (2017). Spatio-temporal connectivity: assessing the amount of reachable habitat in dynamic landscapes
  # useful for landscape graph theory: Decocq et al (2023). Modelling plant community dynamics in changing forest ecosystems: a review

  # check geog is correctly supplied
  if(!exists("geog")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!inherits(geog, "SpatRaster")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!is.lonlat(geog)) {
    stop("geog should be in geographic (long-lat) projection")
  }
  if(any(is.na(geog[]))) {
    stop("NA values present in geog")
  }

  if(!is.numeric(glink) | length(glink) != 1) {
    stop("glink should be one of 4 (Rook's case) or 8 (Queen's case)")
  }
  if(!glink %in% c(4, 8)) {
    stop("glink should be one of 4 or 8")
  }

  if(nlyr(geog) > 1) {

    if(!exists("times")) {
      stop("If there are multiple layers in geog, then times must be specified")
    }
    if(!is.numeric(times) | length(times) != nlyr(geog) + 1) {
      stop("Please supply times as a vector of time bin boundaries with nlayers(geog) + 1 elements")
    }
    if(any(diff(times) > 0)) {
      stop("All elements of times should be positive (i.e. before present) and in descending age order")
    }

    if(!is.numeric(tlink) | length(tlink) != 1) {
      stop("tlink should be one of 1 (forward-in-time), 2 (backward-in-time) or 3 (bidirectional")
    }
    if(!tlink %in% c(1:3)) {
      stop("tlink should be one of 1 (forward-in-time), 2 (backward-in-time) or 3 (bidirectional")
    }

    if(is.null(rotations)) {
      rotations <- list()
      for(i in 1:(nlyr(geog) - 1)) {
        rotations[[i]] <- cbind(from = 1:ncell(geog), to = 1:ncell(geog))
      }

    } else {

      if(!is.list(rotations)) {
        stop("Rotations should be a list")
      }
      if(!all(unlist(lapply(rotations, function(x) {class(x)[1]})) == "matrix")) {
        stop("All elements of rotations should be two-column matrices")
      }
      if(!all(unlist(lapply(rotations, function(x) {mode(x)[1]})) == "numeric")) {
        stop("All elements of rotations should be two-column matrices")
      }
      if(!all(unlist(lapply(rotations, ncol)) == 2)) {
        stop("All elements of rotations should be two-column matrices")
      }
      if(length(rotations) != (nlyr(geog) - 1)) {
        stop("Rotations should have nlayers(geog) - 1 elements")
      }
      if(any(unlist(lapply(rotations, max)) > ncell(geog[[1]]))) {
        stop("One of more of the cell IDs in rotations exceed the layer (ncell) extent of geog")
      }
    }
  }

  if(!is.logical(mask.check) | length(mask.check) != 1) {
    stop("mask.check should be a single logical")
  }
  if(!is.null(klink)) {
    if(length(klink) != 1 | !inherits(klink, "numeric")) {
      stop("If not NULL, klink should be an integer")
    }
    if(!klink %% 1 == 0) {
      stop("If not NULL, klink should be an integer")
    }
  }

  if(is.null(mask)) {

    if(!is.null(mlink)) {
      stop("If mlink is not NULL, then mask must also be supplied")
    }
    mask <- geog
    mask[] <- 1
    add_links <- lapply(1:nlyr(geog), function(x) {cbind(NA, NA)})

  } else {

    if(!inherits(mask, "SpatRaster")) {
      stop("If not NULL, supply mask as a SpatRaster")
    }
    if(any(!unique(mask[]) %in% c(1, NA))) {
      stop("Mask layers can only contain 1 or NA values")
    }
    if(!all(dim(geog) == dim(mask)) | ext(geog) != ext(mask)) {
      stop("geog and mask must all have the same extent, resolution, and number of layers")
    }
  }

  if(!is.null(mlink)) {

    warning("As mlink was supplied, masks were not checked inaccessible regions. TARDIS paths may fail unexpectedly")
    if(!inherits(mlink, "sf")) {
      stop("mlink should be an sf data.frame of linestrings")
    }
    if(!inherits(mlink, "data.frame")) {
      stop("mlink should be an sf data.frame of linestrings")
    }
    if(ncol(mlink) != 5) {
      stop("mlink should match the format of the output of link_mask()")
    }
    if(!all(colnames(mlink) == c("srt", "end", "bin", "distance", "geometry"))) {
      stop("mlink should match the format of the output of link_mask()")
    }
    if(!all(apply(st_drop_geometry(mlink), 2, is.numeric))) {
      stop("One or more columns in mlink is not numeric")
    }
    if(!all(apply(st_drop_geometry(mlink), 2, function(x) {any(!is.na(x))}))) {
      stop("One or more columns in mlink contains NA values")
    }
    if(any(mlink$bin < 1) | any(mlink$bin %% 1 != 0)) {
      stop("Only positive integers are permitted in mlink$bin")
    }
    if(nlyr(geog) == 1) {
      if(length(unique(mlink$bin)) > nlyr(geog)) {
        stop("mlink contains links for more time bins than layers present in geog")
      }
    } else {
      if(any(mlink$bin) > nlyr(geog)) {
        stop("mlink contains bin values exceeding the number of layers present in geog")
      }
    }
    if(!all(as.vector(st_geometry_type(mlink)) == "LINESTRING")) {
      stop("All geometries in mlink should be of type linestring")
    }
    if(!all(table(st_coordinates(mlink)[,3]) == 2)) {
      stop("Each linestring can only contain 2 coordinates (start and end)")
    }
    if(!all(abs(as.vector(ext(geog))) - abs(as.vector(ext(mlink))) >= 0)) {
      stop("The extent of mlink does not fall fully within the extent of mask")
    }

    add_links <- lapply(1:nlyr(geog), function(x) {
      if(x %in% mlink$bin) {
        vals <- extract(boundaries(mask[[x]], directions = 8), vect(mlink[which(mlink$bin == x),]))
        if(0 %in% vals[,2]) {
          stop(paste0("In bin ", x, ", one or more line start/end points do not fall on cells at the edges of islands"))
        }
        if(!all(table(vals[complete.cases(vals),1]) == 2)) {
          stop(paste0("In bin ", x, ", one or more lines intersect non-masked areas other than at their start and end points"))
        }
        crds <- matrix(cellFromXY(mask, st_coordinates(mlink[which(mlink$bin == x),])[,1:2]), ncol = 2, byrow = T)
        matrix(c(t(cbind(crds, crds[,2:1]))), ncol = 2, byrow = T)
      } else {
        cbind(NA, NA)
      }
    })

  } else {

    bar <- rast(lapply(as.list(mask), patches, allowGaps = F, directions = glink))
    if(any(minmax(bar)[2,] > 1) & mask.check == FALSE) {
      warning("mask contains inaccessible regions. TARDIS paths fail unexpectedly. Consider running with mask.check = TRUE")
    }
    if(mask.check) {
      lnk <- link_mask(mask, glink = glink, klink = klink, verbose = verbose)
      add_links <- lapply(1:nlyr(geog), function(x) {
        if(x %in% lnk$bin) {
          crds <- matrix(cellFromXY(mask, st_coordinates(lnk[which(lnk$bin == x),])[,1:2]), ncol = 2, byrow = T)
          matrix(c(t(cbind(crds, crds[,2:1]))), ncol = 2, byrow = T)
        } else {
          cbind(NA, NA)
        }
      })
    }
  }

  # spatial dimensions
  gdat <- c(rows = nrow(geog), cols = ncol(geog), hres = res(geog)[1], vres = res(geog)[2],
            xmin = ext(geog)[1], xmax = ext(geog)[2], ymin = ext(geog)[3], ymax = ext(geog)[4])

  # generate blank graph with full structure and uniform edge weights
  ed <- adjacent(geog, cells = 1:ncell(geog), directions = glink, pairs = T)

  # put into pairwise ordering (each pair of rows forms a directional edge)
  ed <- ed[ed[,1] < ed[,2],]
  ed <- matrix(c(t(cbind(ed, ed[,2:1]))), ncol = 2, byrow = T)

  # get horizontal distances (constant across layers)
  h_dists <- distGeo(xyFromCell(geog, ed[,1]),  xyFromCell(geog, ed[,2]))

  # get edges for each layer
  glinked <- list()
  for(i in 1:nlyr(geog)) {

    if(verbose) {
      cat(paste0("Linking geog [", i, "/", nlyr(geog), "]\r"))
      if(i == nlyr(geog)) {cat("\n")}
    }

    # temporary edgelist to include the mask links
    ed2 <- rbind(ed, na.omit(add_links[[i]]))
    h_dists2 <- c(h_dists, na.omit(distGeo(xyFromCell(geog, add_links[[i]][,1]),  xyFromCell(geog, add_links[[i]][,2]))))

    # get mask
    blocked <- which(is.na(mask[[i]][]))

    # get great circle distances adjusted for elevation by pythagoras
    v_dists <- geog[[i]][][ed2[,2]] - geog[[i]][][ed2[,1]]
    t_dists <- sqrt((h_dists2 ^ 2) + (abs(v_dists) ^ 2))

    # convert blocked cell weights to NA
    t_dists[which(ed2[,1] %in% blocked | ed2[,2] %in% blocked)] <- NA

    # store
    edge <- cbind(from = ed2[,1], to = ed2[,2], hdist = h_dists2, vdist = v_dists, tdist = t_dists)
    glinked[[i]] <- edge[complete.cases(edge),]
  }

  if(!is.null(times)) {
    for(i in 1:length(rotations)) {

      glinked[[i + 1]][,1:2] <- glinked[[i + 1]][,1:2] + (i * ncell(geog))

      ob <- rotations[[i]]
      ob[,1] <- ob[,1] + ((i - 1) * ncell(geog))
      ob[,2] <- ob[,2] + (i * ncell(geog))
      ob[!ob[,1] %in% glinked[[i]][,1], 1] <- NA
      ob[!ob[,2] %in% glinked[[i + 1]], 2] <- NA
      ob <- ob[complete.cases(ob), , drop = F]
      ob <- ob[!duplicated(ob), , drop = F]

      if(nrow(ob) == 0) {
        stop(paste0("No links are available from time layer ", i, ". Check rotations and layer masks"))
      }
      if(tlink == 2) {ob <- ob[,2:1]}
      if(tlink == 3) {ob <- matrix(c(t(cbind(ob, ob[,2:1]))), ncol = 2, byrow = T)}

      glinked[[i]] <- rbind(glinked[[i]], cbind(from = ob[,1], to = ob[,2], hdist = 0, vdist = 0, tdist = 0))
    }
  }

  # create modified cppRouting graph object and return
  if(verbose) {cat("Building graph\n")}

  glinked <- do.call(rbind, glinked)
  src <- as.character(glinked[,1])
  dst <- as.character(glinked[,2])
  nodes <- unique(c(src, dst))
  id = 0:(length(nodes) - 1)
  out <- list(edges = glinked, gdat = gdat, tdat = times, link.mode = c(glink = glink, tlink = tlink),
              tgraph = list(data = NULL, dict = data.frame(ref = nodes, id = id), coords = NULL, nbnode = length(nodes),
                            attrib = list(aux = NULL, cap = NULL, alpha = NULL, beta = NULL), src = src, dst = dst))

  class(out) <- "tardis"
  return(out)
}
