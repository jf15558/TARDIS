#' commute
#'
#' Calculate the commute time between the start and end of a least cost path.
#' The commute time is the degree of 'resistance' between two points in a
#' graph and provides a measure of the practical difficulty of traversing
#' between them, regardless of the optimal route offered by the intervening
#' least cost path. The commute time calculation requires a graph with
#' bidirectional edges. Consequently, the average of the weights for pairs of
#' edges are calculated within each graph layer, rather than considering their
#' potential asymmetries.
#'
#' @param tardis An object of class 'tardis', produced by create_tardis
#' @param weights If not NULL, a vector of weights to be used instead of the
#' geographic distances in tardis. All entries must be >= 0
#' and finite (NaN or Inf), or NA. Typically the output of weight_tgraph
#' @param origin A simple features collection produced by stp, denoting
#' the origin cells for the commutes
#' @param dest As for origin, but for the destination points. Unlike for lcp
#' the directionality in time does not matter
#' @param verbose A logical indicating whether function progress should be
#' reported to the user.
#' @return A list of vectors containing the layer-discrete commute times for
#' each segment of the input least cost paths
#' @import terra Matrix
#' @export

commute <- function(tardis, weights = NULL, origin, dest, verbose = TRUE) {

  #tardis = test2
  #weights = test_weight
  #origin = pts[1:2,]
  #dest = pts[3:4,]
  #verbose = T

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!class(tardis) == "tardis") {
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

  if(!class(origin)[1] == c("sf")) {
    stop("Supply origin as the output of stp")
  }
  if(!class(dest)[1] == c("sf")) {
    stop("Supply dest as the output of stp")
  }
  if(nrow(origin) != nrow(dest)) {
    stop("The number of origin and destination points should be the same (i.e. paired points")
  }

  # get cell id from geographic position and age
  srt <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(origin)) +
    (prod(tardis$gdat[1:2]) * (origin$bin - 1))

  end <- cellFromXY(rast(nrows = tardis$gdat[1], ncols = tardis$gdat[2], ext = ext(tardis$gdat[5:8])), st_coordinates(dest)) +
    (prod(tardis$gdat[1:2]) * (dest$bin - 1))

  # check point accessibility to ensure they come from the correct tardis object
  if(!all(c(srt, end) %in% tardis$edges[,1])) {
    stop("One or more points in paths do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for paths")
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

  # add bidirectional time edges if needed
  if(tardis$link.mode[2] != 3) {

    gl <- which(tardis$edges[,5] != 0)
    tl <- which(tardis$edges[,5] == 0)
    newtl <- c(rbind(tardis$tgraph$src[tl], tardis$tgraph$dst[tl], tardis$tgraph$dst[tl], tardis$tgraph$src[tl]))
    tardis$tgraph$src <- c(tardis$tgraph$src[gl], newtl[seq(1, length(newtl), 2)])
    tardis$tgraph$dst <- c(tardis$tgraph$dst[gl], newtl[seq(2, length(newtl), 2)])
    wts <- c(tardis$edges[gl,5], rep(0, length(tl) * 2))

  } else {
    wts <- tardis$edges[,5]
  }

  # get sparse matrix of pairwise edge weight averages (constant added to avoid the zero weights)
  mns <- as.vector(rep(tapply(wts, rep(1:(length(wts) / 2), each = 2), mean), each = 2)) + 1
  src <- tardis$tgraph$dict$id[match(tardis$tgraph$src, tardis$tgraph$dict$ref)] + 1
  dst <- tardis$tgraph$dict$id[match(tardis$tgraph$dst, tardis$tgraph$dict$ref)] + 1
  adj_mat <- sparseMatrix(i = src, j = dst, x = mns)

  # commute calculation
  srt <- tardis$tgraph$dict$id[match(srt, tardis$tgraph$dict$ref)] + 1
  end <- tardis$tgraph$dict$id[match(end, tardis$tgraph$dict$ref)] + 1
  allcls <- c(srt, end)
  Lr <- Matrix::Diagonal(x = colSums(adj_mat)) - adj_mat
  n <- max(Lr@Dim)
  Lr <- Lr[-n, -n]
  C <- 1e-300 * n
  Lplus <- matrix(ncol = length(allcls), nrow = length(allcls))

  for(i in 1:length(allcls)) {
    ei <- matrix(-C/n, ncol = 1, nrow = n - 1)
    ei[allcls[i], ] <- C - (C/n)
    xi <- solve(Lr, ei)
    xi <- as.vector(xi)
    Lplusallrows <- c(xi - sum(xi/n), (sum(xi)/n))
    Lplus[, i] <- as.vector(Lplusallrows)[allcls]
  }

  Lplus <- Lplus/C
  m1 <- matrix(diag(Lplus), nrow = length(allcls), ncol = length(allcls))
  m2 <- t(matrix(diag(Lplus), nrow = length(allcls), ncol = length(allcls)))
  vol <- sum(adj_mat)
  rdSS <- (-2 * Lplus + m1 + m2) * sum(adj_mat)

  # return
  return(data.frame(path = 1:length(srt), commute = rdSS[cbind(match(srt, allcls), match(end, allcls))]))
}
