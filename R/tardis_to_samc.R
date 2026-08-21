#' tardis_to_samc
#'
#' Create a spatial absorbing markov chain (samc) matrix from a rTARDIS graph
#' layer. This object can then be analysed using functions from the samc R
#' package.
#'
#' @param tardis `tardis` A tardis object
#' @param weights `character`. A character string denoting the weighting scheme to use. By
#' default these are true geographic distances (gdist). Alternatively, the name
#' of a weighting scheme added to the tardis object with weight_tardis().
#' @param layer `numeric`. The layer number to convert to a samc matrix.
#' @param time `numeric`. Instead of a layer number, a numeric age can be
#' supplied. This is resolved to a layer number internally.
#' @param absorption `geoglist` or `numeric`. A `geoglist` with layers of absorption
#' probabilities corresponding to the layers in `tardis`. Otherwise a constant
#' probability can be used (0 by default). All probabilities must be below 1.
#' @param fidelity `geoglist` or `numeric`. A `geoglist` with layers of fidelity
#' probabilities corresponding to the layers in `tardis`. Otherwise a constant
#' probability can be used (0 by default). Probabilities equal to 1 are permitted.
#' @return An object of class `samc`.
#' @import Matrix
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph components
#' @importFrom samc samc
#'
#' @details
#' Internally, tardis_to_sparse() is used to generate the Q matrix component
#' of the samc. In this initial Q matrix, each row corresponds to an accessible
#' cell in tardis, with probabilities calculated as the proportion of the total
#' cost of movements from that cell so that they sum to 1. The absorption values
#' are then subtracted proportionally from each row-wise set of transition
#' probabilities so that transition probabilities retain their proportional
#' relationships with one another, and so that the transition and absorption
#' probabilities for that row sum to one.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(rTARDIS)
#'
#' # load galapagos dataset
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # hex-gridded tardis
#' rasts <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = T, hex = 7)
#' rasts <- link_islands(rasts, klink = NULL)
#' rtd <- build_tardis(rasts)
#'
#' # with constant
#' smc <- tardis_to_samc(tardis = rtd, absorption = 0.001)
#'
#' # with geoglist
#' absr <- rasts
#' not_missing <- which(!is.na(absr$layers[[1]]$layer))
#' absr$layers[[1]][[1]][not_missing,] <- runif(length(not_missing), 0, 0.999)
#' smc <- tardis_to_samc(tardis = rtd, absorption = absr)
#' }

tardis_to_samc <- function(tardis, weights = "gdist", layer = 1, time = NULL,
                           absorption = 0, fidelity = 0) {

  # tardis = rtd1
  # weights = "gdist"
  # absorption = foo1
  # fidelity = 0

  if(!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
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

  if(is.null(tardis$tdat)) {
    tardis$tdat <- c(1.1, 0.9)
    time <- NULL
  }
  if(!is.null(time)) {
    if(!is.atomic(time)) {
      stop("time should be a singe positive number")
    }
    if(!is.numeric(time) | length(time) != 1) {
      stop("time should be a singe positive number")
    }
    if(time < 0 | time > max(tardis$tdat)) {
      stop("time should be a singe positive number in the time range of tardis")
    }
    layer <- sum(which(time > tardis$tdat))
  }
  if(!is.atomic(layer)) {
    stop("layer should be a single positive integer")
  }
  if(!is.numeric(layer) | length(layer) != 1) {
    stop("layer should be a single positive integer")
  }
  if(layer %% 1 != 0 | layer < 0 | layer > (length(tardis$tdat) - 1)) {
    stop("layer should be a single positive integer in the layer range of tardis")
  }

  if(!inherits(absorption, "geoglist") & !is.atomic(absorption)) {
    stop("Absorption should be a geoglist with layers corresponding to tardis, or a positive numeric value")
  }
  if(inherits(absorption, "geoglist")) {
    if(absorption$gdat != tardis$gdat) {
      stop("absorption and tardis spatial parameters do not match")
    }
    ly <- ifelse(inherits(absorption$layers, "SpatRaster"), nlyr(absorption$layers), length(absorption$layers))
    if(ly != (length(tardis$tdat) - 1)) {
      stop("If absorption is a geoglist, its layers must correspond to those in tardis")
    }
  } else {
    if(length(absorption) != 1 | !is.numeric(absorption)) {
      stop("If not a geoglist, absorption should be a single positive numeric")
    }
  }

  if(!inherits(fidelity, "geoglist") & !is.atomic(fidelity)) {
    stop("Fidelity should be a geoglist with layers corresponding to tardis, or a positive numeric value")
  }
  if(inherits(fidelity, "geoglist")) {
    if(fidelity$gdat != tardis$gdat) {
      stop("fidelity and tardis spatial parameters do not match")
    }
    ly <- ifelse(inherits(fidelity$layers, "SpatRaster"), nlyr(fidelity$layers), length(fidelity$layers))
    if(ly != (length(tardis$tdat) - 1)) {
      stop("If fidelity is a geoglist, its layers must correspond to those in tardis")
    }
  } else {
    if(length(fidelity) != 1 | !is.numeric(fidelity)) {
      stop("If not a geoglist, fidelity should be a single positive numeric value between 0 and 1")
    }
  }

  if(length(tardis$tdat) > 2) {
    tardis <- slice_tardis(tardis, layers = c(layer, layer))
    if(inherits(absorption, "geoglist")) {
      absorption <- slice_geoglist(absorption, layers = layer)
    }
    if(inherits(fidelity, "geoglist")) {
      fidelity <- slice_geoglist(fidelity, layers = layer)
    }
  }

  if(inherits(absorption, "geoglist")) {
    absorption <- na.omit(absorption$layers[[1]], field = names(absorption$layers[[1]])[1])
    absorption <- absorption[[1]][,1]
  } else {
    absorption <- rep(absorption, length(unique(tardis$tgraph$src)))
  }
  if(any(absorption >= 1 | absorption < 0)) {
    stop("Absorption values must be greater than or equal to 0 and less than 1")
  }

  if(inherits(fidelity, "geoglist")) {
    fidelity <- na.omit(fidelity$layers[[1]], field = names(fidelity$layers[[1]])[1])
    fidelity <- fidelity[[1]][,1]
  } else {
    fidelity <- rep(fidelity, length(unique(tardis$tgraph$src)))
  }
  if(any(fidelity > 1 | fidelity < 0)) {
    stop("fidelity values must be between 0 and 1")
  }

  # transition matrix
  dm <- tardis_to_sparse(tardis, mode = "transition")
  dm <- summary(dm)
  dm <- dm[order(dm$i),]

  # subtract absorption proportionally from row-wise transition probs
  divs <- tapply(dm$x, dm$i, function(x) {x / sum(x)})
  dm$x <- dm$x - (rep(absorption, table(dm$i)) * unlist(divs))

  # format samc matrix
  dm <- sparseMatrix(i = dm$i, j = dm$j, x = dm$x)

  # add fidelity (still unclear how this works, papers are vague)
  diag(dm) <- fidelity

  # add absorption
  dm <- rbind(cbind(dm, absorption), c(rep(0, ncol(dm)), 1))
  colnames(dm) <- rownames(dm) <- c(unique(tardis$tgraph$src), "absorb")

  # to samc (rip of samc code to avoid annoying messages)
  smc <- samc(data = dm)
  # options = samc:::.validate_options(options)
  #
  # r = nrow(data)
  # c = ncol(data)
  #
  # if(c != r) {
  #   stop("Matrix is not square", call. = FALSE)
  # }
  # if(data[r,c] != 1) {
  #   stop("The last element must be 1", call. = FALSE)
  # }
  # if(sum(data[r,]) != 1) {
  #   stop("Last row must be all zeros with a 1 in the last element", call. = FALSE)
  # }
  # if(!isTRUE(all.equal(Matrix::rowSums(data), rep(1, r), check.names = FALSE))) {
  #   stop("All row sums must be equal to 1", call. = FALSE)
  # } # Use all.equal() to avoid numerical precision issues
  #
  #
  # q_mat = methods::as(methods::as(data[-r,-c], "CsparseMatrix"), "generalMatrix")
  # abs_total = data[-r,c]
  #
  # if(!isTRUE(all.equal(Matrix::rowSums(q_mat) + abs_total, rep(1, length(abs_total)), check.names = FALSE))){
  #   stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues
  # }
  #
  # if(is.null(rownames(q_mat)) & is.null(colnames(q_mat))) {
  #   nm = NULL
  # } else if (!isTRUE(all.equal(rownames(q_mat), colnames(q_mat)))) {
  #   stop("The row and col names of the Q matrix must be identical", call. = FALSE)
  # } else if (any(duplicated(rownames(q_mat)))) {
  #   stop("The row and col names of the Q matrix must be unique", call. = FALSE)
  # } else {
  #   nm = rownames(q_mat)
  #   rownames(q_mat) = NULL
  #   colnames(q_mat) = NULL
  # }
  #
  # q_mat@x = -q_mat@x
  # Matrix::diag(q_mat) = Matrix::diag(q_mat) + 1
  #
  # # TODO The clumps value is a placeholder and needs to be calculated as a safety check for the cond_passage() function
  # samc_obj = methods::new("samc",
  #                         data = methods::new("samc_data", f = q_mat, t_abs = abs_total),
  #                         model = list(name = "RW"),
  #                         crw_map = NULL,
  #                         source = "transition",
  #                         nodes = as.integer(r - 1),
  #                         map = terra::rast(),
  #                         names = nm,
  #                         clumps = -1,
  #                         threads = options$threads,
  #                         override = options$override,
  #                         solver = options$method,
  #                         precision = options$precision)
  #
  # samc_obj@.cache$dgf = numeric(nrow(samc_obj@data@f))
  # samc_obj@.cache$dgf_exists = FALSE
  # samc_obj@.cache$sc = samc:::.solver_cache()

  smc$override = T
  regs <- tardis$edges[,1:2]
  regs[,1] <- match(regs[,1], unique(c(regs)))
  regs[,2] <- match(regs[,2], unique(c(regs)))
  smc@clumps <- components(graph_from_edgelist(regs))$membership
  return(smc)
}
