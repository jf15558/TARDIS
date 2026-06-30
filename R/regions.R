#' regions
#'
#' Classify all cells in the layers of a geoglist by which region they
#' belong to. Optionally, only the boundary cells of these regions may be returned.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param bounds.only `logical`. Should only the boundary cells of regions be
#' returned?
#' @param use.links `logical`. Should any links added to the geoglist also be
#' accounted for when determining cell connectivity and so region membership?
#' This argument will have no effect if the geoglist does not contain any links.
#' @return A geoglist recording the region affinity of the cells in each layer.
#' @import sf terra
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph components
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(TARDIS)
#'
#' gal <- cretaceous()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' rasts <- rast_to_geoglist(gal, gal_m)
#'
#' regs <- regions(rasts)
#' plot(regs)
#' }

regions <- function(geog, bounds.only = T, use.links = F) {
  #
   #geog = rasts
   #bounds.only = F
   #use.links = T

  if(!exists("geog")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }

  if(inherits(geog$layers, "SpatRaster")) {

    ils <- rast(lapply(geog$layers, patches, directions = 8, allowGaps = F))
    if(bounds.only) {
      ils <- rast(lapply(1:nlyr(ils), function(x) {
        foo <- mask(ils[[x]], classify(boundaries(ils[[x]]), cbind(0, NA)))
        names(foo) <- "region"
        foo
      }))
    }

    if(use.links) {
      ils <- rast(lapply(1:length(ils), function(x) {
        lnk <- geog$links[which(geog$links$layer == x)]
        if(length(lnk) != 0) {
          ig <- extract(ils[[x]], lnk)
          ig <- ig[complete.cases(ig),]
          ig2 <- graph_from_edgelist(as.matrix(ig), directed = F)
          ils[[x]] <- classify(ils[[x]], cbind(ig[,1], components(ig2)$membership))
        }
        ils[[x]]
      }))
    }

  } else {

    ils <- lapply(geog$layers, function(z) {

      orig <- length(z)
      z$id <- 1:orig
      z <- na.omit(z, field = names(z)[1])
      bar <- relate(z, z, "intersects", pairs = T)
      z[,1] <- components(graph_from_edgelist(bar))$membership
      if(bounds.only) {
        z <- z[which(table(bar[,1]) != 7)]
      }
      names(z)[1] <- "region"
      baz <- vect(cbind(1:(orig - length(z)), 1, NA, NA, 0), type = "polygon")
      baz$region <- rep(NA, length(baz))
      baz$id <- setdiff(1:orig, z$id)
      baz <- rbind(z, baz)
      baz$region[is.nan(baz$region)] <- NA
      baz <- baz[order(baz$id)]
      baz[,1]
    })

    if(use.links) {
      ils <- lapply(1:length(ils), function(x) {
        lnk <- geog$links[which(geog$links$layer == x)]
        ig <- cbind(ils[[x]]$region[lnk$srt], ils[[x]]$region[lnk$end])
        ig2 <- graph_from_edgelist(ig, directed = F)
        ils[[x]]$region[!is.na(ils[[x]]$region)] <- components(ig2)$membership[ ils[[x]]$region[!is.na(ils[[x]]$region)]]
        ils[[x]]
      })
    }
  }
  if(!inherits(geog$layers, "SpatRaster")) {ils <- svc(ils)}

  geog$layers <- ils
  return(geog)
}
