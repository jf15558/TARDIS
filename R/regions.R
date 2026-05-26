#' regions
#' 
#' Classify all cells in the layers of a geoglist by which region they
#' belong to. Optionally, only the boundary cells of these regions may be returned.
#' 
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param coasts.only `logical`. Should only the boundary cells of regions be
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
#' gal <- cretaceous()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#' rasts <- rast_to_geoglist(gal, gal_m)
#' rasts <- link_islands(rasts, klink = 1)
#' 
#' regs <- regions(geog = rasts)

regions <- function(geog, bounds.only = T, use.links = F) {
  # 
  # geog = rasts
  # bounds.only = F
  # use.links = F
  
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
      ils <- rast(lapply(1:nlyr(ils), function(x) {
        lnk <- geog$links$geometry[which(geog$links$layer == x)]
        if(length(lnk) != 0) {
          ig <- extract(ils[[x]], vect(lnk))
          ig <- ig[complete.cases(ig),]
          ig2 <- graph_from_edgelist(as.matrix(ig), directed = F)
          ils[[x]] <- classify(ils[[x]], cbind(ig[,1], components(ig2)$membership))
        }
        ils[[x]]
      }))
    }
    
  } else {
    
    ils <- lapply(geog$layers, function(z) {
      bar <- st_touches(z)
      bar2 <- do.call(rbind, sapply(1:length(bar), function(x) {rbind(c(x, x), cbind(rep(x, length(bar[[x]])), bar[[x]]))}))
      z[,1] <- components(graph_from_edgelist(bar2))$membership
      colnames(z)[1] <- "region"
      if(bounds.only) {
        z <- z[sapply(bar, length) != 6,]
      }
      z
    })
    
    if(use.links) {
      ils <- lapply(1:length(ils), function(x) {
        lnk <- geog$links$geometry[which(geog$links$layer == x)]
        if(length(lnk) != 0) {
          ig <- t(sapply(st_intersects(lnk, ils[[x]]), function(y) {ils[[x]]$region[y]}))
          ig2 <- graph_from_edgelist(ig, directed = F)
          ils[[x]]$region <- components(ig2)$membership[ils[[x]]$region]
        }
        ils[[x]]
      })
    }
  }
  
  geog$layers <- ils
  return(geog)
}
