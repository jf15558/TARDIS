# SETUP
#
# reset environment
rm(list = ls())
# set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(nngeo)
library(cppRouting)
library(parallel)
library(pbapply)
library(raster)
library(sf)
library(terra)
library(geosphere)
library(igraph)
library(concaveman)
library(ape)
library(Matrix)
library(akima)
library(terra)

# functions
#lapply(list.files("./R", full.names = T), source)

# map data
load("data/Chelonoides_tree.rda")
gal <- terra::rast(lapply(readRDS("data-raw/galap_hr.RData"), terra::rast))
gal <- crop(gal, extent(-92, -88, -2, 1))
gal_m <- classify((gal), rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
gal_l <- classify((gal), rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = T), right = F)
org <- rbind(c(-89, -1.05, 2), c(-89.5, -0.7, 2))
dst <- rbind(c(-91.2, -1, 0), c(-91.6, -0.4, 0))

# tree data
tdat <- Chelonoides_tree$biogeography
rt <- cbind.data.frame(lng = rep(-89.5, 13), lat = rep(-0.9, 13), age = rep(1.5, 13))

cairo_pdf("test.pdf", height = 10, width = 7)
par(mfrow = c(3, 2))
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(gal_m[[1]], axes = F, legend = F)
plot(gal_m[[2]], axes = F, legend = F)
plot(gal_m[[3]], axes = F, legend = F)
plot(gal_m[[4]], axes = F, legend = F)
plot(gal_m[[5]], axes = F, legend = F)
dev.off()

plot(gal_m[[1]])
plot(test3$geometry, add = T)
plot(test4$geometry, add = T, col = scales::alpha(1, 0.2), border = F)
plot(gal_m[[5]])
plot(test3$geometry, add = T)
plot(test4$geometry, add = T, col = scales::alpha(1, 0.2), border = F)



# functions

test <- link_mask(gal_m)
test2 <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
test_weight <- weight_tardis(test2, vars = list(elev = classify(gal, cbind(-Inf, 0, 0))),
                             wfun = function(origin, dest, lnum, ...) {(origin$hdist^2 + abs(origin$vdist)^2)},
                             mfun = function(origin, dest, lnum, ...) {(origin$hdist^2 + abs(origin$vdist)^2) * 10})
rst <- resistance_surface(test2)
lapply(1:nlyr(gal_m), function(x) {plot(gal_m[[x]]); plot(test[[x]], add = T)})


pts <- stp(test2, rbind(org, dst))
test3 <- lcp(tardis = test2, weights = test_weight, pts[1:2,], pts[3:4,])
test4 <- detour(tardis = test2, weights = test_weight, test3, detour = 0.01)
test5 <- commute(tardis = test2, weights = test_weight, pts[1:2,], pts[3:4,])
test6 <- isochrone(tardis = test2, weights = test_weight, pts[3:4,], cost = 1e5)
test7 <- random_walk(tardis = test2, weights = test_weight, pts[3:4,], rwlen = 1e6)

cairo_pdf("test2.pdf", width = 7, height = 10)
par(mfrow = c(2, 1))
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(gal_m[[1]], axes = F, legend = F)
plot(test[[1]], add = T, col = "grey90")
plot(st_union(test4$geometry), add = T, col = scales::alpha(1, 0.1), border = F)
plot(test3$geometry, add = T, col = "red", lwd = 2)
plot(gal_m[[5]], axes = F, legend = F)
plot(test[[5]], add = T, col = "grey90")
plot(test3$geometry, add = T)
plot(st_union(test4$geometry), add = T, col = scales::alpha(1, 0.1), border = F)
plot(test3$geometry, add = T, col = "red", lwd = 2)
dev.off()


pairs <- rbind(tdat[1:13,c("lng", "lat", "age")], rt)
matched <- stp(test2, pairs)

tort <- lcp(tardis = test2, weights = test_weight, origin = matched[14:26,], dest = matched[1:13,])
tort_d <- detour(tardis = test2, weights = test_weight, paths = tort, detour = 0.01)

cairo_pdf("test3.pdf", width = 7, height = 10)
par(mfrow = c(2, 1))
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(gal_m[[2]], axes = F, legend = F)
plot(test[[2]], add = T, col = "grey90")
plot(st_union(tort_d$geometry), add = T, col = scales::alpha(1, 0.1), border = NA)
plot(tort$geometry, add = T, col = "red", lwd = 2)
plot(gal_m[[5]], axes = F, legend = F)
plot(test[[5]], add = T, col = "grey90")
plot(st_union(tort_d$geometry), add = T, col = scales::alpha(1, 0.1), border = NA)
plot(tort$geometry, add = T, col = "red", lwd = 2)
dev.off()


phylo_path <- function(tree, paths, nodes) {

  # tree <- arch_tree
  # paths <- tres4
  # nodes <- c(393, which(arch_tree$tip.label == "Hyperodapedon"))

  if(!exists("tree")) {
    stop("Please supply tree as an object of class 'phylo'")
  }
  if(class(tree) != "phylo") {
    stop("Please supply tree as an object of class 'phylo'")
  }

  if(!exists("paths")) {
    stop("Supply paths as the output of tardis()")
  }
  if(!is.list(paths) | length(paths) != 3) {
    stop("Supply paths as the output of tardis()")
  }
  if(!all(names(paths) %in% c("cells", "paths", "values"))) {
    stop("Supply paths as the output of tardis()")
  }

  if(!exists("nodes")) {
    stop("Please supply nodes as a vector of two or more internal and/or tip node ids")
  }
  if(!is.vector(nodes) | any(is.na(nodes)) | length(nodes) < 2) {
    stop("Please supply nodes as a vector of two or more internal and/or tip node ids")
  }
  if(is.character(nodes)) {
    if(!all(nodes %in% tree$tip.label)) {
      stop("One or more node names is not present in the tree tip labels")
    }
    nodes <- match(nodes, tree$tip.label)
  }
  if(is.numeric(nodes)) {
    if(any(nodes > (Ntip(tree) + Nnode(tree)))) {
      stop("One or more node numbers exceeds the number of nodes in the tree")
    }
  }
  if(length(paths$cells) != Nedge(tree)) {
    stop("The number of paths is not equal to the number of edges in the tree")
  }

  if(any(nodes <= Ntip(tree)) & !((Ntip(tree) + 1) %in% nodes)) {
    nodes <- c(getMRCA(tree, nodes), nodes)
  }
  nodes <- nodes[order(nodes)]
  mrca <- nodes[which(nodes > Ntip(tree))[1]]
  mlist <- list()
  for(i in 1:length(nodes)) {
    pth <- nodepath(tree, mrca, nodes[i])
    if(length(pth) != 1) {
      pth2 <- vector()
      for(j in 1:(length(pth) - 1)) {
        pth2[j] <- which(tree$edge[,1] == pth[j] & tree$edge[,2] == pth[j + 1])
      }
      mlist[[i]] <- pth2
    }
  }
  return(paths$paths[which(paths$paths$path %in% unique(unlist(mlist))),])
}




tardis <- test2
tstep <- seq(tardis$tdat[1], rev(tardis$tdat)[1], -0.1)
i=15
for(i in 1:length(tstep)) {

  trange <- c((sum(tstep[i] <= tardis$tdat) - 1) * prod(tardis$gdat[1:2]) + 1,
               sum(tstep[i] <= tardis$tdat) * prod(tardis$gdat[1:2]))

  # weights
  weights <- tardis$edges[,5]

  # remove NA
  if(any(is.na(weights))) {
    tardis$tgraph$src <- tardis$tgraph$src[which(!is.na(weights))]
    tardis$tgraph$dst <- tardis$tgraph$dst[which(!is.na(weights))]
    tardis$edges <- tardis$edges[which(!is.na(weights)),]
    tardis$tgraph$dict <- tardis$tgraph$dict[which(tardis$tgraph$dict$ref %in% c(tardis$tgraph$src, tardis$tgraph$dst))]
    tardis$tgraph$id <- 0:(length(tardis$tgraph$dict$ref) - 1)
  }
  tardis$tgraph$attrib$aux <- tardis$edges[!is.na(weights),5]
  tardis$edges[,5] <- weights[!is.na(weights)]

  # time slice
  grp <- tardis$tgraph
  inbin <- which(tardis$edges[,1] >= trange[1] & tardis$edges[,1] <= trange[2] & tardis$edges[,5] != 0)
  grp$attrib$aux <- grp$attrib$aux[inbin]
  src <- tardis$tgraph$src[inbin]
  dst <- tardis$tgraph$dst[inbin]
  cst <- tardis$edges[inbin,5]
  grp$data <- data.frame(from = grp$dict$id[match(src, grp$dict$ref)], to = grp$dict$id[match(dst, grp$dict$ref)], dist = cst)
  grp$dict <- grp$dict[which(grp$dict$id %in% grp$data$from),]

  grp$dict$id <- grp$dict$id - grp$data$from[1]
  grp$data$to <- grp$data$to - grp$data$from[1]
  grp$data$from <- grp$data$from - grp$data$from[1]

  # cppRouting graph and sparse matrix
  grp <- grp[1:5]
  adj_mat <- sparseMatrix(i = grp$data$from + 1, j = grp$data$to + 1, x = grp$data$dist, symmetric = F)


  #smp <- sample(grp$dict$ref, 10)
  taxon <- setNames(ceiling(rnorm(10, 50)), smp)
  tmprast <- raster(nrows = tardis$gdat[1], ncols = tardis$gdat[2])
  tmprast <- gal_m[[3]]
  tmprast[as.numeric(names(taxon)) %% prod(tardis$gdat[1:2])] <- 2
  plot(tmprast)

  samprast <- raster(nrows = tardis$gdat[1], ncols = tardis$gdat[2])
  for(j in 1:20) {

    tmprast <- samprast
    loc <- as.numeric(names(taxon)) %% prod(tardis$gdat[1:2])
    tmprast[loc] <- 1
    tmprast <- (clump(tmprast))
    tmprast2 <- boundaries(tmprast)
    pioneer <- names(taxon)[which(loc %in% which(tmprast2[] == 1))]

    disp <- abs(rnorm(length(pioneer), mean = mean(grp$data$dist), sd = 1000))
    dest <- unlist(mapply(x = pioneer, y = disp, function(x, y) {
      inreach <- get_isochrone(grp, x, y)
    }))
    pop <- ceiling(rnorm(length(dest), mean = 50))
    taxon <- tapply(c(taxon, pop), c(names(taxon), dest), sum)

    tmprast <- gal_m[[3]]
    tmprast[as.numeric(names(taxon)) %% prod(tardis$gdat[1:2])] <- 2
    plot(tmprast)
  }


  tmprast <- samprast
  loc <- as.numeric(names(taxon)) %% prod(tardis$gdat[1:2])
  tmprast[loc] <- 1
  tmprast <- (clump(tmprast))
  tmprast2 <- boundaries(tmprast)
  pioneer <- names(taxon)[which(loc %in% which(tmprast2[] == 1))]
  pdm <- get_distance_matrix(grp, pioneer, pioneer)
  pioneer2 <- as.numeric(pioneer) %% prod(tardis$gdat[1:2])

  res <- do.call(rbind, lapply(1:nrow(pdm), function(x) {tapply(pdm[x,], tmprast[pioneer2], min)}))

  pdm <- hclust(as.dist(pdm))
  pops <- cutree(pdm, h = 10000)

}



# other archipelagos here: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13220#geb13220-bib-0037?saml_referrer





source("C:/Users/jf15558/OneDrive - University of Bristol/Documents/func.txt")


# map data (source: https://doi.plutof.ut.ee/doi/10.15156/BIO/786389)
if(TRUE) {
  foo <- lapply(rev(gtools::mixedsort(list.files("galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/", pattern = "txt", full.names = T))), function(x) {


    x <- rev(gtools::mixedsort(list.files("galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/", pattern = "txt", full.names = T)))[1]

dat$v7 <- 1:nrow(dat)

dat <- dat[order(dat$V1, dat$V2, method = "radix", decreasing = c(FALSE, FALSE)),]

    dat <- read.table(x)
    xy <- as.matrix(dat[,1:2])
    v <- dat$V5
    i <- !is.na(v)
    xy <- xy[i,]
    v2 <- v[i]



    #ob <- raster(nrows = length(unique(dat$V2)), ncols = nrow(dat) / length(unique(dat$V2)))
    #ob[] <- dat$V5


    ob <- interp(x = xy[,1], y = xy[,2], z = v2, xo = seq(min(xy[,1]), max(xy[,1]), 0.01),
                 yo = seq(min(xy[,2]), max(xy[,2]), 0.01))
    ob3 <- raster(t(ob$z[,ncol(ob$z):1]), xmn = min(ob$x) - 0.05, xmx = max(ob$x) + 0.05,
                  ymn = min(ob$y) - 0.05, ymx = max(ob$y) + 0.05)
    ob3 <- resample(ob3, raster(res = c(0.01, 0.01), ext = raster::extent(ob3)))
  })
  #ob <- curl::curl_download("ftp://ftp.pmel.noaa.gov/newport/chadwick/galap_bathy/galap.grd", destfile = "galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/mod/modtest.grd")
  ob <- raster("galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/mod/modtest.grd")
  ob <- crop(ob, raster::extent(foo[[1]]))
  ob <- resample(ob, foo[[4]])
  foo2 <- c(foo, ob)
  saveRDS(foo2, "galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/galap_hr.RData")

}


hexblender <- function(color1, color2, alpha = 0.5) {

  # source https://rdrr.io/github/lawine90/elseR/src/R/hexBlender.R

  col1_rgb <- colorspace::hex2RGB(color1)
  col2_rgb <- colorspace::hex2RGB(color2)

  result_rgb <- colorspace::mixcolor(
    alpha,
    colorspace::RGB(col1_rgb@coords[,1],
                    col1_rgb@coords[,2],
                    col1_rgb@coords[,3]),
    colorspace::RGB(col2_rgb@coords[,1],
                    col2_rgb@coords[,2],
                    col2_rgb@coords[,3]))

  result_hex <- grDevices::rgb(result_rgb@coords[,1]*255,
                               result_rgb@coords[,2]*255,
                               result_rgb@coords[,3]*255,
                               maxColorValue = 255)
  return(result_hex)
}
scalerfunc <- function(x) {return((x + abs(min(x, na.rm = T))) / (max(x, na.rm = T) + abs(min(x, na.rm = T))))}
hillshader <- colorRamp(grey(0:100/100))
toposhader <- colorRamp(c("bisque", "burlywood", "lightsalmon3", "lightsalmon4", "white"))
bathshader <- colorRamp(c("cornflowerblue", "lightskyblue1", "lightcyan1"))

# data
maps <- list()
for(i in 1:nlayers(gal)) {

  mp1 <- t(as.matrix(flip(gal[[i]])))

  # topo colours
  bath <- topo <- mp1
  topo[which(topo  < 0)] <- NA
  bath[which(bath >= 0)] <- NA
  topo <- toposhader(scalerfunc(topo))
  bath <- bathshader(scalerfunc(bath))
  bath[which(!complete.cases(bath)),] <- topo[which(complete.cases(topo)),]
  elev <- rgb(bath, maxColorValue = 255)

  # hill shade
  hill <- as.matrix(hillShade(terrain(raster(mp1, crs = crs(gal[[1]])), "slope"),
                              terrain(raster(mp1, crs = crs(gal[[1]])), "aspect")))
  hill[which(!is.finite(hill))] <- 0

  hill <- hillshader(scalerfunc(hill))
  hill <- rgb(hill, maxColorValue = 255)
  hill[which(mp1 < 0)] <- "#FFFFFF"

  # compound colour
  finalcol <- col2rgb(hexblender(hill, elev))

  bar <- raster(ncol = ncol(gal), nrow = nrow(gal))
  bar[] <- 1
  bar2 <- stack(bar, bar, bar)
  bar2[[1]][] <- finalcol[1,]
  bar2[[2]][] <- finalcol[2,]
  bar2[[3]][] <- finalcol[3,]
  bar2 <- flip(bar2)
  print(i)
}


tmp <- bar2
image(tmp, col = NA, axes = F, asp = 1)
plotRGB(tmp, add = T, axes = F, box = F)

