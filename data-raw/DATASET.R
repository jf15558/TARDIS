# SETUP
#
# reset environment
rm(list = ls())
# set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# MAPS
#
# past geography
past <- lapply(rev(gtools::mixedsort(list.files("./", pattern = "txt", full.names = T)))[-5], function(x) {

  # source downloaded manually to wd: https://doi.plutof.ut.ee/doi/10.15156/BIO/786389)
  #dat <- rev(gtools::mixedsort(list.files("galapagos/Karnauskas_etal_FinalModelOutput_Bathymetry/", pattern = "txt", full.names = T)))[1]
  dat <- read.table(x)
  dat$v7 <- 1:nrow(dat)
  dat <- dat[order(dat$V1, dat$V2, method = "radix", decreasing = c(FALSE, FALSE)),]
  xy <- as.matrix(dat[,1:2])
  v <- dat$V5
  i <- !is.na(v)
  xy <- xy[i,]
  v2 <- v[i]

  # interpolate to regular grid
  ob <- interp::interp(x = xy[,1], y = xy[,2], z = v2, xo = seq(min(xy[,1]), max(xy[,1]), 0.01),
               yo = seq(min(xy[,2]), max(xy[,2]), 0.01))
  ob3 <- terra::rast(t(ob$z[,ncol(ob$z):1]))
  terra::ext(ob3) <- c(min(ob$x) - 0.05, max(ob$x) + 0.05, min(ob$y) - 0.05, max(ob$y) + 0.05)
  ob3 <- terra::resample(ob3, rast(res = c(0.01, 0.01), ext = ob3))
})

# present geography
mod <- curl::curl_download("ftp://ftp.pmel.noaa.gov/newport/chadwick/galap_bathy/galap.grd", destfile = "./modern.grd")
mod <- terra::rast("./modern.grd")
mod <- terra::crop(mod, terra::extent(past[[1]]))
mod <- terra::resample(ob, past[[4]])
galapagos <- terra::rast(lapply(c(past, mod), terra::rast))
names(galapagos) <- c("2.25-1.75", "1.75-1.25", "1.25-0.75", "0.75-0.25", "0.25-0")
galapagos <- terra::wrap(galapagos)

# save

foo <- readRDS("galap_hr.RData")
foo <- rast(lapply(foo, terra::rast))
foo <- terra::wrap(foo)
galapagos <- foo
#usethis::use_data(galapagos, internal = F)
writeRaster(galapagos, "../inst/extdata/galapagos.tif")


# PHYLOGENY
#
# modern distr: https://www.nature.com/articles/s42003-022-03483-w
# past origins: https://onlinelibrary.wiley.com/doi/full/10.1111/jzs.12387?saml_referrer
chelonoides <- ape::read.tree(text = "(((becki, darwini), ((donfaustoni, chathamensis), (abingdonii, hoodensis))), (((((guntheri, vicina), (microphyes, vandenburghi)), porteri), niger), duncanensis));")
chelonoides$edge.length <- c(0.45, 1.06, 0.03, 0.03, 0.37, 0.39, 0.33, 0.33, 0.47, 0.25, 0.25,
                     0.52, 0.47, 0.14, 0.28, 0.09, 0.04, 0.04, 0.08, 0.05, 0.05, 0.41, 0.55, 1.02)
chelonoides$root.time <- 1.54

# biogeographic data
tdat <- cbind.data.frame(node = c(chelonoides$tip.label, as.character(14:25)),
                         island =  c("Isa", "Santi", "SantaC", "SanC", "Pin", "Esp", "Isa", "Isa", "Isa", "Isa", "SantaC", "Flo", "Pin",
                                     "Esp-SanC", "SanC", "Santi", "SanC", "SanC", "Esp", "Pin", "SantaC", "SantaC", "Isa", "Isa", "Isa"),
                         lng = c(-91.36, -91.33, -90.3, -89.4, -90.76, -89.65, -90.94,
                                 -91.3, -91.34, -91.1, -90.41, -90.43, -90.66, rep(NA, 12)),
                         lat = c(0,      -0.2,   -0.6, -0.85,  0.6,   -1.38,   -0.75,
                                 -0.96, -0.15, -0.45, -0.69, -1.3, -0.61, rep(NA, 12)),
                         age = chelonoides$root.time - ape::node.depth.edgelength(chelonoides))
chelonoides$biogeography <- tdat

# save
usethis::use_data(chelonoides, overwrite = TRUE)

