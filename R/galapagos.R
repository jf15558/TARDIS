#' galapagos
#'
#' Load a terra 'SpatRaster' object from a source tiff containing four
#' palaeogeographic digital elevation models and one modern digital elevation
#' model of the Galapagos archipelago. Models are in geographic coordinates with
#' elevations in metres above sea level.
#' @export
#'
#' @details Each DEM records the topography and bathymetry of the archipelago
#' at 0.5 million-year intervals from 2 million years ago to the present. As
#' such, each time slice is taken to represent the 0.25 million years either
#' side of these average points (i.e., the 2 Ma reconstruction spans 2.25-1.75
#' Ma).
#'
#' The palaeogeographic reconstructions were interpolated from an irregular grid
#' and additionally contain artefacts. As such, they should only be treated as
#' toy datasets, rather than wholly accurate, and later versions are expected to
#' be released in the future.
#'
#' The present day reconstruction was sourced from the Pacific Marine
#' Environmental Laboratory from data comppiled by William Chadwick, Oregon
#' State University, although it has been downscaled to palaeogeographic
#' resolution and so should again be treated as an example dataset, rather than
#' the basis for a scientifically robust analysis.
#'
#' @source Past DEM: supplied by the original authors from Karnauskas et al
#' (2017). Paleoceanography of the eastern equatorial Pacific over the past 4
#' million years and the geologic origins of modern Galápagos upwelling. Earth
#' and Planetary Science Letters, 460, 22-28
#' @source Modern DEM: documented at https://www.pmel.noaa.gov/eoi/staff/chadwick/galapagos.html,
#' downloaded from the embedded link in the page: ftp://ftp.pmel.noaa.gov/newport/chadwick/galap_bathy/gala
galapagos <- function() {
  rast(system.file("extdata", "galapagos.tif", package = "TARDIS"))
}
