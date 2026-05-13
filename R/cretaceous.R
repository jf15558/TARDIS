#' cretaceous
#'
#' Load a terra 'SpatRaster' object from a source tiff containing two
#' palaeogeographic digital elevation models of the Earth during the Cretaceous
#' Models are in geographic coordinates with elevations in metres above sea level.
#' @export
#'
#' @details Each DEM records the estimated topography and bathymetry of the Earth
#' at 115 and 110 million years in the past. As such, each time slice is taken
#' to represent the 2.5 million years either side of these points (i.e., the
#' 110 Ma reconstruction spans 112.5-107.5 Ma). Both reconstructions come from
#' the PALEOMAP plote rotation model series by Chris Scotese and colleages and
#' were obtained using the `fetch()` function of the `chronosphere` R package
#'
#' @source Kocsis, Á.T., Raja, N.B. 2013 chronosphere: Evolving Earth System
#' Variables. Zenodo. https://doi.org/10.5281/zenodo.3525481.
#' 
#' Scotese, C. R., Vérard, C., Burgener, L., Elling, R. P., & Kocsis,
#' Á. T. (2024). Phanerozoic-scope supplementary material to 'The Cretaceous
#' World: Plate Tectonics, Paleogeography, and Paleoclimate' from the PALEOMAP
#' project (Version v24221). Zenodo. https://doi.org/10.5281/zenodo.10659112.
#' 
#' Scotese, C. R., Vérard, C., Burgener, L., Elling, R. P., & Kocsis, A. T.
#' (2025). The Cretaceous World: Plate Tectonics, Paleogeography, and
#' Paleoclimate. Geological Society, London, Special Publications, 544(1),
#' SP544–2024.

cretaceous <- function() {
  rast(system.file("extdata", "cretaceous.tif", package = "TARDIS"))
}