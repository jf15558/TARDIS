% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slice_tardis.R
\name{slice_tardis}
\alias{slice_tardis}
\title{slice_tardis}
\usage{
slice_tardis(tardis, times = NULL, layers = NULL)
}
\arguments{
\item{tardis}{An object of class 'tardis', produced by create_tardis}

\item{times}{A vector of two positive numbers denoting the desired time range
to subset from tardis. Either this or 'layers' must be specified}

\item{layers}{A vector of two positive numbers denoting the desired range of
layers to subset from tardis. Either this or 'times' must be specified}
}
\value{
An object of class 'tardis', comprising the requested range of layers
}
\description{
Temporally subset a TARDIS graph This function was developed for instances
where successive analyses do not require the entire graph, in which case it
is more efficient to weight and analyse subsets, rather than operate on the
entire graph or create subsets from scratch. A user can subset by specific
geographic layers, or supply a time range which will be used to select the
geographic layers within that range. Note that layers are counted in
decreasing age order, so the oldest time layer will be 1 and so forth.
}
\examples{
#library(terra)
#library(TARDIS)

#gal <- galapagos()
#gal <- crop(gal, extent(-92, -88, -2, 1))
#gal_m <- classify(gal, rcl = matrix(c(-Inf, 0, NA, 0, Inf, 1),
#                                    ncol = 3, byrow = T), right = F)
#gt <- create_tardis(gal, times = c(seq(2.25, 0, -0.5), 0), mask = gal_m)
#gts <- slice_tardis(gt, times = c(1.2, 0))
#gts <- slice_tardis(gt, layers = c(1, 2))
}
