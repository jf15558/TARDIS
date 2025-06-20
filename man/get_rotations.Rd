% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_rotations.R
\name{get_rotations}
\alias{get_rotations}
\title{get_rotations}
\usage{
get_rotations(
  rast,
  ages,
  model = "PALEOMAP",
  out.res = NULL,
  gpath = NULL,
  verbose = TRUE
)
}
\arguments{
\item{rast}{A SpatRaster in longitude-latitude projection}

\item{ages}{A numeric vector of desired reconstruction ages in millions of
years before present, in descending age order. These will correspond to the
layers in the SpatRaster destined for create_tardis().}

\item{model}{The desired plate reconstruction model. See
rgplates::reconstruct() for details.}

\item{out.res}{If you are calculating rotations at a reduced resolution to
save time, then out.res should be the SpatRaster you intent to use with
create_tardis(). This is essential so that the returned cell IDs will conform
to that SpatRaster, rather than to that given by 'rast'.}

\item{gpath}{If NULL, the GPlates API is used. This is very slow for large
numbers of points, so the software call be installed locally for faster
results. If desired gpath is the file path to the local reconstruction
submodule. In case the GPlates executable file is not found at the coded
default location, the full path to the executable (gplates-\if{html}{\out{<ver>}}.exe on
Windows) can be entered, e.g. "C:/gplates_2.3.0_win64/gplates.exe".}

\item{verbose}{Should function progress be reported?}
}
\value{
A list with length(ages) - 1 elements. Each element is a two-column
matrix containing the IDs of geographically homologous cells between
successive reconstruction ages.
}
\description{
Generate a rotation list using commonly-implemented plate rotation models for
use with create_tardis(). The user supplies a SpatRaster with the desired
geographic resolution for the resulting rotation list, along with a vector of
reconstruction ages in descending order. Conceivably, the SpatRaster and age
vector the user intends to supply to create_tardis() could be used here. The
function will take the cell centre coordinates for all grid cells in the
SpatRaster as starting positions at the first age in the vector, then
calculates their geographic positions at the subsequent ages. Calculation
relies on the reconstruct() function from the rgplates package and will be
very slow if high resolution SpatRasters/lengthy time vectors are supplied.
Each SpatRaster layer must have at least some reconstructed grid cell
positions so the age vector cannot be made shorter. Instead, it may be
advisable to use a SpatRaster with a lower geographic resolution than those
the user intends to supply to create_tardis(). This is totally fine as the
rotations list does not need to give reconstructions for every grid cell in
the landscape raster stack, nor will this even be possible as gplates will be
unable to determine some geographic positions due to the limitations of the
plate rotation models themselves. NAs are returned where reconstructions fail,
although this is again not an issue as these will be automatically filtered
by create_tardis(). There is no point performing rotations at a higher
geographic resolution than the SpatRaster stack that you will use with
create_tardis()
#'
}
\examples{
#rast <- rast(nrows = 2, ncols = 2)
#ages <- c(440, 430, 420)
#foo <- get_rotations(rast = rast, ages = ages)

}
