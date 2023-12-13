# SETUP
#
# reset environment
rm(list = ls())
# set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# GENERATE
#usethis::create_package("./")
#usethis::use_git()
#usethis::use_data_raw()
usethis::use_package("nngeo", "cppRouting", "parallel", "pbapply", "raster",
                     "geosphere", "sf", "concaveman", "igraph")
