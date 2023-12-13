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
#usethis::use_package("nngeo")
#usethis::use_package("cppRouting")
#usethis::use_package("parallel")
# usethis::use_package("pbapply")
# usethis::use_package("terra")
# usethis::use_package("geosphere")
# usethis::use_package("sf")
# usethis::use_package("concaveman")
# usethis::use_package("igraph")
# usethis::use_package("Matrix")

devtools::check()
devtools::build()
