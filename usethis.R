# SETUP
#
# reset environment
rm(list = ls())
# set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# GENERATE
usethis::create_package("./")
