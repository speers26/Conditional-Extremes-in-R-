library("devtools")

.rs.restartR()
setwd("~/GitHub/R_packages/cond.extremes")
devtools::document()
devtools::install()
