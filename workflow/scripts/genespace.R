devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

# get mcscanx path from args
path2mcscanx <- commandArgs(trailingOnly = TRUE)[1]

gpar <- init_genespace(wd = "results/", path2mcscanx = path2mcscanx)