#!/usr/bin/env Rscript

args    <- commandArgs(trailingOnly = TRUE)
target  <- args[1]
lib     <- args[2]
threads <- args[3]

options(Ncpus=threads)
.libPaths(lib)

devtools::install_github(
  target,
  dependencies = NA,
  upgrade      = 'never'
)

deps <- c('patchwork', 'DT', 'ggrepel', 'pheatmap')

devtools::install_cran(
  deps,
  dependencies = NA,
  upgrade      = 'never'
)
