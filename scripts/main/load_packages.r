

PKGS <- c('Seurat', 'dplyr', 'ggplot2', 'data.table')

invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files('scripts/utility', full.names=TRUE, pattern='\\.r$'), source))


data_path <- '/data/CARD_singlecell/Brain_atlas/NABEC_multiome/'

setDTthreads(threads=1)

ngbs <- 50