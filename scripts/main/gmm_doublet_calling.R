working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


object.list <- future.apply::future_lapply(snakemake@input[['seurat_object']], readRDS)


m <- rbindlist(lapply(object.list, function(object) {
    m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
}))

# cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, doublet_scores], k=2))
cutoff <- 0.06

m[, `:=` (
    project=rep(snakemake@params[['project_name']], nrow(m)),
    predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
)]

fwrite(x=m, file=snakemake@output[['metadata']])

