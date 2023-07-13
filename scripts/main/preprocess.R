
working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')


reticulate::source_python('scripts/utility/scrublet_py.py')

dataset <- snakemake@params[['dataset']]
batch <- fread(snakemake@input[['datasets']])[sample %chin% dataset, batch]

raw.counts <- Read10X_h5(paste0(data_path, batch, '/Multiome/', dataset, '/outs/raw_feature_bc_matrix.h5'))
filtered.counts <- Read10X_h5(paste0(data_path, batch, '/Multiome/', dataset, '/outs/filtered_feature_bc_matrix.h5'))

adj.matrix <- suppressWarnings(SoupCorrect(raw.counts, filtered.counts, contamination_rate=snakemake@params[['soup_rate']]))
object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)


object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

m <- copy(object@meta.data)
setDT(m, keep.rownames='cells')

m[,
    `:=` (
        sample=dataset,
        batch=batch
        )
]

batch <- m[, batch]
sample <- m[, sample]

names(sample) <- names(batch) <- m[, cells]


doublet_rate <- (ncol(object) / 1000) * 0.008

object <- object %>% 
    
    AddMetaData(metadata=factor(batch), col.name='batch') %>% 
    AddMetaData(metadata=factor(sample), col.name='sample') %>% 
    
    scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate) 


saveRDS(object, snakemake@output[['seurat_object']])






