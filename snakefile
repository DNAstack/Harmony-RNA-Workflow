import pandas

num_workers = 2

project = 'card_sc_brain_atlas'

clustering_algorithm = 3
clustering_resolution = 0.3

input_table = 'input/samples.csv'
datasets = pandas.read_csv(input_table, header=None).loc[:, 0].tolist()[1:]

groups = ['sample', 'batch', 'seurat_clusters']
features = ['doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb']


rule all:
    input:
        'plots/qc_plot1.pdf', 
        'plots/qc_plot2.pdf',
        'output/final_metadata.csv',
        'output/sctype_cluster_annotations.csv',
        expand('plots/{group}_group_umap.pdf', group=groups),
        expand('plots/{feature}_feature_umap.pdf', feature=features),
        'objects/seurat_object_harmony_integrated_neighbors_umap_cluster_07.rds'


rule preprocess:
    input:
        datasets=input_table
    output:
        seurat_object='objects/seurat_object_{dataset}_preprocessed_01.rds'
    params:
        soup_rate=0.20,
        dataset='{dataset}'
    conda:
        'envs/multiome.yml'
    script: 
        'scripts/main/preprocess.R'

rule doublets:
    input:
        seurat_object=expand('objects/seurat_object_{dataset}_preprocessed_01.rds', dataset=datasets)
    output:
        metadata='output/unfiltered_metadata.csv'
    params:
        project_name=project
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/gmm_doublet_calling.R'

rule plot_qc:
    input:
        metadata='output/unfiltered_metadata.csv'
    output:
        plot_1='plots/qc_plot1.pdf', plot_2='plots/qc_plot2.pdf'
    params:
        project_name=project
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/plot_qc_metrics.R'

rule filter:
    input:
        metadata='output/unfiltered_metadata.csv',        
        seurat_object='objects/seurat_object_{dataset}_preprocessed_01.rds'
    output:
        seurat_object='objects/seurat_object_{dataset}_preprocessed_filtered_02.rds'
    conda:
        'envs/multiome.yml'
    script: 
        'scripts/main/filter.R'

rule process:
    input:
        seurat_object='objects/seurat_object_{dataset}_preprocessed_filtered_02.rds'
    output:
        seurat_object='objects/seurat_object_{dataset}_preprocessed_filtered_normalized_03.rds'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/process.R'

rule harmony:
    input:
        expand('objects/seurat_object_{dataset}_preprocessed_filtered_normalized_03.rds', dataset=datasets)
    output:
        'objects/seurat_object_harmony_integrated_04.rds'
    threads:
        num_workers * 4
    shell:
        '''
        source /data/abbass2/mambaforge/bin/activate harmony
        Rscript scripts/main/harmony.R {input} {output} {threads}
        '''
        
rule neighbors:
    input:
        seurat_object='objects/seurat_object_harmony_integrated_04.rds'
    output:
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_05.rds'
    conda:
        'envs/multiome.yml'
    script:
        'scripts/main/find_neighbors.R' 

rule umap:
    input:
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_05.rds'
    output:
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_umap_06.rds'
    conda:
        'envs/multiome.yml'
    script:
        'scripts/main/umap.R'

rule cluster:
    input:
        markers='input/major_cell_type_markers_list.rds',
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_umap_06.rds'
    output:
        plot='plots/major_type_module_umap.pdf',
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_umap_cluster_07.rds'
    conda:
        'envs/multiome.yml'
    params:
        algorithm=clustering_algorithm,
        resolution=clustering_resolution
    threads:
        num_workers * 4
    script:
        'scripts/main/clustering.R'

rule sctype:
    input:
        markers='input/major_cell_type_markers_list.rds',
        seurat_object='objects/seurat_object_harmony_integrated_neighbors_umap_cluster_07.rds'
    output:
        metadata='output/final_metadata.csv'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers * 4
    script:
        'scripts/main/annotate_clusters.R'

rule groups:
    input:
        metadata='output/final_metadata.csv'
    output:
        plot='plots/{group}_group_umap.pdf'
    params:
        groups='{group}'
    conda:
        'envs/multiome.yml'
    script:
        'scripts/main/plot_groups.R'

rule features:
    input:
        metadata='output/final_metadata.csv'
    output:
        plot='plots/{feature}_feature_umap.pdf'
    params:
        features='{feature}'
    conda:
        'envs/multiome.yml'
    script:
        'scripts/main/plot_features.R'


# rule markers:
#     input:
#         'objects/seurat_object_harmony_integrated_neighbors_umap_cluster_07.rds'
#     output:
#         markers='output/seurat_wnn_integrated_cluster_markers.csv',
#         umap='plots/seurat_wnn_seurat_clusters_sample_batch_umap.pdf'
#     conda:
#         'envs/multiome.yml'
#     threads:
#         num_workers * 2
#     script:
#         'scripts/main/markers.R'
