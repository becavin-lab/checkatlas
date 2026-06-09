process PREPROCESS_SCANPY{
    label 'process_preprocess_scanpy'

    input:
    val atlas_info
    path samplesheet

    output:
    val atlas_info, emit: atlas_info

    script:
    out_info = atlas_info.atlas_name + "_PreProcess_Scanpy"
    """
    checkatlas preprocess $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process SUMMARY{
    label 'process_summary'
    tag "${atlas_info.atlas_name}"

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info
    path "${samplesheet}/checkatlas_files/summary/${atlas_info.atlas_name}.tsv",   emit: summary_tsv
    path "${samplesheet}/checkatlas_files/adata/${atlas_info.atlas_name}.tsv",     emit: adata_tsv
    path "${samplesheet}/checkatlas_files/umap/${atlas_info.atlas_name}_umap.png", optional: true, emit: umap_fig
    path "${samplesheet}/checkatlas_files/tsne/${atlas_info.atlas_name}_tsne.png", optional: true, emit: tsne_fig

    script:
    out_info = atlas_info.atlas_name + "_Summary"
    """
    checkatlas summary $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --plot_celllimit ${params.plot_celllimit} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process QC{
    label 'process_qc'
    tag "${atlas_info.atlas_name}"

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info
    path "${samplesheet}/checkatlas_files/qc/${atlas_info.atlas_name}.tsv",          emit: qc_tsv
    path "${samplesheet}/checkatlas_files/violin/${atlas_info.atlas_name}_qc.png",   optional: true, emit: qc_fig

    script:
    out_info = atlas_info.atlas_name + "_QC"
    """
    checkatlas qc $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --qc_display ${params.qc_display} \
        --plot_celllimit ${params.plot_celllimit} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_CLUST{
    label 'process_metric_clust'
    tag "${atlas_info.atlas_name}"

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info
    path "${samplesheet}/checkatlas_files/cluster/${atlas_info.atlas_name}.tsv", optional: true, emit: cluster_tsv

    script:
    out_info = atlas_info.atlas_name + "_Metric_Clust"
    """
    checkatlas metric_cluster $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --obs_cluster ${params.obs_cluster} \
        --metric_cluster ${params.metric_cluster} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_ANNOT{
    label 'process_metric_annot'
    tag "${atlas_info.atlas_name}"

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info
    path "${samplesheet}/checkatlas_files/annotation/${atlas_info.atlas_name}.tsv", optional: true, emit: annot_tsv

    script:
    out_info = atlas_info.atlas_name + "_Metric_Annot"
    """
    checkatlas metric_annot $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --obs_cluster ${params.obs_cluster} \
        --metric_annot ${params.metric_annot} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_DIMRED{
    label 'process_metric_dimred'
    tag "${atlas_info.atlas_name}"

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info
    path "${samplesheet}/checkatlas_files/dimred/${atlas_info.atlas_name}.tsv", optional: true, emit: dimred_tsv

    script:
    out_info = atlas_info.atlas_name + "_Metric_Dimred"
    """
    checkatlas metric_dimred $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --metric_dimred ${params.metric_dimred} \
        --n_jobs ${params.n_jobs} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}
