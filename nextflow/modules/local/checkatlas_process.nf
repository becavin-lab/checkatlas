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
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process SUMMARY{
    label 'process_summary'

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info

    script:
    out_info = atlas_info.atlas_name + "_Summary"
    """
    checkatlas summary $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --plot_celllimit ${params.plot_celllimit} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process QC{
    label 'process_qc'

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info

    script:
    out_info = atlas_info.atlas_name + "_QC"
    """
    checkatlas qc $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --qc_display ${params.qc_display} \
        --plot_celllimit ${params.plot_celllimit} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_CLUST{
    label 'process_metric_clust'

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info

    script:
    out_info = atlas_info.atlas_name + "_Metric_Clust"
    """
    checkatlas metric_cluster $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --obs_cluster ${params.obs_cluster} \
        --metric_cluster ${params.metric_cluster} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_ANNOT{
    label 'process_metric_annot'

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info

    script:
    out_info = atlas_info.atlas_name + "_Metric_Annot"
    """
    checkatlas metric_annot $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --obs_cluster ${params.obs_cluster} \
        --metric_annot ${params.metric_annot} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}

process METRIC_DIMRED{
    label 'process_metric_dimred'

    input:
    val atlas_info
    path samplesheet

    output:
    val out_info, emit: out_info

    script:
    out_info = atlas_info.atlas_name + "_Metric_dimred"
    """
    checkatlas metric_dimred $samplesheet --atlas_name ${atlas_info.atlas_name} \
        --metric_dimred ${params.metric_dimred} \
        ${params.checkatlas_debug ? '--debug' : ''}
    """
}
