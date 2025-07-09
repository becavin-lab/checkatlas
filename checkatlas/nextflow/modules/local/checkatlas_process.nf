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
    checkatlas summary ${atlas_info.atlas_name} $samplesheet
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
    checkatlas qc ${atlas_info.atlas_name} $samplesheet
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
    checkatlas metric_cluster ${atlas_info.atlas_name} $samplesheet
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
    checkatlas metric_annot ${atlas_info.atlas_name} $samplesheet
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
    checkatlas metric_dimred ${atlas_info.atlas_name} $samplesheet
    """
}

