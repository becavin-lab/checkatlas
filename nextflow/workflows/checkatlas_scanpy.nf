/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS_SCANPY } from '../modules/local/checkatlas_process'
include { SUMMARY } from '../modules/local/checkatlas_process'
include { QC } from '../modules/local/checkatlas_process'
include { METRIC_CLUST } from '../modules/local/checkatlas_process'
include { METRIC_ANNOT } from '../modules/local/checkatlas_process'
include { METRIC_DIMRED } from '../modules/local/checkatlas_process'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Checkatlas scanpy WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow CHECKATLAS_SCANPY{
    take:
    atlas_info // dict: atlas_info
    ch_search_path

    main:
    // Preprocess scanpy
    PREPROCESS_SCANPY(atlas_info, ch_search_path)

    // Per-metric channels: 'none' produces an empty channel so the process is skipped
    ch_metric_clust  = params.metric_cluster == 'none' ? Channel.empty() : Channel.from(params.metric_cluster.split(' '))
    ch_metric_annot  = params.metric_annot   == 'none' ? Channel.empty() : Channel.from(params.metric_annot.split(' '))
    ch_metric_dimred = params.metric_dimred  == 'none' ? Channel.empty() : Channel.from(params.metric_dimred.split(' '))

    // Combine each atlas with every metric so every (atlas, metric) pair runs in parallel
    atlas_metric_clust  = PREPROCESS_SCANPY.out.atlas_info.combine(ch_metric_clust)
    atlas_metric_annot  = PREPROCESS_SCANPY.out.atlas_info.combine(ch_metric_annot)
    atlas_metric_dimred = PREPROCESS_SCANPY.out.atlas_info.combine(ch_metric_dimred)

    // Run all checkatlas processes
    SUMMARY(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    QC(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    METRIC_CLUST(atlas_metric_clust, ch_search_path)
    METRIC_ANNOT(atlas_metric_annot, ch_search_path)
    METRIC_DIMRED(atlas_metric_dimred, ch_search_path)

    // Mix all out channels
    scanpy_out = SUMMARY.out.out_info
    scanpy_out = scanpy_out.mix(QC.out.out_info)
    scanpy_out = scanpy_out.mix(METRIC_CLUST.out.out_info)
    scanpy_out = scanpy_out.mix(METRIC_ANNOT.out.out_info)
    scanpy_out = scanpy_out.mix(METRIC_DIMRED.out.out_info)

    emit:
    scanpy_out
}