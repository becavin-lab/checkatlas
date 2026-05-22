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
    // Run all checkatlas processes
    SUMMARY(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    QC(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    METRIC_CLUST(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    METRIC_ANNOT(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    METRIC_DIMRED(PREPROCESS_SCANPY.out.atlas_info, ch_search_path)
    
    // Mix all out channels
    scanpy_out = SUMMARY.out.out_info
    scanpy_out = scanpy_out.mix(QC.out.out_info, METRIC_CLUST.out.out_info)
    scanpy_out = scanpy_out.mix(METRIC_CLUST.out.out_info)
    scanpy_out = scanpy_out.mix(METRIC_ANNOT.out.out_info, METRIC_DIMRED.out.out_info)
    
    emit:
    scanpy_out
}