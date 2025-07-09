/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUMMARY } from '../../modules/local/checkatlas_process'
include { QC } from '../../modules/local/checkatlas_process'
include { METRIC_CLUST } from '../../modules/local/checkatlas_process'
include { METRIC_ANNOT } from '../../modules/local/checkatlas_process'
include { METRIC_DIMRED } from '../../modules/local/checkatlas_process'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN Checkatlas seurat WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHECKATLAS_SEURAT{
    take:
    atlas_info // atlas_info : dict
    ch_search_path
    
    main:
    // Run all checkatlas processes
    SUMMARY(atlas_info, ch_search_path)
    QC(atlas_info, ch_search_path)
    METRIC_CLUST(atlas_info, ch_search_path)
    METRIC_ANNOT(atlas_info, ch_search_path)
    METRIC_DIMRED(atlas_info, ch_search_path)

    // Mix all out channels
    seurat_out = SUMMARY.out.out_info
    seurat_out = seurat_out.mix(QC.out.out_info, METRIC_CLUST.out.out_info)
    seurat_out = seurat_out.mix(METRIC_ANNOT.out.out_info, METRIC_DIMRED.out.out_info)
    
    emit:
    seurat_out
}