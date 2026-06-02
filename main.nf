#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf/checkatlas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/becavin-lab/checkatlas
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
// if (params.validate_params) {
//     validateParameters()
// }

// WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCheckatlas.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/nextflow/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/nextflow/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK } from '../subworkflows/local/input_check'

// Run checkatlas for different altas objects
include { CHECKATLAS_SCANPY } from './nextflow/workflows/checkatlas_scanpy'
include { CHECKATLAS_CELLRANGER } from './nextflow/workflows/checkatlas_cellranger'
include { CHECKATLAS_SEURAT } from './nextflow/workflows/checkatlas_seurat'
include { LIST_SCANPY_ATLASES                     } from './nextflow/workflows/checkatlas_listfiles'
include { LIST_CELLRANGER_ATLASES                     } from './nextflow/workflows/checkatlas_listfiles'
include { LIST_SEURAT_ATLASES                     } from './nextflow/workflows/checkatlas_listfiles'
include { CREATE_REPORT                     } from './nextflow/workflows/checkatlas_multiqc'
include { COPY_MULTIQC_REPORT                     } from './nextflow/workflows/checkatlas_multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from './nextflow/modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './nextflow/modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

def create_atlas_info(LinkedHashMap row) {
    // Function to get list of 
    // Atlas_name,Atlas_type,Atlas_extension,Atlas_path
    def meta = [:]
    meta.atlas_name = row.Atlas_name
    meta.atlas_type = row.Atlas_type
    meta.atlas_extension = row.Atlas_extension
    meta.atlas_path = row.Atlas_path
    
    return meta
}


workflow CHECKATLAS {

    ch_versions = Channel.empty()

    //
    // MODULE: Checkatlas
    // 
    // Create Channel from search path
    myFile = file(params.path)
    ch_search_path = Channel.value(myFile)
    checkatlas_workingdir = file(params.path+"/checkatlas_files/")
    
    // Manage Scanpy atlases
    LIST_SCANPY_ATLASES(ch_search_path)
    
    LIST_SCANPY_ATLASES.out.list_scanpy.splitCsv( header:true, sep:',' )
        .map { create_atlas_info(it) }
        .set { atlas_info_scanpy }
    
    
    CHECKATLAS_SCANPY(atlas_info_scanpy, ch_search_path)
    ch_versions = ch_versions.mix(LIST_SCANPY_ATLASES.out.versions)
    
    // Manage Cellranger atlases
    /* LIST_CELLRANGER_ATLASES(ch_search_path)
    LIST_CELLRANGER_ATLASES.out.list_cellranger.splitCsv( header:true, sep:',' )
        .map { create_atlas_info(it) }
        .set { atlas_info_cellranger }
    CHECKATLAS_CELLRANGER(atlas_info_cellranger, ch_search_path)
    ch_versions = ch_versions.mix(LIST_CELLRANGER_ATLASES.out.versions) */

    // Manage Seurat atlases
    /* LIST_SEURAT_ATLASES(ch_search_path)
    LIST_SEURAT_ATLASES.out.list_seurat.splitCsv( header:true, sep:',' )
        .map { create_atlas_info(it) }
        .set { atlas_info_seurat }
    CHECKATLAS_SEURAT(atlas_info_seurat, ch_search_path)
    ch_versions = ch_versions.mix(LIST_SEURAT_ATLASES.out.versions) */
    
    // Collect all output value of checkatlas processes
    atlases_out = CHECKATLAS_SCANPY.out.scanpy_out
    //atlases_out = atlases_out.mix(CHECKATLAS_CELLRANGER.out.cellranger_out, CHECKATLAS_SEURAT.out.seurat_out)
    atlases_out = atlases_out.collect()

    // Run HTML report creation for QC plots, UMAP and tSNE
    CREATE_REPORT(ch_search_path, atlases_out)

    // Collect all software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCheckatlas.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    
    methods_description    = WorkflowCheckatlas.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    
    MULTIQC (CREATE_REPORT.out.out_info, ch_search_path,
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    COPY_MULTIQC_REPORT(checkatlas_workingdir, MULTIQC.out.report, MULTIQC.out.data)

}

workflow {
    CHECKATLAS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
