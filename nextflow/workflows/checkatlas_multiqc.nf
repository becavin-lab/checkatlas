


process CREATE_REPORT {
    debug true

    input:
    val ch_search_path
    val scanpy_out

    output:
    val out_info, emit: out_info

    script:
    out_info = "HTML reports"
    """
    checkatlas-workflow html_report $ch_search_path
    """
    
}

process COPY_MULTIQC_REPORT{

    input:
    val checkatlas_workingdir
    path report
    path data

    script:
    """
    cp $report ${checkatlas_workingdir}/Checkatlas_MultiQC.html
    if [ ! -d ${checkatlas_workingdir}/multiqc_data/ ]; then
	    mkdir ${checkatlas_workingdir}/multiqc_data/
    fi
    ls ${data}
    cp -R ${data}/* ${checkatlas_workingdir}/multiqc_data
    """

}

