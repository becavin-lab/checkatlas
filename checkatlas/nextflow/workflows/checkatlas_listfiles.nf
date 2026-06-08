
process LIST_SCANPY_ATLASES {
    debug true
    
    input:
    val ch_search_path

    output:
    path "List_scanpy.csv", emit: list_scanpy
    path "versions.yml", emit: versions

    script:
    """
    checkatlas-workflow list_scanpy $ch_search_path
    # copy List .csv to put them iin the scope of nextflow 
    cp ${ch_search_path}/checkatlas_files/List_scanpy.csv List_scanpy.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkatlas: \$(checkatlas --version | sed 's/Checkatlas, version //g')
    END_VERSIONS
    """
}

process LIST_CELLRANGER_ATLASES {
    debug true
    
    input:
    val ch_search_path

    output:
    path "List_cellranger.csv", emit: list_cellranger
    path "versions.yml", emit: versions

    script:
    """
    checkatlas-workflow list_cellranger $ch_search_path
    # copy List .csv to put them iin the scope of nextflow 
    cp ${ch_search_path}/checkatlas_files/List_cellranger.csv List_cellranger.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkatlas: \$(checkatlas --version | sed 's/Checkatlas, version //g')
    END_VERSIONS
    """
}

process LIST_SEURAT_ATLASES {
    debug true
    
    input:
    val ch_search_path

    output:
    path "List_seurat.csv", emit: list_seurat
    path "versions.yml", emit: versions

    script:
    """
    checkatlas-workflow list_seurat $ch_search_path
    # copy List .csv to put them iin the scope of nextflow 
    cp ${ch_search_path}/checkatlas_files/List_seurat.csv List_seurat.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkatlas: \$(checkatlas --version | sed 's/Checkatlas, version //g')
    END_VERSIONS
    """
}


