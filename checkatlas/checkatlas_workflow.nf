nextflow.enable.dsl=2

process summary {
    input:
      path config_file
    script:
      """
      checkatlas-workflow summary ${config_file}
      """
}

process qc {
    input:
      path config_file
    script:
      """
      checkatlas-workflow qc ${config_file}
      """
}

process metric_cluster {
    input:
      path config_file
    script:
      """
      checkatlas-workflow metric_cluster ${config_file}
      """
}

process metric_annot {
    input:
      path config_file
    script:
      """
      checkatlas-workflow metric_annot ${config_file}
      """
}

process metric_dimred {
    input:
      path config_file
    script:
      """
      checkatlas-workflow metric_dimred ${config_file}
      """


}


workflow {
    config_file = channel.fromPath(params.files)
    summary(config_file)
    qc(config_file)
    metric_cluster(config_file)
    metric_annot(config_file)
    metric_dimred(config_file)
}