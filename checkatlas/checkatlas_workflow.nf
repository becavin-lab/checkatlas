nextflow.enable.dsl=2


process checkatlas_workflow {
    input:
      path config_file

    script:
      """
      echo ${config_file};
      checkatlas-workflow ${config_file}
      """
}

workflow {
    config_file = channel.fromPath(params.files)
    checkatlas_workflow(config_file)
}