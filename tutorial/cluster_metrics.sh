#data_path=/data/analysis/data_moccia/checkatlas/tuto/data1
data_path=/data/analysis/data_becavin/checkatlas_test/tuto/data1

#checkatlas_path=/home/moccia/checkatlas
checkatlas_path=/home/becavin/checkatlas


# nextflow
nextflow run ${checkatlas_path}/main.nf --path=${data_path} --outdir=${data_path}

# metric (default)
checkatlas metric_cluster Fetal ${data_path} --debug

# run MultiQC only
multiqc --force ${data_path}