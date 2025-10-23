#data_path=/data/analysis/data_moccia/checkatlas/tuto/data1
data_path=/data/analysis/data_becavin/checkatlas_test/tuto/data1

#checkatlas_path=/home/moccia/checkatlas
#checkatlas_path=/home/becavin/checkatlas


# nextflow
#nextflow run main.nf --path=/data/analysis/data_moccia/checkatlas/tuto/data1/ --outdir=/data/analysis/data_moccia/checkatlas/tuto/data1/

# metric (default)
checkatlas metric_cluster Fetal ${data_path} --debug
