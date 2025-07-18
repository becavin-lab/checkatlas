# nextflow
nextflow run main.nf --path=/data/analysis/data_moccia/checkatlas/tuto/data1/ --outdir=/data/analysis/data_moccia/checkatlas/tuto/data1/

# metric (default)
checkatlas metric_cluster Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug
checkatlas metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug
checkatlas metric_dimred Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug 


# metric

## FETAL

### DimRed

checkatlas metric_dimred Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_dimred=kruskal_stress # on teste les données fetal avec kruskal \


### Cluster 
checkatlas metric_cluster Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=silhouette # on teste les données fetal avec sillhouette

checkatlas metric_cluster Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=davies_bouldin  # on teste les données fetal avec davies-bouldin

checkatlas metric_cluster Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=calinski_harabasz  # on teste les données fetal avec calinski

### Annot
checkatlas metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=vmeasure # on teste les données fetal avec vmeasure

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=rand_index # on teste les données fetal avec rand

checkatlas metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=fowlkes_mallow # on teste les données fetal avec fowlkes-mallow

checkatlas metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=adj_mutual_info # on teste les données fetal avec AMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=adj_rand_index # on teste les données fetal avec ARI \

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=normalized_mutual_info  # on teste le données fetal avec NMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=mutual_info  # on teste le données fetal avec MI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Fetal /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=isolated_f1_score  # on teste le données fetal avec f1


## B-CELL

### DimRed
checkatlas metric_dimred B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_dimred=kruskal_stress # on teste les données b-cells avec kruskal

### Cluster 
checkatlas metric_cluster B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=silhouette # on teste les données b-cells avec sillhouette

checkatlas metric_cluster B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=davies_bouldin  # on teste les données b-cells avec davies-bouldin

checkatlas metric_cluster B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=calinski_harabasz  # on teste les données b-cells avec calinski

### Annot
checkatlas metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=vmeasure # on teste les données b-cells avec vmeasure

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=rand_index # on teste les données b-cells avec rand \ 

checkatlas metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=fowlkes_mallow # on teste les données b-cells avec fowlkes-mallows

checkatlas metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=adj_mutual_info # on teste les données b-cells avec AMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=adj_rand_index # on teste les données fetal avec ARI \

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=normalized_mutual_info  # on teste le données b-cell avec NMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=mutual_info  # on teste le données b-cell avec MI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot B-cells_compartment /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=isolated_f1_score  # on teste le données b-cell avec f1

## TABULA

### DimRed
checkatlas metric_dimred Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_dimred=kruskal_stress # on teste les données Tabula avec kruskal

### Cluster 
checkatlas metric_cluster Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=silhouette # on teste les données Tabula avec sillhouette

checkatlas metric_cluster Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=davies_bouldin  # on teste les données Tabula avec davies-bouldin

checkatlas metric_cluster Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_cluster=calinski_harabasz  # on teste les données Tabula avec calinski

### Annot
checkatlas metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=vmeasure # on teste les données Tabula avec vmeasure

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=rand_index # on teste les données Tabula avec rand

checkatlas metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=fowlkes_mallow # on teste les données Tabula avec fowlkes-mallows

checkatlas metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 --debug --metric_annot=adj_mutual_info # on teste les données Tabula avec AMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=adj_rand_index # on teste les données tabula avec ARI \

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=normalized_mutual_info  # on teste le données tabula avec NMI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=mutual_info  # on teste le données tabula avec MI

python /home/moccia/checkatlas/checkatlas/__main__.py metric_annot Tabula_Sapiens_Endothelial /data/analysis/data_moccia/checkatlas/tuto/data1 \
    --debug --metric_annot=isolated_f1_score  # on teste le données tabula avec f1