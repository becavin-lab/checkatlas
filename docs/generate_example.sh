cd /data/analysis/data_becavin/checkatlas_test/tuto/

# Run checkatlas for all data folder
nextflow run nf-core-checkatlas -r dev --path data1/ -queue-size 20
nextflow run nf-core-checkatlas -r dev --path data2/ -queue-size 20
nextflow run nf-core-checkatlas -r dev --path data3/ -queue-size 20
nextflow run nf-core-checkatlas -r dev --path data4/ -queue-size 20
nextflow run nf-core-checkatlas -r dev --path data5/ -queue-size 20

# Copy Multiqc reports
chk_path="/home/becavin/checkatlas/"
examples_path=${chk_path}"docs/examples/"

mkdir CheckAtlas_example_1
rm -rf CheckAtlas_example_1/*
cp -R data1/checkatlas_files/*.html CheckAtlas_example_1/
cp -R data1/checkatlas_files/multiqc_data/ CheckAtlas_example_1/
mv CheckAtlas_example_1/ ${examples_path}/

mkdir CheckAtlas_example_2
rm -rf CheckAtlas_example_2/*
cp -R data2/checkatlas_files/*.html CheckAtlas_example_2/
cp -R data2/checkatlas_files/multiqc_data/ CheckAtlas_example_2/
mv CheckAtlas_example_2/ ${examples_path}/

mkdir CheckAtlas_example_3
rm -rf CheckAtlas_example_3/*
cp -R data3/checkatlas_files/*.html CheckAtlas_example_3/
cp -R data3/checkatlas_files/multiqc_data/ CheckAtlas_example_3/
mv CheckAtlas_example_3/ ${examples_path}/

mkdir CheckAtlas_example_4
rm -rf CheckAtlas_example_4/*
cp -R data4/checkatlas_files/*.html CheckAtlas_example_4/
cp -R data4/checkatlas_files/multiqc_data/ CheckAtlas_example_4/
mv CheckAtlas_example_4/ ${examples_path}/

mkdir CheckAtlas_example_5
rm -rf CheckAtlas_example_5/*
cp -R data5/checkatlas_files/*.html CheckAtlas_example_5/
cp -R data5/checkatlas_files/multiqc_data/ CheckAtlas_example_5/
mv CheckAtlas_example_5/ ${examples_path}/

# Copy jupyter lab files
mkdir jupyter/
cp *.ipynb jupyter/ 
mv jupyter/ /home/becavin/checkatlas/docs/examples/