#!/bin/bash
#SBATCH --mem 32G
#SBATCH -c 1
#SBATCH --time=2-00:00:00
#SBATCH --account=spider2




source /home/anneaa/miniconda3/etc/profile.d/conda.sh
conda activate assay_copy

directory_is=$1
file_is=$2

echo $directory_is " and " $file_is

chmod a+x "assay_new_setup/derived_data_and_code/methylation/dss_run.r"
Rscript "assay_new_setup/derived_data_and_code/methylation/dss_run.r" $directory_is $file_is

