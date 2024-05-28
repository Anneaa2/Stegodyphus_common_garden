#!/bin/bash
#SBATCH --mem 16G
#SBATCH -c 1
#SBATCH --time=0-12:00:00
#SBATCH --account spider2


source /home/anneaa/miniconda3/etc/profile.d/conda.sh
conda activate assay

directory_is=$1
file_is=$2
extension_is=${file_is#*.}

chmod a+x "assay_new_setup/derived_data_and_code/methylation/r_calc_weighted_mlvl.r"
"assay_new_setup/derived_data_and_code/methylation/r_calc_weighted_mlvl.r" $directory_is $file_is $extension_is