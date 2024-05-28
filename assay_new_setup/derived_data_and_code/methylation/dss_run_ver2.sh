#!/bin/bash
#SBATCH --mem 128G
#SBATCH -c 1
#SBATCH --time=2-10:00
#SBATCH --account spider2


source /home/anneaa/miniconda3/etc/profile.d/conda.sh
conda activate assay

directory_is=$1
file_is=$2

cd $directory_is
echo $directory_is " and " $file_is

chmod a+x "/faststorage/project/spider2/assay_study_anne/methylation/dss_run_ver2.r"
"/faststorage/project/spider2/assay_study_anne/methylation/dss_run_ver2.r" $directory_is $file_is