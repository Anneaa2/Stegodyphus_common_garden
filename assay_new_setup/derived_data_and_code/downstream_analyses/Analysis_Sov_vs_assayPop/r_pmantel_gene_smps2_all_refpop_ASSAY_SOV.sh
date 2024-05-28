#!/bin/bash
#SBATCH --mem 16G
#SBATCH -c 1
#SBATCH --time=0-11:00:00
#SBATCH --account spider2


directory_is=$1
echo "reading from directory" $directory_is

env_dir=$2
runname=$3
SMPS_GENE_FILE=$4 ## OBS
SNP_all_file=$5 ## OBS
outdir=$6

source /home/anneaa/miniconda3/etc/profile.d/conda.sh
#source ~/miniconda3/etc/profile.d/conda.sh
conda activate sources_of_var

#main="/genewise_smps_env/"
#results="/results"
#graph="/graphics"
#matrices="/smp_matrices"

#[ ! -d "$directory_is$main" ] && mkdir -p "$directory_is$main"
#[ ! -d "$directory_is$main$results" ] && mkdir -p "$directory_is$main$results"
#[ ! -d "$directory_is$main$graph" ] && mkdir -p "$directory_is$main$graph"
#[ ! -d "$directory_is$main$matrices" ] && mkdir -p "$directory_is$main$matrices"

chmod a+x "$directory_is/r_pmantel_gene_smps2_all_refpop_ASSAY_SOV.r"

"$directory_is/r_pmantel_gene_smps2_all_refpop_ASSAY_SOV.r" $directory_is $env_dir $runname $SMPS_GENE_FILE $SNP_all_file $outdir
