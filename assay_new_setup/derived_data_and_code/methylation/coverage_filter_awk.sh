#!/bin/bash
#SBATCH --mem 10G
#SBATCH -c 1
#SBATCH --time=0-12:00:00
#SBATCH --account spider2

infile=$1
outfile=$2
lowcut=$3

awk -v lowcut=$lowcut '(($5+$6)>=lowcut)&&(($5+$6)<=32){print}' $infile > $outfile


