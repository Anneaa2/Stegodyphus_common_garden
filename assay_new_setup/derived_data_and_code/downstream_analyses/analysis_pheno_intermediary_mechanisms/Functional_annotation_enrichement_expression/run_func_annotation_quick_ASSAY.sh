#!/bin/bash -l
#SBATCH --mem 64G
#SBATCH --time=2-00:00:00
#SBATCH --get-user-env
#SBATCH --account spider2

## This script downloads the new nucleotide database fasta file (nt.fsa) from NCBI in a gzipped format, unzips it and builts a local database.
## Run it like this: sbatch /faststorage/project/spider2/blast_db/Update_blast_db.sh OUTPUT_FOLDER
# eg. sbatch /faststorage/project/spider2/blast_db/Update_blast_db.sh /faststorage/project/spider2/blast_db/nt_db

conda activate eggnog

# unhash if you want to rebuild databases.
#download_eggnog_data.py -y

input_fasta=$1	# including path
output_folder=$2 #"/faststorage/project/spider2/blast_db/nt_db/"
run_name=$3

echo "Running functional annotation."

emapper.py -m diamond --dbmem --itype CDS -i $input_fasta -o $run_name --output_dir $output_folder



#--index_chunks 
#--block_size
#check diamonds page