#!/bin/bash
#SBATCH --mem 10G
#SBATCH -c 1
#SBATCH --time=0-02:00:00
#SBATCH --account spider2

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

#directory=$1
#file=$2
file=$1
#output_directory=$3
output_directory=`dirname $file`
output_directory=${output_directory/cov_files/DSS_input_files}
extension=${file##*.}
base_file=`basename $file`

 # file output: gene, start position, total nr. reads, nr. methylated reads, number of c sites in gene
sort -s -k 11,11 $file | awk '
		 BEGIN{FS="\t";OFS="\t";m=0;u=0;count=0;pos=0;getline;gene=$11;m=m+$5;u=u+$6;pos=$8;count++}
		 $11!=gene{print gene,pos,m+u,m,count;m=0;u=0;gene=$11;pos=0;count=0;m=m+$5;u=u+$6;pos=$8;count++;next}
		 {m=m+$5;u=u+$6;count++}
		 END{print gene,pos,m+u,m,count;m=0;u=0;pos=0;count=0}' > $output_directory/${base_file/.$extension/_inputDSS.$extension".count"}

# m = methylated , u = unmethylated, pos= position of starting cytosine (base, not cytosine)

	# Cut the count column awai for DSS
cut -d$'\t' -f 1,2,3,4 $output_directory/${base_file/.$extension/_inputDSS.$extension".count"} > $output_directory/${base_file/.$extension/_inputDSS.$extension"f"}
	# Add header for files
sed -i '1s/^/chr\tpos\tN\tX\n/' $output_directory/${base_file/.$extension/_inputDSS.$extension"f"}		
