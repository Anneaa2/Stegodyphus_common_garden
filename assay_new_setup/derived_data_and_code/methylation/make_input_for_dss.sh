#!/bin/bash
#SBATCH --mem 10G
#SBATCH -c 1
#SBATCH --time=0-12:00:00
#SBATCH --account spider2

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

directory=$1
file=$2
extension=${file##*.}

 # file output: gene, start position, total nr. reads, nr. methylated reads, number of c sites in gene
sort -s -k 13,13 $directory/$file | awk '
		 BEGIN{FS="\t";OFS="\t";m=0;u=0;count=0;pos=0;getline;gene=$13;m=m+$5;u=u+$6;pos=$10;count++}
		 $13!=gene{print gene,pos,m+u,m,count;m=0;u=0;gene=$13;pos=0;count=0;m=m+$5;u=u+$6;pos=$10;count++;next}
		 {m=m+$5;u=u+$6;count++}
		 END{print gene,pos,m+u,m,count;m=0;u=0;pos=0;count=0}' > $directory/${file/.$extension/_inputDSS.$extension".count"}		

# m = methylated , u = unmethylated, pos= position of starting cytosine

	# Cut the count column awai for DSS
cut -d$'\t' -f 1,2,3,4 $directory/${file/.$extension/_inputDSS.$extension".count"}		> $directory/${file/.$extension/_inputDSS.$extension"f"}		

	# Add header for files
sed -i '1s/^/chr\tpos\tN\tX\n/' $directory/${file/.$extension/_inputDSS.$extension"f"}		
