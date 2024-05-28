#!/bin/bash
#SBATCH --mem 16G
#SBATCH -c 1
#SBATCH --time=0-12:00:00
#SBATCH --account=spider2



directory=$1
file=$2

extension=${file##*.}

sort -s -k 13,13 $directory/$file | awk '
		 BEGIN{FS="\t";OFS="\t";m=0;u=0;count=0;getline;gene=$13;m=m+$5;u=u+$6;count++}
		 $13!=gene{print gene,m,u,count;m=0;u=0;gene=$13;count=0;m=m+$5;u=u+$6;count++;next}
		 {m=m+$5;u=u+$6;count++}
		 END{print gene,m,u,count;m=0;u=0;count=0}' > $directory/${file/.$extension/.$extension.sum}		