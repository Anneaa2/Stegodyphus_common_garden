#!/bin/bash
#SBATCH --mem 8G
#SBATCH -c 1
#SBATCH --time=0-11:00:00
#SBATCH --account=spider2



#directory=$1
#file=$2
file=$1

extension=${file##*.}

#sort -s -k 11,11 $directory/$file | awk '
sort -s -k 11,11 $file | awk '
		 BEGIN{FS="\t";OFS="\t";m=0;u=0;count=0;getline;gene=$11;m=m+$5;u=u+$6;count++}
		 $11!=gene{print gene,m,u,count;m=0;u=0;gene=$11;count=0;m=m+$5;u=u+$6;count++;next}
		 {m=m+$5;u=u+$6;count++}
		 END{print gene,m,u,count;m=0;u=0;count=0}' > ${file/.$extension/.$extension.sum}		
		 #END{print gene,m,u,count;m=0;u=0;count=0}' > $directory/${file/.$extension/.$extension.sum}		