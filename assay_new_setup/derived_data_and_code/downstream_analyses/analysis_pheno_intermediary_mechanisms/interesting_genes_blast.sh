

# conda activate assay_copy ; cd /home/anneaa/spider2/faststorage/assay_study_anne/pheno_intermediary_mechanisms 

# transcriptome of Assay RNA: /faststorage/project/spider2/assay_study_anne/expression/assembling/merged/RNA_all_merged_ref.gtf


# reference anotation of genome:  "/faststorage/project/spider2/AAL_20190722_SOV_correlations/january_run/snp_smp_data/genome_ano/v2_stringtie_good.gff"

genome="/faststorage/project/spider2/dumicola_genome/assemblies/S_dumicola_genome.fasta"
ctmax_genes=`cut -f2 -d "," CTmax_Exp_subset_fits_patterns.csv | tail -n +2 | sed "s/\"//g"` # get only second column #remove first line
ccr_genes=` cut -f2 -d "," CCR_Exp_subset_fits_patterns.csv | tail -n +2 | sed "s/\"//g"` # get only second column #remove first line

# get posistion of intereesting gene:
for gene_id in $ctmax_genes
	#nake variable for gene_id_input
	do
	echo $gene_id
	genome_pos1=`grep -m 1 $gene_id < "/faststorage/project/spider2/assay_study_anne/expression/assembling/merged/RNA_all_merged_ref.gtf"|rev|cut -f2 -d"\""|rev` #;genome_pos1=\"$genome_pos1\"
	#echo $genome_pos1
	if [ $gene_id == $genome_pos1 ]
	then 
		echo "TRUE" 
	else 
		echo "FALSE" 
	fi
	genome_pos=`grep -m 1 $gene_id < "/faststorage/project/spider2/assay_study_anne/expression/assembling/merged/RNA_all_merged_ref.gtf" | cut -f1,4,5 ` # outputs SEQ_XXXX START_pos END_pos
	genome_pos=`echo $genome_pos | sed '0,/ /{s/ /:/}'`
	genome_pos=`echo $genome_pos | sed '0,/ /{s/ /-/}'`
	echo $genome_pos
	samtools faidx $genome $genome_pos > gene_seqs_specific_RNA/temp.fsa # pick out sequence and write to file as fasta
	sed "s/>.*/&:$gene_id/" gene_seqs_specific_RNA/temp.fsa > gene_seqs_specific_RNA/${gene_id/g/ctmax_nt_g}.fsa #modify fasta to add gene_id to header
done

# combine files to add to blast webpage as a multiple fasta file.
rm gene_seqs_specific_RNA/allSeqs_ctmax_nt*
counter=0 ; counter2=0
for file in `ls gene_seqs_specific_RNA/ctmax_nt_g*.fsa`
do
	counter=$((counter+1)) ; echo $file
	number_files=`ls gene_seqs_specific_RNA/ctmax_nt_g*.fsa|wc -l`
	if [[ $((counter%7)) -eq 1 ]] ; # if lefover from division by N equals 1
		then echo "yes"; 
		counter2=$((counter2+1))
		filename=gene_seqs_specific_RNA/allSeqs_ctmax_nt_$counter2.fsa
	fi
	cat $file >> $filename
done 

# do the same for chill coma recovery
for gene_id in $ccr_genes
	#nake variable for gene_id_input
	do
	echo $gene_id
	genome_pos1=`grep -m 1 $gene_id < "/faststorage/project/spider2/assay_study_anne/expression/assembling/merged/RNA_all_merged_ref.gtf"|rev|cut -f2 -d"\""|rev` #;genome_pos1=\"$genome_pos1\"
	#echo $genome_pos1
	if [ $gene_id == $genome_pos1 ]
	then 
		echo "TRUE" 
	else 
		echo "FALSE" 
	fi
	genome_pos=`grep -m 1 $gene_id < "/faststorage/project/spider2/assay_study_anne/expression/assembling/merged/RNA_all_merged_ref.gtf" | cut -f1,4,5 ` # outputs SEQ_XXXX START_pos END_pos
	genome_pos=`echo $genome_pos | sed '0,/ /{s/ /:/}'`
	genome_pos=`echo $genome_pos | sed '0,/ /{s/ /-/}'`
	echo $genome_pos
	samtools faidx $genome $genome_pos > gene_seqs_specific_RNA/temp.fsa
	sed "s/>.*/&:$gene_id/" gene_seqs_specific_RNA/temp.fsa > gene_seqs_specific_RNA/${gene_id/g/ccr_nt_g}.fsa
done


rm gene_seqs_specific_RNA/allSeqs_ccr_nt*
counter=0 ; counter2=0
for file in `ls gene_seqs_specific_RNA/ccr_nt_g*.fsa`
do
	counter=$((counter+1)) ; echo $file
	number_files=`ls gene_seqs_specific_RNA/ccr_nt_g*.fsa|wc -l`
	if [[ $((counter%7)) -eq 1 ]] ; # if lefover from division by N equals 1
		then echo "yes"; 
		counter2=$((counter2+1))
		filename=gene_seqs_specific_RNA/allSeqs_ccr_nt_$counter2.fsa
	fi
	cat $file >> $filename
done 


rm gene_seqs_specific_RNA/temp.fsa


#tried this, did not work.
	#/home/anneaa/projects/ncbi-blast-2.9.0+/bin/blastn -db nt -query nt.fsa -out test.out -remote
	# It also didn't work getting this command to work:
		#perl /home/anneaa/projects/ncbi-blast-2.9.0+/bin/update_blastdb.pl --showall ## could be modified to download formatted database. Cannot do that from cluster, apparently.
# Building local db from downloaded nt database using commands below. Should / could be worked into a script
	# cd /home/anneaa/spider2/faststorage/blast_db
	# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
	# gunzip nt.gz
	# makeblastdb -in nt.gz -out nt_db.gz -parse_seqids -dbtype nucl


# Now blast to local nt database and proceed with GO/network?




#RID_ID=`grep RID test.out`
#RID: X3R7GAUS014
#blast_formatter –rid $RID_ID –out test.tab –outfmt 7 
#blast_formatter –rid $RID_ID –out test.xml –outfmt 5 