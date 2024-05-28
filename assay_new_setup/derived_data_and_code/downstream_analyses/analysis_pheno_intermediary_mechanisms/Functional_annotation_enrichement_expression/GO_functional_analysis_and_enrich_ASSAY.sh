


# this is the master script for running the funcitonal and go enrichement on Sources of Variation data



############## Load data and get genes that are strongly correlated. grep the fasta sequence from the genome.

# This functional annotation of all Stegodyphus dumicola files was done in the Sources_of_variation project.
	# paper: https://onlinelibrary.wiley.com/doi/10.1111/mec.16696
# Here I will reuse the same annotation in the ASSAY project (also on Stegodyphus dumicola)
	# Below is a decribtion (script) of how it was done

############################################################################

## IN BASH
cd "faststorage/project/spider2/AAL_20190722_SOV_correlations/"
# reference anotation of genome:  "/faststorage/project/spider2/AAL_20190722_SOV_correlations/january_run/snp_smp_data/genome_ano/v2_stringtie_good.gff"
genome="../dumicola_genome/assemblies/S_dumicola_genome.fasta"


##### make fasta with ALL genes and gene names in fasta header
### Make fasta with only gene seqs in fasta:
genome="../dumicola_genome/assemblies/S_dumicola_genome.fasta"
grep "transcript" < "/january_run/snp_smp_data/genome_ano/v2_stringtie_good.gff" > "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_transcripts_temp.gff"
grep -v "gp" < "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_transcripts_temp.gff" > "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_temp.gff"
	
rm "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/test_awk.gff"
awk 'BEGIN{OFS="\t";FS=" "}
	{split($9, array, "="); gene_id=array[3]; gene_pos1=$1; gene_pos2=$4; gene_pos3=$5;
	print gene_pos1,$2,gene_id, gene_pos2,gene_pos3,$6,$7,$8,$9 >> "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/test_awk.gff"
	}' "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_temp.gff"
# 37619 after uniq geneID
# 39298 before uniq --> thus some with several transcripts mapping to same gene.
bedinput="functional_analysis/gene_names_strong_correlations_pc/gene_fastas/test_awk.gff"
bedtools getfasta -name -fi $genome -bed $bedinput -fo "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all.fsa" 

rm gene_fastas/*temp*


## split into smaller pieces to run annotation faster:

seqtk seq "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all.fsa"  > "functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single.fsa" #change from multiline fasta to single line.
grep ">" functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single.fsa |wc -l
seqkit rmdup -s < functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single.fsa > functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa
grep ">" functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa |wc -l
cd functional_analysis/gene_names_strong_correlations_pc/gene_fastas
pyfasta split -n 100 functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa
cd "faststorage/project/spider2/AAL_20190722_SOV_correlations/"

		# retrieving the position of the gene in the genome using the genome annotation gff
		# using samtools to grep the fasta sequence from the genome
		# concatenate all fastas for each correlation group into one fasta file for functional annotation


############### Run functional annotation using Eggnog mapper

conda activate eggnog

input_file="functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all.fsa"
		### gene/header count = 39298
input_file="functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa"
		### gene/header count after removing doubles = 38321
					## This file was made in the Sources of Variation study, containing all positions of genes in the Dumicola genome.
###  RUn for ALL annotated genes
#sbatch run_func_annotation.sh $input_file $outfile_folder "all_genes"


# run all single annotations (but on split file):
#cd /faststorage/project/spider2/AAL_20190722_SOV_correlations/functional_analysis
outfile_folder="functional_analysis/func_anno_output"
for file in  functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa.[0-9][0-9]
do
	echo $file
	sbatch functional_analysis/run_func_annotation_quick.sh $file $outfile_folder $file"_all_genes_split"
done

# apparently it didn't get the outfolder right.
mv functional_analysis/gene_names_strong_correlations_pc/gene_fastas/*_all_genes_split* functional_analysis/func_anno_output/

# now cat the files together
	#skipping first 5 lines (except in first file), as well as last 4 lines.
#cd /faststorage/project/spider2/AAL_20190722_SOV_correlations/functional_analysis/func_anno_output/
	# first cat all files
cat functional_analysis/func_anno_output/*_all_genes_split*annotations > functional_analysis/func_anno_output/concat_all_genes_split_annotations_all
		### we loose some genes/lines in functional annotation (down to 36699) (some of these are header lines, to be removed below, so number is smaller)
		
cat	functional_analysis/func_anno_output/*_all_genes_split*hits > functional_analysis/func_anno_output/concat_all_genes_split_hits_all
cat	functional_analysis/func_anno_output/*_all_genes_split*seed_orthologs > functional_analysis/func_anno_output/concat_all_genes_split_seed_orthologs_all
	# now remove lines with # starting from line 6
sed '6,${/^#/d;}' functional_analysis/func_anno_output/concat_all_genes_split_annotations_all > functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_rem_temp # skips first five lines, then look for lines/strings starting with "#", the "d" refers to sed -delete
	# NO LINES TO REMOVE: sed '6,${/^#/d;}' concat_all_genes_split_hits_all > concat_all_genes_split_hits_all_rem # skips first five lines, then look for lines/strings starting with "#", the "d" refers to sed -delete
sed '7,${/^#/d;}' functional_analysis/func_anno_output/concat_all_genes_split_seed_orthologs_all > functional_analysis/func_anno_output/concat_all_genes_split_seed_orthologs_all_rem_temp # skips first six lines (delete from line 7), then look for lines/strings starting with "#", the "d" refers to sed -delete


# simplify file before input to R
	# only keeping GO and KEGG?
awk 'FS="\t", OFS="\t" {print $1,$2,$3,$10,$11,$12,$13,$14,$15,$16,$18}' functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_rem_temp > functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_reduced_temp
sed '1,4{/^#/d}' functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_reduced_temp > functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_reduced_R # among first 4 lines, delete occurances of lines starting with #
		### after removing header/end lines count is = 35908 --> so we loose about 2000 genes that did not get a functional annotation (2413)
	# this file can be entered into R
####
#Do similarly for 

# now the combined annotation file should be good.

# when concatenation is done, delete all individual mapping files
rm functional_analysis/func_anno_output/gene_positions_all_single_unique.fsa.[0-9][0-9]_all_genes_split.emapper*
rm functional_analysis/func_anno_output/*temp
# also remove split gene fastas
rm functional_analysis/gene_names_strong_correlations_pc/gene_fastas/gene_positions_all_single_unique.fsa.[0-9][0-9]



############################################################################

# Here I will reuse the same annotation in the ASSAY project (also on Stegodyphus dumicola)

############################################################################


# employ this script to run the actual functional enrichment
"assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/func_enrichement_ASSAY.r"
# obs, cannot be run by source or rscript, as it lacks shebang line. add and it will run fine, I think.