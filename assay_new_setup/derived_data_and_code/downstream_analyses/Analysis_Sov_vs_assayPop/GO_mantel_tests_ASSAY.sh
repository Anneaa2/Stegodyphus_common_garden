## ---------------- Make Fst files for both assay and SOV reduced

# done using this script:
assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/Fst_calc_gene_region_ASSAY_SOV.r


## ----------------- Run correlations for all datasets (Assay 15, 19, 23, 25, 29, and SOV_wild)



working_script_directory="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop"
env_dir="../AAL_20190722_SOV_correlations/environmental_dist_matrixes/reduce_env_june20_imp/"
runname="PCA"
directory_datafiles="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/WML_files"
snp_all_file="/faststorage/project/spider2/AAL_20190722_SOV_correlations/january_run/snp_smp_data/fst_snp/snp_genes_all_simple.fst"
outdirectory="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/"

mkdir -p $outdirectory
for SMPfile in `ls $directory_datafiles | grep "gene_pairwise_fst_smps_d10D3"`
do
    echo $directory_datafiles/$SMPfile
    sbatch $working_script_directory/r_pmantel_gene_smps2_all_refpop_ASSAY_SOV.sh $working_script_directory $env_dir $runname $directory_datafiles/$SMPfile $snp_all_file $outdirectory
done



# get top percent of null correlation (above 99.99 % of random correlation)
    # run it in R

# srun --mem=10g --pty /bin/bash
# conda activate sources_of_var; R
# source("figures_for_paper/pval_dist_improved_ASSAY.r")



# Run pval_dist_improved_ASSAY.r
    # Just by entering it. there is a bunch that likely will not run...
    # Only run top part.

# In sum, it did not reveal results that could be used. Maybe the analysis is still not the right one.

### ---------- Next idea ----------------##

# Run DSS using SOV data with 4 pops as time point 1,
    # ASSAY data for one temp (25) as time point 2.
    # Thus check differential methylation between those two at the gene level. 
    # I should have the files for that already.
    # Run this monday.










# Previously suggested analysis

## -------------------------------- Make files ready for input into DSS for modelling

# OBS DSS cannot be run on so few data points (it also makes great sense not fitting a line with 4 points)

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

#cd /faststorage/project/spider2/assay_study_anne/methylation/
chmod a+x make_input_for_dss.sh

directory_datafiles=../Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region

outdirectory=assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_dssrunSov_vs_assayPop/dss_files_SOV
working_script_directory=assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_dssrunSov_vs_assayPop

mkdir -p $outdirectory
for file in `ls $directory_datafiles | grep ".d10D30" | grep ".gene.cpg$"`
do
    echo $directory_datafiles/$file
    #echo $file
    sbatch $working_script_directory/make_input_for_dss_SOV.sh $directory_datafiles $file $outdirectory	
done
	# Outfile (_inputDSS.cpg) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads
	# Outfile (_inputDSS.cpg.counts) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads	count_of_Cs_in_region
	

# Run binominal model in R

#srun --mem=10g --pty /bin/bash
#conda activate assay; R


mv make_input_for_dss.sh*out outfiles/make_input_outs/

# ------------- R code to run DSS for methylation data per region


#chmod a+x dss_run.sh

#for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/*region`
#do 
	#for file in `ls $directory | grep ".d5D32_" | grep "_inputDSS.cpgf$"|head -n1`
	#do
	#	echo /faststorage/project/spider2/assay_study_anne/methylation/$directory/$file $file
	#	sbatch dss_run.sh /faststorage/project/spider2/assay_study_anne/methylation/$directory/ $file
	#done


# I cannot get dss to run on four populations with no replicates
    # It will not fit a linear model due to lacking degrees of freedom (which makes sence) --> thus I cannot make similar analyses across assay and SOV
    # Could I use a threshold for the pairwise values in distance matrix?
        # Such as: at least one population difference needs to be above .3
        # /home/anneaa/spider2/faststorage/AAL_20190722_SOV_correlations/correlation_results/genewise_SMPS_res10/dist_matrices
            #remove G and N rows
    # Or - I should take all correlations to climatic parameters (before I only took temperature related parameters),
    # and assume that these gradients make up a population signal. Then use top correlations such as in the SOV paper.
        # This means rerunning the correlation analyses without N and G
        # OR rerunning correlations that give a result in SOV, removing N and G
    # or: we turn it around and do SOV correlations to PCs for both assay (each temp) and SOV data subset.



for file in `ls $outdirectory | grep ".d10D30." | grep "_inputDSS.cpgf$"|head -n1`
do
    echo $outdirectory $file
    sbatch $working_script_directory/dss_run_SOV.sh $outdirectory $file	
done


mv dss_run.sh*out outfiles/dss_run_outs/
