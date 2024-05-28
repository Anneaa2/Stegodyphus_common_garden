
# In this script I investigate the number of overlapping genes between Sources of variation and the common garden study.
# I only look at genes whoose methylation pattern correlates to temperature parameters.

conda activate assay_copy; R


#setwd("/faststorage/project/spider2/assay_study_anne/analysis_meth_SOV_candidates")
# import and subset candidate list from SOV:
candidate_file <- "/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/only_temp_prec/genewise_SMPS_res/temp_prec_partial_mantel_res_genes_smps_refpop_ALL.csv"

candidates <- read.csv(candidate_file, header=T, row.names=1)
NROW(unique(candidates$smp_vector))  
	# number of genes total in Sources of variation study, where methylation was found in some pops. ie. non-zero methylation levels.
	#[1] 31830

dim(candidates)
		#[1] 1623329       6
	# subset environmental vectors to only be the ones with temperature
candidates <- candidates[grep("Temp|T_", candidates$env_vector),]
dim(candidates)
		#[1] 1400519       6
		
# import file with correelation coefficients threshold:
thresholds <- read.csv("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/only_temp_prec/N_above_corcoef_tempprec_9999percentile_mantel.csv", header=T, row.names=1)
thresholds <- thresholds[grep("smp",rownames(thresholds)),]

# env_var=levels(candidates$env_vector)[1]
for( env_var in levels(candidates$env_vector) ) {
		# get threshold from thresholds dataframe:
	threshold_number <- thresholds[          grep(env_var, rownames(thresholds)),    which(colnames(thresholds) == "Perm_0.9999")    ]
		# subset candidates dataframe to only have the given env.
	sub_candidates <- candidates[grep(env_var, candidates$env_vector),]
		# subset sub candidates, to have only the ones passing threshold
	sub_candidates <- sub_candidates[which(sub_candidates$r_vector > threshold_number),]
		# make dataframe: env_vector, candidate_gene, 	r_vector, below_correlation_coeficient_of
	if ( nrow(sub_candidates) == 0 ) {	next	}
	if(env_var == levels(candidates$env_vector)[1] ){
		df_candidate <- data.frame(sub_candidates$env_vector, sub_candidates$smp_vector, sub_candidates$r_vector, threshold_number		)
	}else{
		df_new <- data.frame(sub_candidates$env_vector, sub_candidates$smp_vector, sub_candidates$r_vector, threshold_number		)
		df_candidate <- rbind(df_candidate, df_new)
	}	
}

colnames(df_candidate) <- c("environmental_parameter","gene","correlation_coefficient","cor_coeff_threshold(99.99th_percentile_of_null_dist)")

NROW(unique(df_candidate$gene))
	#  4582 		....genes show strong correlation to *any* temperature related parameters
write.table(df_candidate, file="assay_new_setup/derived_data_and_code/downstream_analyses/analysis_meth_SOV_candidates/Candidate_genes_dataframe.txt", sep="\t", row.names=F)
write.table(unique(df_candidate$gene), file="assay_new_setup/derived_data_and_code/downstream_analyses/analysis_meth_SOV_candidates/Candidate_genes_unique_list.txt", sep="\t", row.names=F)

candidate_vec <- gsub("gene_smp_","",as.character(unique(df_candidate$gene)))
	# These are genes with strong correlation to temperature parameters ( > 99.99th percentile)
	#  4582 		....genes show strong correlation to *any* temperature related parameters

##############################################################################################################
# Get genes in common garden with population response (since that is all there is):
##############################################################################################################

meth_dat <- read.table(file="assay_new_setup/derived_data_and_code/methylation/smps/gene_region/DSS_output/DSS_result_addit_d10D32_gene_CPG_population_significant.txt",sep="\t", row.names=1,header=T)
meth_dat_sub <- meth_dat[which(meth_dat$fdrs < .1),]
nrow(meth_dat_sub)
	# [1] 1537


meth_level_dat <- read.table(file="assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/wmlvl_allpops_d10D32_gene_cpg.txt",sep="\t", row.names=1,header=T)
nrow(meth_level_dat)	# number of genes, where we detected methylated reads in common garden study.
# [1] 34971					... as opposed to 31830 from sources of variation.


common_genes_SOV_Assay <- intersect(rownames(meth_dat_sub),candidate_vec)

NROW(common_genes_SOV_Assay)
write.table(common_genes_SOV_Assay,file="assay_new_setup/derived_data_and_code/downstream_analyses/analysis_meth_SOV_candidates/common_genes_genenames.txt", sep="\t", row.names=F,col.names=F)
#[1] 264
meth_level_dat_sub <- meth_level_dat[which(rownames(meth_level_dat) %in% common_genes_SOV_Assay),]
write.table(meth_level_dat_sub,file="assay_new_setup/derived_data_and_code/downstream_analyses/analysis_meth_SOV_candidates/common_genes_wmlvl_cpg_allpops.txt", sep="\t", row.names=T,col.names=T)


# 264/1537*100
# [1]	17.17632			... 17 % og genes with pop effect in methylation also had high correlation in SOV.

#264/34971*100
#[1] 0.7549112				.. .75 % of genes from common garden was also found in SOV.

# conclusion: 	1) methylation does not drive plastic response in temperaure acclimation (of subadults)
		# 		2) 264 / 1537 are possibly responsive to long term temperature stimuli based on candidate list from Sources_of_variation study.
		#				a) The 264 genes are not plastic in subadult stage (or they would have shown acclimation temperature response)
		#				b) May be developmentally plastic (and therefore show up as population response in our common garden)
		






