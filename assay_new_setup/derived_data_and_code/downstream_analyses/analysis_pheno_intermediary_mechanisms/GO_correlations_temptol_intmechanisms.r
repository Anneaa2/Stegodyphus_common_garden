# Script to correlate phenotypes (heat tol, cold tol) with intermediary mechanisms
# written by Anne Aagaard 2022

#srun --mem=8g --pty /bin/bash
conda activate assay_copy ;R

#import libraries
#library(Hmisc) # for negation %nin% function
library(scales) # for alpha function
library(emmeans)
source("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/functions_for_pattern_finding.r")
#############################################
############ 		get data
#############################################

## heat tolerance
	ToD_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/ToD/ToD_nest-temp_mean.csv", header=T,stringsAsFactors=T) 
		# /faststorage/project/spider2/assay_study_anne/phenotypes/data/ToD/ToD_nest-temp_mean.csv
		#full_mod <- lm(ToD_mean ~ Population*Treatment + Mass_mean + Fed_mean + D_21_mean , data=ToD_final)
		#selec_mod <- step(full_mod)
	ToD_final$poptemp <-paste0(substr(as.character(ToD_final$Population),start=1,stop=1),ToD_final$Treatment)
	ToD_final_mean <- aggregate((ToD_final[1:nrow(ToD_final),"ToD_mean"]), by=list(ToD_final$poptemp), mean,na.rm=T)

		#get slopes for ctmax to see direction
	Betta_slope <- lm_slope(ToD_final[which(ToD_final$Population == "Betta"),"ToD_mean"], ToD_final[which(ToD_final$Population == "Betta"),"Treatment"])
	Otavi_slope <- lm_slope(ToD_final[which(ToD_final$Population == "Otavi"),"ToD_mean"], ToD_final[which(ToD_final$Population == "Otavi"),"Treatment"])
	Karasburg_slope <- lm_slope(ToD_final[which(ToD_final$Population == "Karasburg"),"ToD_mean"], ToD_final[which(ToD_final$Population == "Karasburg"),"Treatment"])
	Stampriet_slope <- lm_slope(ToD_final[which(ToD_final$Population == "Stampriet"),"ToD_mean"], ToD_final[which(ToD_final$Population == "Stampriet"),"Treatment"])
	slope_vec_CTmax <- c(Betta_slope,Otavi_slope,Karasburg_slope,Stampriet_slope)	; names(slope_vec_CTmax) <- c("Betta_slope","Otavi_slope","Karasburg_slope","Stampriet_slope")
	slope_vec_CTmax <- slope_vec_CTmax[order(slope_vec_CTmax,decreasing=T)]
		# get means for ctmax to see direction
	Betta_mean <- mean(ToD_final[which(ToD_final$Population == "Betta"),"ToD_mean"],na.rm=T)
	Otavi_mean <- mean(ToD_final[which(ToD_final$Population == "Otavi"),"ToD_mean"],na.rm=T)
	Karasburg_mean <- mean(ToD_final[which(ToD_final$Population == "Karasburg"),"ToD_mean"],na.rm=T)
	Stampriet_mean <- mean(ToD_final[which(ToD_final$Population == "Stampriet"),"ToD_mean"],na.rm=T)
	mean_vec_CTmax <- c(Betta_mean,Otavi_mean,Karasburg_mean,Stampriet_mean)	; names(mean_vec_CTmax) <- c("Betta_mean","Otavi_mean","Karasburg_mean","Stampriet_mean")
	mean_vec_CTmax <- mean_vec_CTmax[order(mean_vec_CTmax,decreasing=T)]
	
		# Get whether slopes were different from 0 in ctmax:
	ctmax_trends_0diff <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_trends_different_from_zero.csv",row.names=1)
	ctmax_zero_pops <- ctmax_trends_0diff[ which(ctmax_trends_0diff$p.value>0.05), "Population"]


	CCR_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/CCR/CCR_nest-temp_mean.csv", header=T,stringsAsFactors=T) 
		# /faststorage/project/spider2/assay_study_anne/phenotypes/data/CCR/CCR_nest-temp_mean.csv
		#full_mod <- lm(CCR_mean ~ Population*Treatment + Mass_mean + Fed + D_21 , data=CCR_final)
		#selec_mod <- step(full_mod)
	CCR_final$poptemp <-paste0(substr(as.character(CCR_final$Population),start=1,stop=1),CCR_final$Treatment)
	CCR_final_mean <- aggregate((CCR_final[1:nrow(CCR_final),"CCR_mean"]), by=list(CCR_final$poptemp), mean,na.rm=T)

	Betta_slope <- lm_slope(CCR_final[which(CCR_final$Population == "Betta"),"CCR_mean"], CCR_final[which(CCR_final$Population == "Betta"),"Treatment"])
	Otavi_slope <- lm_slope(CCR_final[which(CCR_final$Population == "Otavi"),"CCR_mean"], CCR_final[which(CCR_final$Population == "Otavi"),"Treatment"])
	Karasburg_slope <- lm_slope(CCR_final[which(CCR_final$Population == "Karasburg"),"CCR_mean"], CCR_final[which(CCR_final$Population == "Karasburg"),"Treatment"])
	Stampriet_slope <- lm_slope(CCR_final[which(CCR_final$Population == "Stampriet"),"CCR_mean"], CCR_final[which(CCR_final$Population == "Stampriet"),"Treatment"])
	slope_vec_CCR <- c(Betta_slope,Otavi_slope,Karasburg_slope,Stampriet_slope)	; names(slope_vec_CCR) <- c("Betta_slope","Otavi_slope","Karasburg_slope","Stampriet_slope")
	slope_vec_CCR <- slope_vec_CCR[order(slope_vec_CCR,decreasing=T)]

	Betta_mean <- mean(CCR_final[which(CCR_final$Population == "Betta"),"CCR_mean"],na.rm=T)
	Otavi_mean <- mean(CCR_final[which(CCR_final$Population == "Otavi"),"CCR_mean"],na.rm=T)
	Karasburg_mean <- mean(CCR_final[which(CCR_final$Population == "Karasburg"),"CCR_mean"],na.rm=T)
	Stampriet_mean <- mean(CCR_final[which(CCR_final$Population == "Stampriet"),"CCR_mean"],na.rm=T)
	mean_vec_CCR <- c(Betta_mean,Otavi_mean,Karasburg_mean,Stampriet_mean)	; names(mean_vec_CCR) <- c("Betta_mean","Otavi_mean","Karasburg_mean","Stampriet_mean")
	mean_vec_CCR <- mean_vec_CCR[order(mean_vec_CCR,decreasing=T)]
	
	ccr_trends_0diff <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_trends_different_from_zero.csv",row.names=1)
	ccr_zero_pops <- ccr_trends_0diff[ which(ccr_trends_0diff$p.value>0.05), "Population"]

##############################################################################
############ 		Correlation heat tol vs. plasticity in EXPRESSION
##############################################################################

## expression - genes with which effect
		#for expression data (also fdr - why two different cutoffs) earlier: 0.05
	padj_threshold = 0.1
	padj_threshold = 0.05
		# import expression data:
	exp_file = "assay_new_setup/derived_data_and_code/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt"
		#"/faststorage/project/spider2/assay_study_anne/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt"
	exp_dat <- read.table(exp_file, header=T,       row.names=1, stringsAsFactors=T)
	exp_dat_Temp <- exp_dat[	which( exp_dat[,which( colnames(exp_dat) == "padj_T" )] < padj_threshold	)	,	]
		# only working with those genes that show a plastic signal (temperature response)
	write.table(exp_dat_Temp$region_name, paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/functional_enrichement/gene_names_w_temperature_effect_",padj_threshold,".tab"),row.names=F,col.names=F)
		# export this data to make functional annotation and enrichement for genes with similar patterns as heat and cold tolerances
# expression counts
	## OBS: use vst counts, since they are transformed to have stabilized variance, evening out the variance inflation caused by really low or high count genes.
		# Deseq2 recommends it for clustering or visualization purposes. The scale is log2
		# from deseq2: The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance 
				# 		on the mean, particularly the high variance of the logarithm of count data when the mean is low Both VST 
				#		and rlog use the experiment-wide trend of variance over mean, in order to transform the data to remove the 
				#		experiment-wide trend. Note that we do not require or desire that all the genes have exactly the same 
				#		variance after transformation. Indeed, in a figure below, you will see that after the transformations the 
				#		genes with the same mean do not have exactly the same standard deviations, but that the experiment-wide 
				#		trend has flattened. It is those genes with row variance above the trend which will allow us to cluster 
				#		samples into interesting groups.
	counts_exp_file <- "assay_new_setup/derived_data_and_code/expression/results/data_for_downstream/norm_counts_vst_scaled_2022.csv"
		# "/faststorage/project/spider2/assay_study_anne/expression/results/data_for_downstream/norm_counts_vst_scaled_2022.csv"
		#counts_exp_file <- "/faststorage/project/spider2/assay_study_anne/expression/results/data_for_downstream/norm_counts_2022.csv"
	counts_exp <- read.csv(file=counts_exp_file, header=T,stringsAsFactors=T) 
	nrow(counts_exp) == nrow(exp_dat)
		## OBS IF FALSE
		# get only plastic genes counts
	plastic_counts_names <- intersect(exp_dat_Temp[,1],counts_exp[,1])
	plastic_counts <- counts_exp[counts_exp[,1] %in% plastic_counts_names,]
		# meta data for expression: 
	meta_exp <- read.table(file="assay_new_setup/derived_data_and_code/expression/analyzing/RNAextradction_samples_results022618_dosimp.txt",sep="\t", row.names=1,header=T)
	
		#/faststorage/project/spider2/assay_study_anne/expression/analyzing/RNAextradction_samples_results022618_dosimp.txt
	meta_exp$pop_temp <- gsub("N|E","O",meta_exp$pop_temp)
	meta_exp$population <- gsub("N|E","O",meta_exp$population)

		# order expression files to be pop/temp group means.
	rownames(plastic_counts) <- plastic_counts[,1]	; plastic_counts <- plastic_counts[,-1]
	plastic_counts_mean <- aggregate(t(plastic_counts[1:nrow(plastic_counts),]), by=list(meta_exp$pop_temp), mean)
	rownames(plastic_counts_mean) <- plastic_counts_mean[,1]		; plastic_counts_mean <- plastic_counts_mean[,-1]
	#plastic_counts_mean <- data.frame(t(plastic_counts_mean))


source("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/functions_for_pattern_finding.r")

# time consuming script
source("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/correlations_temptol_intmechanisms_expression.r")

' For padj 0.05:
[1] "270 genes fit the improved CTmax pattern (Nov22), out of 10435 genes."
[1] "65 genes fit the improved CTmax pattern (Nov22), out of 10435 genes."
'
' For padj 0.1:
[1] "274 genes fit the improved CTmax pattern (Nov22), out of 11807 genes."
[1] "71 genes fit the improved CCRTemp pattern (Nov22), out of 11807 genes."
'


##############################################################################
############ 		Correlation heat tol vs. plasticity in METABOLITES 
##############################################################################
source("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/functions_for_pattern_finding.r")

source("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/correlations_temptol_intmechanisms_metabolites.r")

##############################################################################
############ 		Correlation heat tol vs. plasticity in MICROBIOME
##############################################################################

# It makes no sense to look here, since we have none that respond plastically to temperature treatment.

##############################################################################
############ 		Correlation heat tol vs. plasticity in METYLATION
##############################################################################
# I will look for whether population pattern is similar for slopes in temp tolerances and means in methylation.
library(tidyr)

count_meth <- read.table(file="assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/wmlvl_allpops_d10D32_gene_cpg.txt",sep="\t", row.names=1,header=T)
colnames(count_meth) <- gsub("N|E","O",colnames(count_meth))
nrow(count_meth)
# Significantly differentially methylated genes dependent on population (below fdr 0.1, and below p val 0.05)
signif_meth <- read.table("assay_new_setup/derived_data_and_code/methylation/smps/gene_region/DSS_output/DSS_result_addit_d10D32_gene_CPG_population_significant.txt", header=T)
# It makes no sense to look here, since we have none that respond plastically to temperature treatment.

#subset counts to only contain significant
count_meth <- count_meth[	which(rownames(count_meth) %in% signif_meth$chr)	,	]
nrow(count_meth)

