### R code inspired from vignette source 'GOstatsForUnsupportedOrganisms.Rnw'

# srun --mem=16g --pty /bin/bash
conda activate assay_func_enrich ; R


###################################################
### code chunk number 1: Library import
###################################################
library("AnnotationForge")
library("GOstats")
library("GSEABase")

fdr_thresh=0.1

###################################################
### code chunk number 2: Read in genelists and functional annotation
###################################################

funct_annot_all_genes <- "assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/functional_annotation_S_dum/concat_all_genes_split_annotations_all_reduced_R"
		# copied from: "/faststorage/project/spider2/AAL_20190722_SOV_correlations/functional_analysis/func_anno_output/concat_all_genes_split_annotations_all_reduced_R"
#funct_annot_all_genes <- "/faststorage/project/spider2/AAL_20190722_SOV_correlations/functional_analysis/gene_names_strong_correlations_pc/gene_fastas/concat_all_genes_split_annotations_all_reduced_R"
all_genes_annot1 <- read.csv(funct_annot_all_genes,header=T,sep="\t")

# import lists of genes from individual tests
files_correlations <- list.files(path="assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/",pattern="*_Exp_fits_slope_patterns.*\\.tab", full.names=T)
files_correlations <- files_correlations[grep("PERMUTED", files_correlations, invert=T)]
files_correlations <- files_correlations[grep(fdr_thresh, files_correlations, invert=T)]
	# in pattern field: dot means any character, while asterix means any number of the preceding character
#"analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/*_Exp_fits_slope_patterns.tab"

###################################################
### code chunk number 3: get list of gene names for universe and individual tests.
###################################################


gene_names <- sapply(strsplit(as.character(all_genes_annot1[,1]),"::"),"[[",1)
all_genes_annot <- data.frame(gene_names,all_genes_annot1)

NROW(unique(sort(gene_names))) == NROW(gene_names) # should be true, otherwise there are duplicate gene names?
	# if false, the duplicates will automatically be removed at a later point.

gene_names_correaltions_list <- list()
for( file in files_correlations	){
	message(file)
	datas <- read.table(file,header=F)
	gene_names_correaltions_list <- append(gene_names_correaltions_list,datas)
	new_name_cor <- paste( sapply(strsplit(basename(file),"\\.|_"), "[[",1), sapply(strsplit(basename(file),"\\.|_"), "[[",5),	sep="_")
	names(gene_names_correaltions_list)[NROW(gene_names_correaltions_list)] <- new_name_cor 
}
	
	
###################################################
### code chunk number 4: construct annotation df and geneset
###################################################

# get first go-term (I assume they are ranked)
first_goTerm <- sapply(strsplit(as.character(all_genes_annot$GOs),","),"[[",1)
all_genes_annot1 <- data.frame(first_goTerm, "ISS", gene_names)	# Why ISS?
goFrame=GOFrame(all_genes_annot1,organism="Stegodyphus_dumicola")
goAllFrame=GOAllFrame(goFrame)
	# setup a gsea gene set collection
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
	# the annotation is within this feature


###################################################
### code chunk number 5: Set up list of functional analysis parameter objects
###################################################

## CONTINUE FROM HERE.
# make universe into only temperature genes. 
genes_temperature_effect <- read.table("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/functional_enrichement/gene_names_w_temperature_effect.tab",header=F)
genes_temperature_effect <- sort(genes_temperature_effect[,1])
# set parameters
universe_gene_ids <- as.character(genes_temperature_effect)	#gene_names
test_gene_ids <- sample(genes_temperature_effect,30)
ontology_type="MF" # "BP" or "MF"
test_direction_is="over"
### make main parameter object
params <- GSEAGOHyperGParams(	name="test_enrichement", geneSetCollection=gsc, 
				geneIds=test_gene_ids,
				universeGeneIds=universe_gene_ids, 
				ontology=ontology_type, pvalueCutoff=0.05, conditional=F,
				testDirection=test_direction_is)#,GeneSetsMaxSize=11000 )
		### OBS duplicate genes in environment is remvoved at this point.
		
# Define ontology abbreviations: 
	# BP = Biological process; 
	# CC = Cellular component; 
	# MF = Molecular function

###################################################
### code chunk number 6: test if gene_subsets are all in gene universe
###################################################


#for( gene_list_nr in 1:NROW(gene_names_correaltions_list) ){
#	list_name <- names(gene_names_correaltions_list[gene_list_nr])
#	print(list_name)
#	# print genes not found in universe
#		#print(sort(as.character(gene_names_correaltions_list[[gene_list_nr]])) %in% sort(universe_gene_ids))
#	print(sort(gene_names_correaltions_list[[gene_list_nr]])[which((sort(as.character(gene_names_correaltions_list[[gene_list_nr]])) %in% sort(universe_gene_ids))==F)])
#}


###################################################
### code chunk number 7: while setting up parameter objects, remove instances of genes not in universe
### 						Also run hyperGTest wihtin loop
###################################################

	# run an lapply hyperGtest --> not working since some of the list elements fails. (probably does not fail for these two datasets, but I just run the same setup as before.)
		# just ran the hyperGTest within the loop. results are called below in lapply.

list_of_test_parameter_objects <- list()
list_of_test_parameter_objects_failing <- list()
list_of_test_results <- list()	
for( gene_list_nr in 1:NROW(gene_names_correaltions_list) ){
	list_name <- names(gene_names_correaltions_list[gene_list_nr])
	print(list_name)
		### removing genes not un universe gene set: (or more correctly, keeping genes found within universe)
	test_gene_ids1 <- as.character(gene_names_correaltions_list[[gene_list_nr]])
	test_gene_ids <- test_gene_ids1[ which(test_gene_ids1 %in% universeGeneIds(params)) ]
	#print(test_gene_ids)
	print(NROW(universeGeneIds(params)))
	if( NROW(test_gene_ids1)-NROW(test_gene_ids) > 0 ){
		message( "Removing ", NROW(test_gene_ids1)-NROW(test_gene_ids),"/",NROW(test_gene_ids1)," (",round((NROW(test_gene_ids1)-NROW(test_gene_ids))/NROW(test_gene_ids1)*100,1),"%)"," genes from test-set, since they are not available in the functional gene universe.")
	}
	geneIds(params) <- test_gene_ids
	print(NROW(geneIds(params)))
		#change original parameter object to have test set of genes
	params_over <- params
	testDirection(params_over) <- "over"
	
	list_of_test_parameter_objects <- append(list_of_test_parameter_objects, params_over)
	names(list_of_test_parameter_objects)[NROW(list_of_test_parameter_objects)] <- paste(list_name,"Over",sep="_")
	# make a parameter object that tests overrepresentation and add to list, lastly change name to reflect correlation

		### Running enrichement test and storing result in list
	list_of_test_results <- append(list_of_test_results, hyperGTest(params_over))
	names(list_of_test_results)[NROW(list_of_test_results)] <- paste(list_name,"Over",sep="_")
}
 #If you want to test for underrepresented GOs, change the term in original params setting and rerun.

results_sum_list <- lapply(list_of_test_results, summary)
	# Understand the count and size parameter to understand if all genes are actually tested? despite the universe being much smaller than it really is?
		#are genes tested in smaller random sets?
		# https://support.bioconductor.org/p/53218/
# maybe remove genes where size and count == 1.


'[1] "CCR_patterns"
[1] "CTmax_patterns"
'

		### Writing results for each comparison to tab separated files 
mapply(function(dname, data) 
   write.table(data, sep="\t", file = paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/functional_enrichement/",
   												ontology_type,"_ontology_",test_direction_is,"REP_", dname, "_", fdr_thresh, ".tab"), row.names = FALSE), 
   paste0("expression_",names(results_sum_list)), results_sum_list)
# it is okay it writes NULL

# Write results to common dataframe for easier interpretation
rm(resultsframe);rm(resultsframe1)
for( result_nr in 1:NROW(results_sum_list) ){
	res_name <- names(results_sum_list[result_nr])
	if( nrow( results_sum_list[[result_nr]] ) > 0 ){
		resultsframe1 <- results_sum_list[[result_nr]]
		resultsframe1$Ont_type_short <- ontology_type
		if( "BP" %in% resultsframe1$Ont_type_short ){
			resultsframe1$Ont_type_long <- "Biological Process"
		}else if( "MF" %in% resultsframe1$Ont_type_short ){
			resultsframe1$Ont_type_long <- "Molecular Function"
		}  
		
		resultsframe1$Data_type <- sapply( strsplit(res_name,"_"), "[[", 1)
		#resultsframe1$Env_axis <- sapply( strsplit(res_name,"_"), "[[", 2)
		resultsframe1$Rep_type <- sapply( strsplit(res_name,"_"), "[[", 3)
	}
	if( nrow( results_sum_list[[result_nr]] ) == 0 ){next}
	if(result_nr==1){
		resultsframe <- resultsframe1
	}else if( !exists("resultsframe") ) {
		resultsframe <- resultsframe1
	}else{
		resultsframe <- rbind(resultsframe,resultsframe1)
	}
}

#if( "refpop" %in% resultsframe$Data_type ) resultsframe[which(resultsframe$Data_type == "refpop"),"Data_type"] <- "WML_PopCorr"
#if( "smps" %in% resultsframe$Data_type ) resultsframe[which(resultsframe$Data_type == "smps"),"Data_type"] <- "WML_CisCorr"
#if( "snp" %in% resultsframe$Data_type ) resultsframe[which(resultsframe$Data_type == "snp"),"Data_type"] <- "SNPs_PopCorr"
#if( "BP" %in% resultsframe$Ontology_type ) resultsframe[which(resultsframe$Ontology_type == "snp"),"Ontology_type"] <- "Biological Process"
#if( "MF" %in% resultsframe$Ontology_type ) resultsframe[which(resultsframe$Ontology_type == "snp"),"Ontology_type"] <- "Molecular function"

resultsframe <- resultsframe[,c(8,9,10,11,1:7)]
nrow(resultsframe)
resultsframe <- resultsframe[-which((resultsframe$Count == 1) & (resultsframe$Size ==1)),]
nrow(resultsframe)
# here - remove if count and size == 1.

write.table(resultsframe, file = paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/Functional_annotation_enrichement_expression/functional_enrichement/",ontology_type,"_ontology_",test_direction_is,"REP_expression_tempTolPatterns","_" ,fdr_thresh, ".tab"), row.names = FALSE,sep="\t")
