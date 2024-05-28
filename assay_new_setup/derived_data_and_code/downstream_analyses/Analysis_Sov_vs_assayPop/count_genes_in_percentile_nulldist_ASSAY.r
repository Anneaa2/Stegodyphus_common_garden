#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript

## For PCA results colnames

# Count number of genes in real distribution above XX percentile of the permuted distribution

list_name_overlap = list()
list_name_overlap_coef_pval = list()
if (analysis_type == "mantel"){
	for( datadist in names(listofdists)[grep("RealCor_COR", names(listofdists))] ){
		message("\n",datadist)
		secdatadist_name <- sapply(strsplit(datadist,"RealCor_COR_"),"[[",2)
		names_togrepfrom <- names(listofdists)[	grep(secdatadist_name, names(listofdists))]
		header_name <- secdatadist_name# sapply(strsplit(secdatadist_name,"Cor_"),"[[",2)
		for( secdatadist in names_togrepfrom[	grep("RandCor_COR",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
			message(secdatadist)
			perctile_sec <- quantile(listofdists[[secdatadist]],probs=percentile_is)
			#for( thirddist in names_togrepfrom[	grep("RandCor_NO",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
			#	message(thirddist)
			#	perctile_third <- quantile(listofdists[[thirddist]],probs=percentile_is)
				if (datadist == names(listofdists)[grep("RealCor_COR", names(listofdists))][1] ){
					rownamesare <- header_name			
						# how many genes/asv pass both significance and coef threhold
					passvector_both <- c(which(listofdists[[datadist]]>perctile_sec),which(listofp[[datadist]]<.05))
					passvector_both_sec <- c(which(listofdists[[secdatadist]]>perctile_sec),which(listofp[[secdatadist]]<.05))
					passvector_SIGNIF <- which(listofp[[datadist]]<.05)
					passvector_SIGNIF_sec <- which(listofp[[secdatadist]]<.05)
					theoretical_exp_N <- NROW(listofdists[[datadist]])*(1-percentile_is)
						# add data to frame
					dataset_above_perctile <- c(	NROW(which(listofdists[[datadist]]>perctile_sec)),
													NROW(which(listofdists[[secdatadist]]>perctile_sec)),				
													perctile_sec,		
													#NROW(which(listofdists[[datadist]]>perctile_third)),
													#NROW(which(listofdists[[thirddist]]>perctile_third)),												
													#perctile_third,
													#NROW(passvector_both[which(duplicated(passvector_both))]),
													#NROW(passvector_both_sec[which(duplicated(passvector_both_sec))]),
													#NROW(passvector_SIGNIF),
													#NROW(passvector_SIGNIF_sec),
													theoretical_exp_N,
													NROW(listofdists[[datadist]])
													)
					list_name_overlap[[datadist]] <- list_overlap[[datadist]][		which(listofdists[[datadist]]>perctile_sec)		]
					list_name_overlap_coef_pval[[datadist]] <- droplevels(list_overlap[[datadist]][passvector_both[which(duplicated(passvector_both))]])
					
				}else{
					rownamesare <- append(rownamesare, header_name)
					passvector_both <- c(which(listofdists[[datadist]]>perctile_sec),which(listofp[[datadist]]<.05))
					passvector_both_sec <- c(which(listofdists[[secdatadist]]>perctile_sec),which(listofp[[secdatadist]]<.05))
					passvector_SIGNIF <- which(listofp[[datadist]]<.05)
					passvector_SIGNIF_sec <- which(listofp[[secdatadist]]<.05)
					theoretical_exp_N <- NROW(listofdists[[datadist]])*(1-percentile_is)
					dataset_above_perctile <- rbind( dataset_above_perctile,	c(	NROW(which(listofdists[[datadist]]>perctile_sec)),			
																					NROW(which(listofdists[[secdatadist]]>perctile_sec)),	
																					perctile_sec,		
																					#NROW(which(listofdists[[datadist]]>perctile_third)),		
																					#NROW(which(listofdists[[thirddist]]>perctile_third)),		
																					#perctile_third,
																					#NROW(passvector_both[which(duplicated(passvector_both))]),
																					#NROW(passvector_both_sec[which(duplicated(passvector_both_sec))]),
																					#NROW(passvector_SIGNIF),
																					#NROW(passvector_SIGNIF_sec),
																					theoretical_exp_N,
																					NROW(listofdists[[datadist]])
																				)	)
					list_name_overlap[[datadist]] <- list_overlap[[datadist]][		which(listofdists[[datadist]]>perctile_sec)		]
					overlapis <- which(duplicated(c(which(listofdists[[datadist]]>perctile_sec),which(listofp[[datadist]]<.05))))
					list_name_overlap_coef_pval[[datadist]] <- droplevels(list_overlap[[datadist]][passvector_both[which(duplicated(passvector_both))]])
				}
				#message( "Real distribution: N( sign ): ",NROW(which(listofdists[[datadist]]<.05)))
				#message( "Permuted distribution: N( sign ): ",NROW(which(listofdists[[secdatadist]]<.05)))
				#message( "Permuted distribution, non corrected: N( sign ): ",NROW(which(listofdists[[thirddist]]<.05)))
			#}	
		}	
	}	
	colnames_dataset <- c(	paste0("N_Real>",percentile_is,"_Perm"),
							paste0("N_Perm>",percentile_is,"_Perm"),
							paste0("Perm_",percentile_is),
							#paste0("N_Real>",percentile_is,"_Perm_NOcor"),
							#paste0("N_Perm_NOcor>",percentile_is,"_Perm_NOcor"),
							#paste0("Perm_NOcor_",percentile_is),
							#paste0("N_Real>",percentile_is,"_Perm_&pass_Signif"),
							#paste0("N_Perm>",percentile_is,"_Perm_&pass_Signif"),
							#paste0("N_Real<Signif"),
							#paste0("N_Perm<Signif"),
							paste0("Theoretical_exp_perctile"),
							paste0("N_correlations_run")
							)
}



						
dataset_above_perctile <- as.data.frame(dataset_above_perctile)
colnames(dataset_above_perctile) <- colnames_dataset ;		rownames(dataset_above_perctile) <- rownamesare
dataset_above_perctile
#barplot(dataset_above_perctile[,1],names.arg=rownamesare)
percentile_is1 <- sapply(strsplit(as.character(percentile_is),"\\."),"[[",2)
write.csv(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/N_above_",percentile_is1,"percentile","_",analysis_type,".csv"), dataset_above_perctile)


############################################################################################
# data frame for plotting the number of genes that pass coef threshold:
############################################################################################
percent_below_axis <- .15 	# legend placement: fraction below x axis ie. fraction of y axis below 0
dataset_above_perctile_sub <- data.frame(	dataset_above_perctile[1:5,1],
											dataset_above_perctile[6:10,1],
											dataset_above_perctile[11:15,1],
											dataset_above_perctile[16:20,1],
											dataset_above_perctile[21:25,1],
											dataset_above_perctile[26:30,1])
rownames(dataset_above_perctile_sub) <- c("PC1","PC2","PC3","PC4","PC5")
colnames(dataset_above_perctile_sub) <- unique(sapply(strsplit(rownames(dataset_above_perctile),"_PC"),"[[",1))
write.csv(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/df_N_above_",percentile_is1,"percentile","_",analysis_type,".csv"), dataset_above_perctile_sub)


	# data frame for plotting the number of genes that pass coef threshold in Permuted distribution (EXPECTED)
dataset_exp_above_perctile_sub <- data.frame(	dataset_above_perctile[1:5,2],
											dataset_above_perctile[6:10,2],
											dataset_above_perctile[11:15,2],
											dataset_above_perctile[16:20,2],
											dataset_above_perctile[21:25,2],
											dataset_above_perctile[26:30,2])
rownames(dataset_exp_above_perctile_sub) <- c("PC1","PC2","PC3","PC4","PC5")
colnames(dataset_exp_above_perctile_sub) <- unique(sapply(strsplit(rownames(dataset_above_perctile),"_PC"),"[[",1))
write.csv(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/df_N_above_inNULLdisttrib_",percentile_is1,"percentile","_",analysis_type,".csv"), dataset_exp_above_perctile_sub)
	


############################################################################################
# data frame for plotting the number of genes that pass coef threshold in Permuted distribution (EXPECTED)
############################################################################################
dataset_exp_above_perctile_theory <- data.frame(	dataset_above_perctile[1:5,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )],
											dataset_above_perctile[6:10,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )],
											dataset_above_perctile[11:15,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )],
											dataset_above_perctile[16:20,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )],
											dataset_above_perctile[21:25,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )],
											dataset_above_perctile[26:30,	which(colnames(dataset_above_perctile) == "Theoretical_exp_perctile" )])
rownames(dataset_exp_above_perctile_theory) <- c("PC1","PC2","PC3","PC4","PC5")
colnames(dataset_exp_above_perctile_theory) <- unique(sapply(strsplit(rownames(dataset_above_perctile),"_PC"),"[[",1))
write.csv(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/df_N_EXP_above_",percentile_is1,"percentile","_",analysis_type,".csv"), dataset_exp_above_perctile_theory)

############################################################################################
# data frame with number of correlations run
############################################################################################
dataset_N_correlations <- data.frame(		dataset_above_perctile[1:5,	which(colnames(dataset_above_perctile) == "N_correlations_run" )],
											dataset_above_perctile[6:10,	which(colnames(dataset_above_perctile) == "N_correlations_run" )],
											dataset_above_perctile[11:15,	which(colnames(dataset_above_perctile) == "N_correlations_run" )],
											dataset_above_perctile[16:20,	which(colnames(dataset_above_perctile) == "N_correlations_run" )],
											dataset_above_perctile[21:25,	which(colnames(dataset_above_perctile) == "N_correlations_run" )],
											dataset_above_perctile[26:30,	which(colnames(dataset_above_perctile) == "N_correlations_run" )])
rownames(dataset_N_correlations) <- c("PC1","PC2","PC3","PC4","PC5")
colnames(dataset_N_correlations) <- unique(sapply(strsplit(rownames(dataset_above_perctile),"_PC"),"[[",1))
write.csv(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/df_N_correlations_",percentile_is1,"percentile","_",analysis_type,".csv"), dataset_N_correlations)



#if ( analysis_type	== "mrm" ) {
#	dataset_exp_above_perctile_theory <- dataset_exp_above_perctile_theory[,	-which( is.na(colnames(dataset_exp_above_perctile_theory)))]
	#dataset_exp_below_SIG_sub	<- dataset_exp_below_SIG_sub[,	-which( is.na(colnames(dataset_exp_below_SIG_sub)))]
	#dataset_below_SIG_sub	<- dataset_below_SIG_sub[,	-which( is.na(colnames(dataset_below_SIG_sub)))]
	#dataset_exp_above_perctile_andSIG_sub	<- dataset_exp_above_perctile_andSIG_sub[,	-which( is.na(colnames(dataset_exp_above_perctile_andSIG_sub)))]
	#dataset_above_perctile_andSIG_sub	<- dataset_above_perctile_andSIG_sub[,	-which( is.na(colnames(dataset_above_perctile_andSIG_sub)))]
#	dataset_exp_above_perctile_sub	<- dataset_exp_above_perctile_sub[,	-which( is.na(colnames(dataset_exp_above_perctile_sub)))]
#	dataset_above_perctile_sub	<- dataset_above_perctile_sub[,	-which( is.na(colnames(dataset_above_perctile_sub)))]
#	dataset_N_correlations	<- dataset_N_correlations[,	-which( is.na(colnames(dataset_N_correlations)))]
#}