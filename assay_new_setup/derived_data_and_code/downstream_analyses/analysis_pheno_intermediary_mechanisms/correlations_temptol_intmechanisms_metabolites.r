#!/usr/bin/env Rscript

##############################################################################
############ 		Correlation heat tol vs. plasticity in metabolites
##############################################################################

	# this script uses functions from the functions_script to investigate
		# which features show intensity level in with similar (mean) paterns
		# as the phenotype Heat tolerance (ctmax) and Cold tolerance (CCR)

# files with intesities of significant metabolites. the main effect is found in the name (temp/pop), as is the treatment type (AQ/ORG & CTmax/CCR)
intensity_files <- dir("assay_new_setup/derived_data_and_code/metabolites/data",pattern="Signi_features_.{3}_|Signi_features_.{4}_",full.names=T,recursive=T)
'[1] "metabolites/data/LCMS_tobias/LCMS_CCR_results/Signi_features_pop_LCMS_CCR.csv"           
 [2] "metabolites/data/LCMS_tobias/LCMS_CCR_results/Signi_features_temp_LCMS_CCR.csv"          
 [3] "metabolites/data/LCMS_tobias/LCMS_CTmax_results/Signi_features_pop_LCMS_CTmax.csv"       
 [4] "metabolites/data/LCMS_tobias/LCMS_CTmax_results/Signi_features_temp_LCMS_CTmax.csv"      
 [5] "metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/Signi_features_pop_NMR_AQ_CCR.csv"       
 [6] "metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/Signi_features_temp_NMR_AQ_CCR.csv"      
 [7] "metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/Signi_features_pop_NMR_AQ_CTmax.csv"   
 [8] "metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/Signi_features_temp_NMR_AQ_CTmax.csv"  
 [9] "metabolites/data/NMR_Kirsten/NMR_ORG_CCR_results/Signi_features_temp_NMR_ORG_CCR.csv"    
[10] "metabolites/data/NMR_Kirsten/NMR_ORG_CTmax_results/Signi_features_temp_NMR_ORG_CTmax.csv"
'
metaCOLOR_files <- dir("assay_new_setup/derived_data_and_code/metabolites/data",pattern="Signi_features_effectANDcolor_.*",full.names=T,recursive=T)
'
[1] "metabolites/data/LCMS_tobias/LCMS_CCR_results/Signi_features_effectANDcolor_pop_LCMS_CCR.csv"           
 [2] "metabolites/data/LCMS_tobias/LCMS_CCR_results/Signi_features_effectANDcolor_temp_LCMS_CCR.csv"          
 [3] "metabolites/data/LCMS_tobias/LCMS_CTmax_results/Signi_features_effectANDcolor_pop_LCMS_CTmax.csv"       
 [4] "metabolites/data/LCMS_tobias/LCMS_CTmax_results/Signi_features_effectANDcolor_temp_LCMS_CTmax.csv"      
 [5] "metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/Signi_features_effectANDcolor_pop_NMR_AQ_CCR.csv"       
 [6] "metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/Signi_features_effectANDcolor_temp_NMR_AQ_CCR.csv"      
 [7] "metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/Signi_features_effectANDcolor_pop_NMR_AQ_CTmax.csv"   
 [8] "metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/Signi_features_effectANDcolor_temp_NMR_AQ_CTmax.csv"  
 [9] "metabolites/data/NMR_Kirsten/NMR_ORG_CCR_results/Signi_features_effectANDcolor_temp_NMR_ORG_CCR.csv"    
[10] "metabolites/data/NMR_Kirsten/NMR_ORG_CTmax_results/Signi_features_effectANDcolor_temp_NMR_ORG_CTmax.csv"'

################################################
### Look at patterns similar to  phenotypes  ###
################################################
# For plastic genes:
	# in relation to CTmax:
		# look for these patterns:
			# 1:	slope of Betta > Otavi > Stampriet > Karasburg >= 0
			# 2:	slope of Otavi > Betta > Stampriet > Karasburg >= 0

	# in relation to CCR:
		# look for these patterns:
			# 1:	slope of  Stampriet & Karasburg > Otavi & Betta > -.05
			
# get only temp files (files containing plastic metabolites) :
intensity_files1 <- intensity_files[grep("temp",intensity_files)]
metaCOLOR_files1 <- metaCOLOR_files[grep("temp",metaCOLOR_files)]
#metabo_type="NMR"

list_of_metabodata <- list()
list_of_metabodata_PERMUTED <- list()
list_of_cor_vectors <- list()
for ( metabolite_file in intensity_files1){
	treatment_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")))[2]
	metabolite_type <-  unlist(         strsplit(basename(metabolite_file), "\\.|_"))[4]
		# get the right files
	if( metabolite_type == "NMR" ){
		solute_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")))[3]
		corresp_metafile1 <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( solute_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
		corresp_metafile <- corresp_metafile1[ grep( treatment_type, corresp_metafile1)	]
	}else if( metabolite_type == "LCMS" ){
		corresp_metafile <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( treatment_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
	}
		#read in the files
	metabo_dat <- read.csv(metabolite_file, header=T,row.names=1)
	meta_data <- read.csv(corresp_metafile,header=T,row.names=1)
		#change column names
	colnames(meta_data) <- c("temperature","population","tempetature_color","pop_color")
		# if LCMS change some of the meta data structure
	if(metabolite_type == "LCMS"){	meta_data$population <- substr(meta_data$population,1,1)	}
	meta_data$pop_temp <- paste0(meta_data$population,meta_data$temperature)

	# Now make the setup for correlations as below:
	intensity_poptemp_mean <- aggregate(metabo_dat, by=list(meta_data$pop_temp), mean)
	rownames(intensity_poptemp_mean) <- intensity_poptemp_mean[,1]		; intensity_poptemp_mean <- intensity_poptemp_mean[,-1]

	# correlate all cols, which genes/metabolites do temperature tolerance correlate to.
	if( metabolite_type == "NMR" ){
		if( treatment_type == "CCR" ){	cor_vec_metabo <- sapply(1:ncol(intensity_poptemp_mean), function(i) cor(intensity_poptemp_mean[,i], CCR_final_mean[,2]))	}
		if( treatment_type == "CTmax" ){	cor_vec_metabo <- sapply(1:ncol(intensity_poptemp_mean), function(i) cor(intensity_poptemp_mean[,i], ToD_final_mean[,2]))	}
		name_name <- paste0(metabolite_type,"_",treatment_type,"_",solute_type)
	}else if( metabolite_type == "LCMS" ){
		if( treatment_type == "CCR" ){	
			CCR_final_mean1 <- CCR_final_mean[ grep("19|29", CCR_final_mean$Group.1) , ]
			cor_vec_metabo <- sapply(1:ncol(intensity_poptemp_mean), function(i) cor(intensity_poptemp_mean[,i], CCR_final_mean1[,2]))
			}
		if( treatment_type == "CTmax" ){	
			ToD_final_mean1 <- ToD_final_mean[ grep("19|29", ToD_final_mean$Group.1) , ]
			cor_vec_metabo <- sapply(1:ncol(intensity_poptemp_mean), function(i) cor(intensity_poptemp_mean[,i], ToD_final_mean1[,2]))
			}
		name_name <- paste0(metabolite_type,"_",treatment_type)
	}
	names(cor_vec_metabo) <- colnames(intensity_poptemp_mean)	
	list_of_cor_vectors[[which(intensity_files1==metabolite_file)]] <- cor_vec_metabo
	names(list_of_cor_vectors)[which(intensity_files1==metabolite_file)] <- name_name
	
	metabo_dat <- as.data.frame(t(metabo_dat))
		# this as data.frame conversion fixed an annoying issue. Now the functions are similar across metabolites and expression. as it should be.
	print(paste(name_name,"Now the script employs functions to calculate population wise slopes, intercept and means for each gene. This may be slow."))
	
	trendlines_data_met <- lapply(1:nrow(metabo_dat), function(i) trendline_pop_specific_dat_func(i, dataset=metabo_dat, meta_dataset=meta_data, meta_treatment_column_name = "temperature", meta_pop_column_name = "population"))
	names(trendlines_data_met)<- paste0(name_name,"_", rownames(metabo_dat))
	
		# Run permutation of data for false positive count
	metabo_dat_PERMUTED <- lapply(1:nrow(metabo_dat), function(i) sample(metabo_dat[i,]))
	metabo_dat_PERMUTED <- as.data.frame(do.call(rbind, metabo_dat_PERMUTED))
	trendlines_data_met_PERMUTED <- lapply(1:nrow(metabo_dat_PERMUTED), function(i) trendline_pop_specific_dat_func(i, dataset=metabo_dat_PERMUTED, meta_dataset=meta_data,meta_treatment_column_name="temperature",meta_pop_column_name="population"))
	names(trendlines_data_met_PERMUTED)<- paste0(name_name,"_", rownames(metabo_dat_PERMUTED))

	if( length(list_of_metabodata) == 0 ){
		list_of_metabodata <- trendlines_data_met
		list_of_metabodata_PERMUTED <- trendlines_data_met_PERMUTED
	}else{
		list_of_metabodata <- append(list_of_metabodata, trendlines_data_met)
		list_of_metabodata_PERMUTED <- append(list_of_metabodata_PERMUTED, trendlines_data_met_PERMUTED)
	}
}




################################################################################################
### Use functions to find genes with similar expression patterns to CTmax
################################################################################################

print("Now the script employs functions to find patterns in slopes similar to the pattern for CTmax.")
# look for this pattern
print(slope_vec_CTmax)
print(ctmax_trends_0diff)

# subset to only CTmax data
list_of_metabodata_CTmax <- list_of_metabodata[grep("CTmax", names(list_of_metabodata))]
# Use NEW functions to find genes with similar expression patterns to CTmax
ctmax_fit_pattern_improved <- lapply(list_of_metabodata_CTmax, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=NULL,	threshold_for_slope=0.02,	pop_vector_zero_trend=ctmax_zero_pops, use_THRESH_for_slopeNotZero=F ))
	# finds genes where K has slope not different from zero, and below (absolute) threshold
	# while S,O,B has slopes different from zero, but in same direction and above/below threshold
NROW(which(sapply(ctmax_fit_pattern_improved,"[")==T))
## OBS: 
	# 6 without threshold on slopes different from zero (5 with threshold)
	# should it not go the other way, when I have a smaller the threshold?
	# No, because, one of the criteria is also that ALL slopes of the other pops (than K) should all be either below the -threshold or above the threshold.
		# loosening the threshold makes more genes fullfill this criteria, than it sort away based on the K slope. 
CTMAX_onlyFIT <- list_of_metabodata_CTmax[which(sapply(ctmax_fit_pattern_improved,"[")==T)]
print(paste0(NROW(CTMAX_onlyFIT)," Metabolites fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved)," metabolites."))
write.table(names(CTMAX_onlyFIT),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/CTmax_metab_fits_slope_patterns_.tab"), row.names=F,col.names=F)
'"6 Metabolites fit the improved CTmax pattern (Nov22), out of 239 metabolites."'

	# Check number of false positive on this permuted dataset (permuted intensity values per metabolite/peak)
list_of_metabodata_CTmax_PERMUTED <- list_of_metabodata_PERMUTED[grep("CTmax", names(list_of_metabodata_PERMUTED))]
ctmax_fit_pattern_improved_PERMUTED <- lapply(list_of_metabodata_CTmax_PERMUTED, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=NULL,	threshold_for_slope=0.02,	pop_vector_zero_trend=ctmax_zero_pops, use_THRESH_for_slopeNotZero=F ))
NROW(which(sapply(ctmax_fit_pattern_improved_PERMUTED,"[")==T))
CTMAX_onlyFIT_PERMUTED <- list_of_metabodata_CTmax_PERMUTED[which(sapply(ctmax_fit_pattern_improved_PERMUTED,"[")==T)]
print(paste0(NROW(CTMAX_onlyFIT_PERMUTED)," PERMUTED Metabolites fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved_PERMUTED)," metabolites. (Estimate of False Positives)"))
'"0 PERMUTED Metabolites fit the improved CTmax pattern (Nov22), out of 239 metabolites. (Estimate of False Positives)"'

#############################
# now plot all genes to see if it is true no pattern matches
#############################

print("Now the genes fitting the general pattern is plotted - so we can see if it makes sense.")
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,6,4,2)]	# B, E/O,	S, K, 
pop_vector <- c("B","O","S","K")
popcols <- popcols[order(pop_vector)]
pop_vector <- pop_vector[order(pop_vector)]
if( any((pop_vector) == "B") ) {pop_vector[which(pop_vector == "B")] <- "Betta"}
if( any((pop_vector) == "K") ) {pop_vector[which(pop_vector == "K")] <- "Karasburg"}
if( any((pop_vector) == "S") ) {pop_vector[which(pop_vector == "S")] <- "Stampriet"}
if( any((pop_vector) == "O") ) {pop_vector[which(pop_vector == "O")] <- "Otavi"}
		# loop to separate plots in files of 20 each
for ( metabolite_file in intensity_files1){
	treatment_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")	)	)[2]
	metabolite_type <-  unlist(         strsplit(basename(metabolite_file), "\\.|_")    )[4]
	if( treatment_type == "CCR"){
		next
	}
		# get the right files
	if( metabolite_type == "NMR" ){
		solute_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")	)	)[3]
		corresp_metafile1 <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( solute_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
		corresp_metafile <- corresp_metafile1[ grep( treatment_type, corresp_metafile1)	]
		# subset slopes dataset to fit current data:
		list_of_metabodata_SUB <- list_of_metabodata[grep(metabolite_type,names(list_of_metabodata))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(treatment_type,names(list_of_metabodata_SUB))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(solute_type,names(list_of_metabodata_SUB))]
		name_data <- paste0(metabolite_type,"_",treatment_type,"_",solute_type)
	}else if( metabolite_type == "LCMS" ){
		corresp_metafile <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( treatment_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
		# subset slopes dataset to fit current data:
		list_of_metabodata_SUB <- list_of_metabodata[grep(metabolite_type,names(list_of_metabodata))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(treatment_type,names(list_of_metabodata_SUB))]
		name_data <- paste0(metabolite_type,"_",treatment_type)
	}
		#read in the files
	metabo_dat <- read.csv(metabolite_file, header=T,row.names=1)
	meta_data <- read.csv(corresp_metafile,header=T,row.names=1)
		#change column names
	colnames(meta_data) <- c("temperature","population","tempetature_color","pop_color")
	metabo_dat <- as.data.frame(t(metabo_dat))
	rownames(metabo_dat) <- paste0(name_data,"_", rownames(metabo_dat))
		# if LCMS change some of the meta data structure
	if(metabolite_type == "LCMS"){	meta_data$population <- substr(meta_data$population,1,1)	}
	meta_data$pop_temp <- paste0(meta_data$population,meta_data$temperature)

	#Plot the intensities of ALL metabolites
	counter=0
	name_vector <- sprintf("%s_%02d", paste0(name_data,"_metabolites_patternsfit_nov22_slope_ALL"), 1:100)
	for (rows in seq(NROW(list_of_metabodata_SUB))	){
		if ( (rows %in%  seq(0, NROW(list_of_metabodata_SUB), 20) ) || (rows == NROW(list_of_metabodata_SUB) )  	){
			# if rows are 20, 40 60 etc, or if it is the last row
			if( rows != 1){ # if different from 1
				counter=counter+1
				png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[counter],".png"),
							width=8.5,height=5.5,units="in",res=400)
					#Open picture file
				#message(rows)
				if (rows == NROW(list_of_metabodata_SUB)){
						#If last row
					seq_start <- NROW(list_of_metabodata_SUB) - (	rows%%seq(0, NROW(list_of_metabodata_SUB), 20)[	NROW(seq(0, NROW(list_of_metabodata_SUB), 20))	]	) + 1
				}else{
					seq_start <- seq(0, NROW(list_of_metabodata_SUB), 20) [ 	which( seq(0, NROW(list_of_metabodata_SUB), 20) == rows) 	] - 19			}
						#get starting number (rownumber to plot)
				#message(seq_start)			
				number_vector <- seq_start:rows
				if( NROW(number_vector) <= 5){
					par(mfrow=c( ceiling((NROW(number_vector)/5)) , NROW(number_vector)),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				}else{
					par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				}
				#par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(3,3,2.5,1), oma=c(2,2,2,0))
				for(each in number_vector){
					plot_count_is=which(number_vector == each)
					plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=list_of_metabodata_SUB, plot_count=plot_count_is, ylab_text="Metabolite intensity")	
						#gene number/ feature number/row number in the slopes dataset
					}
				#message(name_vector[counter]) # if I need to check numbers.
				par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
    			plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
				legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
				dev.off()
	}	}	}

		# plot the intensities of only significant metabolites.
	
	CTMAX_onlyFIT_SUB <- CTMAX_onlyFIT[grep(name_data, names(CTMAX_onlyFIT))]
	if( length(CTMAX_onlyFIT_SUB) > 0 ){
		counter=0
		name_vector <- sprintf("%s_%02d", paste0(name_data,"_metabolites_patternsfit_nov22_slope_FIT"), 1:100)
		for (rows in seq(NROW(CTMAX_onlyFIT_SUB))	){
			if( length(seq(NROW(CTMAX_onlyFIT_SUB))) < 20	){
				counter=counter+1
				png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[1],".png"),
							width=8.5,height=5.5,units="in",res=400)
				par(mfrow=c( ceiling(length(seq(NROW(CTMAX_onlyFIT_SUB)))/4) , length(seq(NROW(CTMAX_onlyFIT_SUB))) / ceiling(length(seq(NROW(CTMAX_onlyFIT_SUB)))/4)),
						mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				for(each in seq(NROW(CTMAX_onlyFIT_SUB))){
					plot_count_is=which(seq(NROW(CTMAX_onlyFIT_SUB)) == each)
					plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=CTMAX_onlyFIT_SUB , plot_count=plot_count_is, ylab_text="Metabolite intensity")			
							#gene number/ feature number/row number in the slopes dataset
				}
				par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
				plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
				legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
				dev.off()
				next

			}else if ( rows %in%  seq(0, NROW(CTMAX_onlyFIT_SUB), 20) || rows == NROW(CTMAX_onlyFIT_SUB) 	){
				# if rows are 20, 40 60 etc, or if it is the last row
				if( rows != 1){ # if different from 1
					counter=counter+1
					png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[counter],".png"),
								width=8.5,height=5.5,units="in",res=400)
						#Open picture file
					#message(rows)
					if (rows == NROW(CTMAX_onlyFIT_SUB)){
							#If last row
						seq_start <- NROW(CTMAX_onlyFIT_SUB) - (	rows%%seq(0, NROW(CTMAX_onlyFIT_SUB), 20)[	NROW(seq(0, NROW(CTMAX_onlyFIT_SUB), 20))	]	) + 1
					}else{
						seq_start <- seq(0, NROW(CTMAX_onlyFIT_SUB), 20) [ 	which( seq(0, NROW(CTMAX_onlyFIT_SUB), 20) == rows) 	] - 19			}
							#get starting number (rownumber to plot)
					#message(seq_start)			
					number_vector <- seq_start:rows
					if( NROW(number_vector) <= 5){
						par(mfrow=c( ceiling((NROW(number_vector)/5)) , NROW(number_vector)),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
					}else{
						par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
					}
					#par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(3,3,2.5,1), oma=c(2,2,2,0))
					for(each in number_vector){
						plot_count_is=which(number_vector == each)
						plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=CTMAX_onlyFIT_SUB, plot_count=plot_count_is, ylab_text="Metabolite intensity")			
							#gene number/ feature number/row number in the slopes dataset
						}
					#message(name_vector[counter]) # if I need to check numbers.
					par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
					plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
					legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
					dev.off()
		}	}	}
	}
}











################################################################################################
###  CCR! Use functions to find genes with similar expression patterns to CCR
################################################################################################

print("Now the script employs functions to find patterns in slopes similar to the pattern for CTmax.")
# look for this pattern
print(slope_vec_CCR)
print(ccr_trends_0diff)
print(ccr_zero_pops)
# subset to only CCR data
list_of_metabodata_CCR <- list_of_metabodata[grep("CCR", names(list_of_metabodata))]
# Use NEW functions to find genes with similar expression patterns to CTmax
CCR_fit_pattern_improved <- lapply(list_of_metabodata_CCR, function(i) fit_CCR_improved_filter(genewise_dataframe=i,	threshold_for_mean=NULL,	threshold_for_slope=.1,	pop_vector_zero_trend=ccr_zero_pops, use_THRESH_for_slopeNotZero=F ))
#OBS:
	# finds genes where K has slope not different from zero, and below (absolute) threshold
	# while S,O,B has slopes different from zero, but in same direction and above/below threshold
NROW(which(sapply(CCR_fit_pattern_improved,"[")==T))
## OBS: 1 without threshold on slope for slopes different from 0.
CCR_onlyFIT <- list_of_metabodata_CCR[which(sapply(CCR_fit_pattern_improved,"[")==T)]
print(paste0(NROW(CCR_onlyFIT)," Metabolites fit the improved CCRTemp pattern (Nov22), out of ",NROW(CCR_fit_pattern_improved)," metabolites."))
write.table(names(CCR_onlyFIT),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/CCR_metab_fits_slope_patterns_.tab"), row.names=F,col.names=F)
'"1 Metabolites fit the improved CCRTemp pattern (Nov22), out of 197 metabolites."'

# Estimate false positives:
list_of_metabodata_CCR_PERMUTED <- list_of_metabodata_PERMUTED[grep("CCR", names(list_of_metabodata_PERMUTED))]
CCR_fit_pattern_improved_PERMUTED <- lapply(list_of_metabodata_CCR_PERMUTED, function(i) fit_CCR_improved_filter(genewise_dataframe=i,	threshold_for_mean=NULL,	threshold_for_slope=.1,	pop_vector_zero_trend=ccr_zero_pops, use_THRESH_for_slopeNotZero=F ))
NROW(which(sapply(CCR_fit_pattern_improved_PERMUTED,"[")==T))
CCR_onlyFIT_PERMUTED <- list_of_metabodata_CCR_PERMUTED[which(sapply(CCR_fit_pattern_improved_PERMUTED,"[")==T)]
print(paste0(NROW(CCR_onlyFIT_PERMUTED)," PERMUTED Metabolites fit the improved CCRTemp pattern (Nov22), out of ",NROW(CCR_fit_pattern_improved_PERMUTED)," metabolites. (False positive estimate)"))
'"3 PERMUTED Metabolites fit the improved CCRTemp pattern (Nov22), out of 197 metabolites. (False positive estimate)"'
#############################
# now plot all genes to see if it is true no pattern matches
#############################

print("Now the genes fitting the general pattern is plotted - so we can see if it makes sense.")
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,6,4,2)]	# B, E/O,	S, K, 
pop_vector <- c("B","O","S","K")
popcols <- popcols[order(pop_vector)]
pop_vector <- pop_vector[order(pop_vector)]
if( any((pop_vector) == "B") ) {pop_vector[which(pop_vector == "B")] <- "Betta"}
if( any((pop_vector) == "K") ) {pop_vector[which(pop_vector == "K")] <- "Karasburg"}
if( any((pop_vector) == "S") ) {pop_vector[which(pop_vector == "S")] <- "Stampriet"}
if( any((pop_vector) == "O") ) {pop_vector[which(pop_vector == "O")] <- "Otavi"}
		# loop to separate plots in files of 20 each
for ( metabolite_file in intensity_files1){
	treatment_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")	)	)[2]
	metabolite_type <-  unlist(         strsplit(basename(metabolite_file), "\\.|_")    )[4]
	if( treatment_type == "CTmax"){
		next
	}
		# get the right files
	if( metabolite_type == "NMR" ){
		solute_type <- rev(	unlist(		strsplit(basename(metabolite_file), "\\.|_")	)	)[3]
		corresp_metafile1 <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( solute_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
		corresp_metafile <- corresp_metafile1[ grep( treatment_type, corresp_metafile1)	]
		# subset slopes dataset to fit current data:
		list_of_metabodata_SUB <- list_of_metabodata[grep(metabolite_type,names(list_of_metabodata))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(treatment_type,names(list_of_metabodata_SUB))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(solute_type,names(list_of_metabodata_SUB))]
		name_data <- paste0(metabolite_type,"_",treatment_type,"_",solute_type)
	}else if( metabolite_type == "LCMS" ){
		corresp_metafile <- metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] [grep( treatment_type, metaCOLOR_files1[grep(metabolite_type,metaCOLOR_files1)] )]
		# subset slopes dataset to fit current data:
		list_of_metabodata_SUB <- list_of_metabodata[grep(metabolite_type,names(list_of_metabodata))]
		list_of_metabodata_SUB <- list_of_metabodata_SUB[grep(treatment_type,names(list_of_metabodata_SUB))]
		name_data <- paste0(metabolite_type,"_",treatment_type)
	}
		#read in the files
	metabo_dat <- read.csv(metabolite_file, header=T,row.names=1)
	meta_data <- read.csv(corresp_metafile,header=T,row.names=1)
		#change column names
	colnames(meta_data) <- c("temperature","population","tempetature_color","pop_color")
	metabo_dat <- as.data.frame(t(metabo_dat))
	rownames(metabo_dat) <- paste0(name_data,"_", rownames(metabo_dat))
		# if LCMS change some of the meta data structure
	if(metabolite_type == "LCMS"){	meta_data$population <- substr(meta_data$population,1,1)	}
	meta_data$pop_temp <- paste0(meta_data$population,meta_data$temperature)

	#Plot the intensities of ALL metabolites 
	counter=0
	name_vector <- sprintf("%s_%02d", paste0(name_data,"_metabolites_patternsfit_nov22_slope_ALL"), 1:100)
	for (rows in seq(NROW(list_of_metabodata_SUB))	){
		if( rows == NROW(list_of_metabodata_SUB) && rows %in% c(seq(0, NROW(list_of_metabodata_SUB), 20)+1) ){
				counter=counter+1
				metabodat_sub_list <- list_of_metabodata_SUB[rows]
				png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[counter],".png"),
							width=8.5,height=5.5,units="in",res=400)
				par(mfrow=c( ceiling(length(seq(NROW(metabodat_sub_list)))/4) , length(seq(NROW(metabodat_sub_list))) / ceiling(length(seq(NROW(metabodat_sub_list)))/4)),
						mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				for(each in seq(NROW(metabodat_sub_list))){
					plot_count_is=which(seq(NROW(metabodat_sub_list)) == each)
					plot_genes_with_pattern_improved(gene_nr=each, intensity_data=metabo_dat, meta_data_frame=meta_data, meta_dat_column_name="temperature",cols_vector=popcols, dataset_slopes_intercept=metabodat_sub_list , plot_count=plot_count_is, ylab_text="Metabolite intensity")			
							#gene number/ feature number/row number in the slopes dataset
							# OBS SOME ERROR HERE: not names(list_of_metabodata_SUB) %in% rownames(metabo_dat) when metabolite_data=intensity_files[1]
				}
				par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
				plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
				legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
				dev.off()
				next
		}else if( rows %in%  seq(0, NROW(list_of_metabodata_SUB), 20) || rows == NROW(list_of_metabodata_SUB)  	){
			# if rows are 20, 40 60 etc, or if it is the last row
			if( rows != 1){ # if different from 1
				counter=counter+1
				png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[counter],".png"),
							width=8.5,height=5.5,units="in",res=400)
					#Open picture file
				#message(rows)
				if (rows == NROW(list_of_metabodata_SUB)){
						#If last row
					seq_start <- NROW(list_of_metabodata_SUB) - (	rows%%seq(0, NROW(list_of_metabodata_SUB), 20)[	NROW(seq(0, NROW(list_of_metabodata_SUB), 20))	]	) + 1
				}else{
					seq_start <- seq(0, NROW(list_of_metabodata_SUB), 20) [ 	which( seq(0, NROW(list_of_metabodata_SUB), 20) == rows) 	] - 19			}
						#get starting number (rownumber to plot)
				#message(seq_start)			
				number_vector <- seq_start:rows
				#c( ceiling(length(seq(NROW(metabodat_sub_list)))/4) , length(seq(NROW(metabodat_sub_list))) / ceiling(length(seq(NROW(metabodat_sub_list)))/4))
				if( NROW(number_vector) <= 5){
					par(mfrow=c( ceiling((NROW(number_vector)/5)) , NROW(number_vector)),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				}else{
					par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
				}
				for(each in number_vector){
					plot_count_is=which(number_vector == each)
					plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=list_of_metabodata_SUB, plot_count=plot_count_is, ylab_text="Metabolite intensity")			
						#gene number/ feature number/row number in the slopes dataset
				}
				#message(name_vector[counter]) # if I need to check numbers.
				par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
				plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
				legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
				dev.off()
	}	}	}

		# plot the intensities of only significant metabolites.
	
	CCR_onlyFIT_SUB <- CCR_onlyFIT[grep(name_data, names(CCR_onlyFIT))]
	if( length(CCR_onlyFIT_SUB) > 0 ){
		counter=0
		name_vector <- sprintf("%s_%02d", paste0(name_data,"_metabolites_patternsfit_nov22_slope_FIT"), 1:100)
		for (rows in seq(NROW(CCR_onlyFIT_SUB))	){
			if( length(seq(NROW(CCR_onlyFIT_SUB))) < 20	){
				counter=counter+1
				png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[1],".png"),
							width=8.5,height=5.5,units="in",res=400)
				par(mfrow=c( ceiling(length(seq(NROW(CCR_onlyFIT_SUB)))/4) , ceiling( length(seq(NROW(CCR_onlyFIT_SUB))) / ceiling(length(seq(NROW(CCR_onlyFIT_SUB)))/4) ) ),
						mar=c(2.1,2,2,1), oma=c(2,2.5,2,0) )
				for(each in seq(NROW(CCR_onlyFIT_SUB))){
					plot_count_is=which(seq(NROW(CCR_onlyFIT_SUB)) == each)
					plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=CCR_onlyFIT_SUB, plot_count=plot_count_is, ylab_text="Metabolite intensity")			
							#gene number/ feature number/row number in the slopes dataset
				}
				par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
				plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
				legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
				dev.off()
				next

			}else if ( rows %in%  seq(0, NROW(CCR_onlyFIT_SUB), 20) || rows == NROW(CCR_onlyFIT_SUB) 	){
				# if rows are 20, 40 60 etc, or if it is the last row
				if( rows != 1){ # if different from 1
					counter=counter+1
					png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/metabolites_vs_temperature_tolerances/",name_vector[counter],".png"),
								width=8.5,height=5.5,units="in",res=400)
						#Open picture file
					#message(rows)
					if (rows == NROW(CCR_onlyFIT_SUB)){
							#If last row
						seq_start <- NROW(CCR_onlyFIT_SUB) - (	rows%%seq(0, NROW(CCR_onlyFIT_SUB), 20)[	NROW(seq(0, NROW(CCR_onlyFIT_SUB), 20))	]	) + 1
					}else{
						seq_start <- seq(0, NROW(CCR_onlyFIT_SUB), 20) [ 	which( seq(0, NROW(CCR_onlyFIT_SUB), 20) == rows) 	] - 19			}
							#get starting number (rownumber to plot)
					#message(seq_start)			
					number_vector <- seq_start:rows
					if( NROW(number_vector) <= 5){
						par(mfrow=c( ceiling((NROW(number_vector)/5)) , NROW(number_vector)),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
					}else{
						par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
					}
					#par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(3,3,2.5,1), oma=c(2,2,2,0))
					for(each in number_vector){
						plot_count_is=which(number_vector == each)
						plot_genes_with_pattern_improved(gene_nr=each,intensity_data=metabo_dat,meta_data_frame=meta_data,meta_dat_column_name="temperature",cols_vector=popcols,dataset_slopes_intercept=CCR_onlyFIT_SUB, plot_count=plot_count_is, ylab_text="Metabolite intensity")			
							#gene number/ feature number/row number in the slopes dataset
						}
					#message(name_vector[counter]) # if I need to check numbers.
					par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
					plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
					legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
					dev.off()
		}	}	}
	}
}


print(paste0(NROW(CTMAX_onlyFIT)," Metabolites fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved)," metabolites."))
print(paste0(NROW(CCR_onlyFIT)," Metabolites fit the improved CCR pattern (Nov22), out of ",NROW(CCR_fit_pattern_improved)," metabolites."))
print("Done")

'
[1] "6 Metabolites fit the improved CTmax pattern (Nov22), out of 239 metabolites."
[1] "1 Metabolites fit the improved CCR pattern (Nov22), out of 197 metabolites."
'


'
# histogram of correlations ## OBS ADD type to name and only plot relevant treatment (CCR/CTmax)
	pdf(file="hist_corcoef_tolerances_plastic_metabolites.pdf", width=7, height=5)
		# potentially plot four graphs: only plastic responses : (data saved to list)
			# AQ CTmax and CCR
			# ORG CTmax and CCR
	par(mfrow=c(2,2),cex=.7)
	lapply(names(list_of_cor_vectors), function(x) hist(list_of_cor_vectors[[x]], main = x, xlab="Correlation coefficient"))
	dev.off()

	png(file="hist_corcoef_tolerances_plastic_metabolites.png", width=7, height=5,units="in",res=400)
	par(mfrow=c(2,2),cex=.7)
	lapply(names(list_of_cor_vectors), function(x) hist(list_of_cor_vectors[[x]], main = x, xlab="Correlation coefficient"))
	dev.off()
	'