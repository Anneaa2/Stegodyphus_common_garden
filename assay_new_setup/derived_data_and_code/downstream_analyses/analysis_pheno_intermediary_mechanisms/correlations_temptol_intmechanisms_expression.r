#!/usr/bin/env Rscript

##############################################################################
############ 		Correlation heat tol vs. plasticity in expression
##############################################################################

	# this script uses functions from the functions_script to investigate
		# which genes show expression level in with similar (mean) paterns
		# as the phenotype Heat tolerance (ctmax) and Cold tolerance (CCR)

#----------------

## Correlation heat tol vs. plasticity in expression
	# correlate all cols, which genes expression level do temperature tolerance correlate to.
	print("First making simple correlations.")
	cor_vec_exp_CCR <- sapply(1:ncol(plastic_counts_mean), function(i) cor(plastic_counts_mean[,i], CCR_final_mean[,2]))
	cor_vec_exp_Ctmax <- sapply(1:ncol(plastic_counts_mean), function(i) cor(plastic_counts_mean[,i], ToD_final_mean[,2]))

	
	#make named vectors
	names(cor_vec_exp_CCR) <- colnames(plastic_counts_mean)
	names(cor_vec_exp_Ctmax) <- colnames(plastic_counts_mean)
		
	# subset to high absolute correlations
	cutoff_correlation <- .75
	cor_vec_exp_Ctmax_sub <- cor_vec_exp_Ctmax[which(abs(cor_vec_exp_Ctmax) > cutoff_correlation)]
	cor_vec_exp_CCR_sub <- cor_vec_exp_CCR[which(abs(cor_vec_exp_CCR) > cutoff_correlation)]

		
# histogram of correlations
		# when plottet looks bimodal
	pdf(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/hist_corcoef_tolerances_expression_",padj_threshold,".pdf"), width=7, height=5)
	par(mfrow=c(1,2),cex=.7)
	hist(cor_vec_exp_Ctmax)
	hist(cor_vec_exp_CCR)
	dev.off()
'
	png(file=paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/hist_corcoef_tolerances_expression_",padj_threshold,".png"), width=7, height=5,units="in",res=400)
	par(mfrow=c(1,2),cex=.7)
	hist(cor_vec_exp_Ctmax)
	hist(cor_vec_exp_CCR)
	dev.off()'
	# I guess this makes sense, since these genes show plastic response, and therefore the ones that does not respond (that is 0 correlations) would not show up in this plot.
	# Also, the expression level of a gene could respond either negatively or positively as a function of temperature acclimation, so almost equal amounts of one or the other (pos/neg corr) should be expected.



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

#Using functions from script to calculate pop specific slopes
	print("Now the script employs functions to calculate population wise slopes, intercept and means for each gene. This may be slow.")
	
	
	trendlines_data <- lapply(1:nrow(plastic_counts), function(i) trendline_pop_specific_dat_func(i, dataset=plastic_counts, meta_dataset=meta_exp,meta_treatment_column_name="Temperature",meta_pop_column_name="population"))
	names(trendlines_data)<- rownames(plastic_counts)
	
	plastic_counts_permuted <- lapply(1:nrow(plastic_counts), function (i) sample(plastic_counts[i,]))
	plastic_counts_permuted <- as.data.frame(do.call(rbind, plastic_counts_permuted))

	trendlines_data_permuted <- lapply(1:nrow(plastic_counts), function(i) trendline_pop_specific_dat_func(i, dataset=plastic_counts_permuted, meta_dataset=meta_exp,meta_treatment_column_name="Temperature",meta_pop_column_name="population"))
	names(trendlines_data_permuted)<- rownames(plastic_counts_permuted)


################################################################################################
### Use functions to find genes with similar expression patterns to CTmax
################################################################################################

print("Now the script employs functions to find patterns in slopes similar to the pattern for CTmax.")
# look for this pattern
print(slope_vec_CTmax)
print(ctmax_trends_0diff)
print(paste0("Populations with slopes NOT DIFFERENT from zero: ",ctmax_zero_pops))

# Use NEW functions to find genes with similar expression patterns to CTmax
ctmax_fit_pattern_improved <- lapply(trendlines_data, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=4.5,	threshold_for_slope=0.02,	pop_vector_zero_trend=ctmax_zero_pops, use_THRESH_for_slopeNotZero=T ))
	# finds genes where K has slope not different from zero, and below (absolute) threshold
	# while S,O,B has slopes different from zero, but in same direction and above/below threshold
NROW(which(sapply(ctmax_fit_pattern_improved,"[")==T))  # 274 with threshold 0.2 for slopes different from zero, 320 without threshold for slopes different from zero (keeping it for slopes equal to zero)
## OBS: why when I set slope_threshold=0.03, er der 212,
#  		when I set it to threshold=0.02, er der 274
	# should it not go the other way, when I have a smaller the threshold?
	# No, because, one of the criteria is also that ALL slopes of the other pops (than K) should all be either below the -threshold or above the threshold.
		# loosening the threshold makes more genes fullfill this criteria, than it sort away based on the K slope. 
CTMAX_onlyFIT <- trendlines_data[which(sapply(ctmax_fit_pattern_improved,"[")==T)]
print(paste0(NROW(CTMAX_onlyFIT)," genes fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved)," genes."))
write.table(names(CTMAX_onlyFIT),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/CTmax_Exp_fits_slope_patterns_",padj_threshold,".tab"), row.names=F,col.names=F)
	'"274 genes fit the improved CTmax pattern (Nov22), out of 11807 genes."'


	# running with permutation to get a false positve estimate
ctmax_fit_pattern_improved_PERMUTED <- lapply(trendlines_data_permuted, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=4.5,	threshold_for_slope=0.02,	pop_vector_zero_trend=ctmax_zero_pops, use_THRESH_for_slopeNotZero=T ))
NROW(which(sapply(ctmax_fit_pattern_improved_PERMUTED,"[")==T))  # 274 with threshold 0.2 for slopes different from zero, 320 without threshold for slopes different from zero (keeping it for slopes equal to zero)
CTMAX_onlyFIT_PERMUTED <- trendlines_data_permuted[which(sapply(ctmax_fit_pattern_improved_PERMUTED,"[")==T)]
print(paste0(NROW(CTMAX_onlyFIT_PERMUTED)," PERMUTED genes fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved_PERMUTED)," genes. (N FALSE POSITIVEs)"))
	'0 PERMUTED genes fit the improved CTmax pattern (Nov22), out of 11807 genes. (N FALSE POSITIVEs)'
write.table(names(CTMAX_onlyFIT_PERMUTED),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/CTmax_Exp_fits_slope_patterns_PERMUTED_",padj_threshold,".tab"), row.names=F,col.names=F)


		######################
		##### Plot them	(only non-permuted)
		######################
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
counter=0
name_vector <- sprintf("%s_%02d", "CTmax_expression_patternsfit_nov22_slope_FIT", 1:100)
for (rows in seq(NROW(CTMAX_onlyFIT))	){
	if ( rows %in%  seq(0, NROW(CTMAX_onlyFIT), 20)  || rows == NROW(CTMAX_onlyFIT)   	){
		# if rows are 20, 40 60 etc, or if it is the last row
		if( rows != 1){ # if different from 1
			counter=counter+1
			png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/",name_vector[counter], "_", padj_threshold, ".png"),
						width=8.5,height=5.5,units="in",res=200)
				#Open picture file
			#message(rows)
			if (rows == NROW(CTMAX_onlyFIT)){
					#If last row
				#seq_start <- NROW(CTMAX_onlyFIT) - (	rows%%seq(0, NROW(CTMAX_onlyFIT), 20)[	NROW(seq(0, NROW(CTMAX_onlyFIT), 20))	]	) + 1
				seq_start <- NROW(CTMAX_onlyFIT) - (	rows%%seq(0, NROW(CTMAX_onlyFIT), 20)[	NROW(seq(0, NROW(CTMAX_onlyFIT), 20))	]	) 
			}else{
				seq_start <- seq(0, NROW(CTMAX_onlyFIT), 20) [ 	which( seq(0, NROW(CTMAX_onlyFIT), 20) == rows) 	] - 19			}
					#get starting number (rownumber to plot)
			#message(seq_start)			
			
			number_vector <- seq_start:rows
			par(mfrow=c( ceiling(NROW(number_vector)/5) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
			for(each in number_vector){
				plot_count_is=which(number_vector == each)
				plot_genes_with_pattern_improved(gene_nr=each,intensity_data=plastic_counts,meta_data_frame=meta_exp,meta_dat_column_name="Temperature",cols_vector=popcols,dataset_slopes_intercept=CTMAX_onlyFIT, plot_count=plot_count_is)			
				}
			#message(name_vector[counter]) # if I need to check numbers.
			par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
    		plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
			legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
			dev.off()
}	}	}

# if you want to check gene slopes etc. use
	# trendlines_data[[ which( names( trendlines_data ) == "g3" ) ]]


#sapply()
			#sapply(number_vector, function(i) plot_genes_with_pattern_improved(gene_nr = i,intensity_data=plastic_counts, meta_data_frame=meta_exp, meta_dat_column_name="Temperature", cols_vector=popcols, dataset_slopes_intercept=CTMAX_onlyFIT))
			#function(i) plot_genes_with_pattern_improved(intensity_data=plastic_counts, meta_data_frame=meta_exp, meta_dat_column_name="Temperature", gene_nr=i, cols_vector=popcols, dataset_slopes_intercept=CTMAX_onlyFIT))
			
################################################################################################
### Use functions to find genes with similar expression patterns to CCR
################################################################################################

print("Similar analyses are done for genes with patterns fitting the CCR pattern.")
# look for this pattern
print(slope_vec_CCR)
print(ccr_trends_0diff)
print(paste0("Populations with slopes NOT DIFFERENT from zero: ",ccr_zero_pops))

ccr_fit_pattern_improved <- lapply(trendlines_data, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=4.5,	threshold_for_slope=0.02,	pop_vector_zero_trend=ccr_zero_pops, use_THRESH_for_slopeNotZero=T ))
NROW(which(sapply(ccr_fit_pattern_improved,"[")==T))
## OBS: Play with threshold based on graphs
	# 0.02 --> 87 (as opposed to 78 with the thresholds for slopes different than zero.)
	# We keep it as the threshold for CTmax (0.02)
	
CCR_onlyFIT <- trendlines_data[which(sapply(ccr_fit_pattern_improved,"[")==T)]
print(paste0(NROW(CCR_onlyFIT)," genes fit the improved CCRTemp pattern (Nov22), out of ",NROW(ccr_fit_pattern_improved)," genes."))
write.table(names(CCR_onlyFIT),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/CCR_Exp_fits_slope_patterns_",padj_threshold,".tab"), row.names=F,col.names=F)
'"71 genes fit the improved CCRTemp pattern (Nov22), out of 11807 genes."'
	# Running with PERMUTATION to get an estimate of false positives
ccr_fit_pattern_improved_PERMUTED <- lapply(trendlines_data_permuted, function(i) fit_ctmax_improved_filter(genewise_dataframe=i,	threshold_for_mean=4.5,	threshold_for_slope=0.02,	pop_vector_zero_trend=ccr_zero_pops, use_THRESH_for_slopeNotZero=T ))
NROW(which(sapply(ccr_fit_pattern_improved_PERMUTED,"[")==T))
CCR_onlyFIT_PERMUTED <- trendlines_data_permuted[which(sapply(ccr_fit_pattern_improved_PERMUTED,"[")==T)]
print(paste0(NROW(CCR_onlyFIT_PERMUTED)," PERMUTED genes fit the improved CCRTemp pattern (Nov22), out of ",NROW(ccr_fit_pattern_improved_PERMUTED)," genes. Estimating N false positives."))
write.table(names(CCR_onlyFIT_PERMUTED),paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/CCR_Exp_fits_slope_patterns_PERMUTED_",padj_threshold,".tab"), row.names=F,col.names=F)
'"2 PERMUTED genes fit the improved CCRTemp pattern (Nov22), out of 11807 genes. Estimating N false positives."'

		######################
		##### Plot them	
		######################
print("Now the genes fitting the general pattern is plotted - so we can see if it makes sense.")
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,6,4,2)]	# B, E/O,	S, K, 
pop_vector <- c("B","O","S","K")
popcols <- popcols[order(pop_vector)]
pop_vector <- pop_vector[order(pop_vector)]
if( any((pop_vector) == "B") ) { pop_vector[which(pop_vector == "B")] <- "Betta" }
if( any((pop_vector) == "K") ) { pop_vector[which(pop_vector == "K")] <- "Karasburg" }
if( any((pop_vector) == "S") ) { pop_vector[which(pop_vector == "S")] <- "Stampriet" }
if( any((pop_vector) == "O") ) { pop_vector[which(pop_vector == "O")] <- "Otavi" }
		# loop to separate plots in files of 20 each
counter=0
name_vector <- sprintf("%s_%02d", "CCR_expression_patternsfit_nov22_slope_FIT", 1:100)
for (rows in seq(NROW(CCR_onlyFIT))	){
	if ( rows %in%  seq(0, NROW(CCR_onlyFIT), 20) || rows == NROW(CCR_onlyFIT)  	){
		# if rows are 20, 40 60 etc, or if it is the last row
		if( rows != 1){ # if different from 1
			counter=counter+1
			png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/expression_vs_temperature_tolerances/",name_vector[counter],"_",padj_threshold,".png"),
						width=8.5,height=5.5,units="in",res=400)
				#Open picture file
			message(rows)
			if (rows == NROW(CCR_onlyFIT)){
					#If last row
				seq_start <- NROW(CCR_onlyFIT) - (	rows%%seq(0, NROW(CCR_onlyFIT), 20)[	NROW(seq(0, NROW(CCR_onlyFIT), 20))	]	) + 1
			}else{
				seq_start <- seq(0, NROW(CCR_onlyFIT), 20) [ 	which( seq(0, NROW(CCR_onlyFIT), 20) == rows) 	] - 19			}
					#get starting number (rownumber to plot)
			message(seq_start)			
			
			number_vector <- seq_start:rows
			par(mfrow=c( ceiling((NROW(number_vector)/5)) , 5),		mar=c(2.1,2,2,1), oma=c(2,2.5,2,0))
			for(each in number_vector){
				plot_count_is=which(number_vector == each)
				plot_genes_with_pattern_improved(gene_nr=each,intensity_data=plastic_counts,meta_data_frame=meta_exp,meta_dat_column_name="Temperature",cols_vector=popcols,dataset_slopes_intercept=CCR_onlyFIT, plot_count=plot_count_is)			
				}
			message(name_vector[counter])
			par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
    		plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
			legend("top", legend=pop_vector, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
			dev.off()
}	}	}

print(paste0(NROW(CTMAX_onlyFIT)," genes fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved)," genes."))
print(paste0(NROW(CCR_onlyFIT)," genes fit the improved CCRTemp pattern (Nov22), out of ",NROW(ccr_fit_pattern_improved)," genes."))
print("Done")

'
	print(paste0(NROW(CTMAX_onlyFIT)," genes fit the improved CTmax pattern (Nov22), out of ",NROW(ctmax_fit_pattern_improved)," genes."))
[1] "274 genes fit the improved CTmax pattern (Nov22), out of 11807 genes."
	> print(paste0(NROW(CCR_onlyFIT)," genes fit the improved CTmax pattern (Nov22), out of ",NROW(ccr_fit_pattern_improved)," genes."))
[1] "71 genes fit the improved CTmax pattern (Nov22), out of 11807 genes."
'