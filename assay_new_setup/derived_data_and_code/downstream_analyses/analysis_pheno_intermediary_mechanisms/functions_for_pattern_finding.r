#!/usr/bin/env Rscript
##### possibly modify to go into assay_copy env.


# function for calculating population specific mean, intercept and slope, and test whether the slope is significantly different that 0
	# returns a data frame with 5 columns and the number of populations as rows.
		#columns are "Row.names" (basicallt popylations), "mean", "intercept", "slope", "p.value" of the slope.
	# inputs are
		# gene_nr (ie. row number in imput dataset, fx. in an lapply function)
		# the dataset (dataframe with samples in columns and genes in rows)
		# metadataset indicating which samples belong to which population and which treatment
		# meta_treatment_column_name = the name of the column with treatment data
		# meta_pop_column_name = the name of the column with population data
	trendline_pop_specific_dat_func <- function(gene_nr, dataset, meta_dataset, meta_treatment_column_name="temperature",meta_pop_column_name="population"){
		data_sub_gene <- t(dataset[gene_nr,])
		meta_sub_treatment_metacol <- meta_dataset[,which(colnames(meta_dataset) == meta_treatment_column_name)]
		meta_sub_pop_metacol <- meta_dataset[,which(colnames(meta_dataset) == meta_pop_column_name)]
		model_linear <- lm(data_sub_gene ~ meta_sub_treatment_metacol*meta_sub_pop_metacol)
			#summary(model_linear)
		# Get information on whether trends are different from zero
			#summary(emtrends(model_linear, ~meta_sub_pop_metacol,var="meta_sub_treatment_metacol"),infer=c(TRUE,TRUE),null=0)
		significant_trend_zero <- summary(emtrends(model_linear, ~meta_sub_pop_metacol,var="meta_sub_treatment_metacol"),infer=c(TRUE,TRUE),null=0)[,c(1,8)]
		# get pop specific intercepts (and slopes, since they are similar to emtrends slopes)
		for(pop_is in sort(levels(as.factor(meta_sub_pop_metacol)))){
			pop_abLine_values<- lm(data_sub_gene[which(meta_sub_pop_metacol==pop_is),]~meta_sub_treatment_metacol[which(meta_sub_pop_metacol==pop_is)])$coefficient
			if( sort(levels(as.factor(meta_sub_pop_metacol)))[1] == pop_is ){
				pop_vec <- pop_is
				mean_calc <- mean(data_sub_gene[which(meta_sub_pop_metacol==pop_is),], na.rm=T)
				intercept_df <- data.frame(mean_calc,pop_abLine_values[1],pop_abLine_values[2])
			}else{
				pop_vec <- append(pop_vec,pop_is)
				mean_calc <- mean(data_sub_gene[which(meta_sub_pop_metacol==pop_is),], na.rm=T)
				intercept_df <- rbind(intercept_df,c(mean_calc, pop_abLine_values[1],pop_abLine_values[2]))
			}
		}
		rownames(intercept_df) <- pop_vec
		colnames(intercept_df) <- c("mean","intercept","slope")
		# Merge intercept, trends and whether trends are significant
		trendlines_data <- merge(intercept_df,significant_trend_zero,by.x=0,by.y=1)
		return(trendlines_data)
		}


# smart function for negating "found in" (%in%), thus not found in
'%!in%' <- function(x,y)!('%in%'(x,y))

	pop_subset_counts_improved <- function(intensity_data, meta_data_frame, meta_dat_column_name="temperature", pop, pop_nr, gene_name, jitter_move_vec){ 
			gene_number_int <- grep(paste0(gene_name,"$"),rownames(intensity_data))
			data_subset <- intensity_data[	gene_number_int	,	which(meta_data_frame$population == pop)	]
			treat_subset <- meta_data_frame[	which(meta_data_frame$population == pop)	,	which(colnames(meta_data_frame) == meta_dat_column_name)	]
			treat_subset <- treat_subset+jitter_move_vec[pop_nr]
			data.frame(treat_subset,t(data_subset))
		}


		
	plot_genes_with_pattern_improved <- function(gene_nr, intensity_data, meta_data_frame, meta_dat_column_name="temperature", cols_vector, dataset_slopes_intercept,ylab_text="Expression\n(VST scaled counts)",plot_count){	
		jitter_move=seq(-.5,.5,length.out=4)
		xlab_text <- expression(paste("Acclimation Temperature (", degree,"C)"))
		#get gene name from number in dataset_slipes_intercept
		gene_name_is <- names(dataset_slopes_intercept[gene_nr])
		#get pop column
		pop_vector <- dataset_slopes_intercept[[gene_nr]]$Row.names
		par_mfrowSets <- par()$mfrow
		if( par_mfrowSets[1] == 1){	# if there is only one row.
			xlab_is <- xlab_text		# all will have xlab
			if ( plot_count == 1 ) {		# if plot is in the first column
				ylab_is <- ylab_text
			}else{
				ylab_is <- NA
			}
		}else if( par_mfrowSets[1] > 1 ) {# if there are multiple rows
			last_row_plot_numbers <- seq( (par_mfrowSets[1]-1)*par_mfrowSets[2]+1, par_mfrowSets[1]*par_mfrowSets[2] )
			if ( plot_count %in% c(1,6,11,16,21) ) {		# if plot is in the first row.
				ylab_is <- ylab_text
			}else{			# if plot is not in first row --> no ylab
				ylab_is <- NA
			}

			if ( plot_count %in% last_row_plot_numbers ) {		# if plot is in the last row.
				xlab_is <- xlab_text
			}else{
				xlab_is <- NA
			}
		}
		for( pop in sort(pop_vector) ){
				#subset gene intensity(expression data)
			intensity_sub_gene <- intensity_data[which(rownames(intensity_data)==gene_name_is),]
				#set y limits
			ylims <- c(min(intensity_sub_gene,na.rm=T),max(intensity_sub_gene,na.rm=T))
			pop_number <- which(pop_vector == pop)
					
			subset_is <- pop_subset_counts_improved(intensity_data=intensity_data, meta_data_frame=meta_data_frame, meta_dat_column_name=meta_dat_column_name, pop=pop, pop_nr=pop_number, gene_name=gene_name_is, jitter_move_vec=jitter_move)
				 # looks like this: RNA (rownames), x_value (temp+jitter), y_value (expression count vstd)
			
			slopes_sub <- dataset_slopes_intercept[[gene_name_is]]
			interceptee	<- 	slopes_sub[	which(slopes_sub$Row.names == pop)	,	]$intercept
			slopee	<- 	slopes_sub[	which(slopes_sub$Row.names == pop)	,	]$slope

			message(gene_name_is, ":    ",pop,": ",interceptee," + ",slopee)
			
			if ( pop_number == 1 ){
				par(new=F, cex.main=.9)
				plot(jitter(subset_is[,1],factor=.15),subset_is[,2], ylab=NA, col=cols_vector[pop_number],
					xlab=NA, main=colnames(subset_is)[2], xaxt="n",xpd=NA,	xlim=c(14.5,29.5), ylim=ylims)
				axis(1, at=c(15,17,19,21,23,25,27,29), labels=NA, xpd=NA, col.tick="grey75") # Add grey ticks
				axis(1, at=c(15, 19, 23, 29), labels=c(15,19,23,29), xpd=NA)
				axis(1, at=c(25), labels=c(25), xpd=NA)
				abline(interceptee,slopee, col=cols_vector[pop_number])
				title(xlab=xlab_is, xpd=NA)
				title(ylab=ylab_is, xpd=NA, line=2.5)
				if( plot_count == 1){
					#legend_y <- max(subset_is[,2],na.rm=T) + (( max(subset_is[,2],na.rm=T)	- min(subset_is[,2],na.rm=T)	)	*1.3 )
					#legend(35,legend_y, legend=pop_vector, lty=1, lwd=1.1, col=cols_vector, xpd=NA, bty="n",horiz=T, cex=1.1)
				}
				
			}else{
				par(new=T)
				plot(jitter(subset_is[,1],factor=.15),subset_is[,2], col=cols_vector[pop_number],xlab=NA,ylab=NA,
					xaxt="n",yaxt="n",xlim=c(14.5,29.5), ylim=ylims)
				abline(interceptee,slopee, col=cols_vector[pop_number])
			}
		}
	}
	
# function for getting the genes that share patterns with CTmax:
	# returns whether all the filters are TRUE or FALSE for the given dataset (genewise)
fit_ctmax_improved_filter <- function(genewise_dataframe, threshold_for_mean=4.5, threshold_for_slope=NULL,pop_vector_zero_trend=ctmax_zero_pops, use_THRESH_for_slopeNotZero=T){
			#fix used subsets/variables
				# make population vector with slopes not different from zero into 1 letter vector(factor)
		pops_zero_trend_1Letter <- substr(pop_vector_zero_trend,1,1)
		# main part of function
		truth_vector <- vector()

			# all means above threshold	
		if( is.numeric(threshold_for_mean) ){
			truth_vector <- append(truth_vector,all(genewise_dataframe$mean > threshold_for_mean))
		}
			# all pops in vector must have slopes not different from zero: (in this case only karasburg should have pvalue above 0.05)

		truth_vector <- append(truth_vector,	all(genewise_dataframe[which(genewise_dataframe$Row.names %in% pops_zero_trend_1Letter),]$p.value	> 0.05) 	)
				# and absolute value of slope must also be below threshold for zero slopes
		if( is.numeric(threshold_for_slope) ){
				truth_vector <- append(truth_vector,	all(abs(genewise_dataframe[which(genewise_dataframe$Row.names %in% pops_zero_trend_1Letter),]$slope)	< threshold_for_slope)	) 
		} #&&
			# all pops not in the vector mentioned above, must have slopes different from zero (ergo pvalue below 0.05)
		truth_vector <- append(truth_vector,	all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$p.value	< 0.05) 	)#&&

		#must be either:
			# all population with a slope different from 0, should have a positive slope above set threshold
		if( is.numeric(threshold_for_slope) ){
			if( use_THRESH_for_slopeNotZero ){
				truth_vector <- append(truth_vector,	(	all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	> threshold_for_slope) ||
															all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	< (-threshold_for_slope)) )
				)
			}else if( use_THRESH_for_slopeNotZero == F ){
				truth_vector <- append(truth_vector,	(	all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	> 0) ||
														all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	< 0) )
				)
			}
		}
		if(	length(truth_vector)>0 ){
			all(truth_vector) 
		}
	}





fit_CCR_improved_filter <- function(genewise_dataframe, threshold_for_mean=4.5, threshold_for_slope=NULL,pop_vector_zero_trend=ccr_zero_pops, use_THRESH_for_slopeNotZero=T){
			#fix used subsets/variables
				# make population vector with slopes not different from zero into 1 letter vector(factor)
		pops_zero_trend_1Letter <- substr(pop_vector_zero_trend,1,1)
		# main part of function
		truth_vector <- vector()		# append all logical end-statements to vector to evaluate in the end.
			# all means above threshold	
		if (is.numeric(threshold_for_mean)){
			truth_vector <- append(truth_vector,		all(genewise_dataframe$mean > threshold_for_mean) 	)
		}
			# all pops in vector must have slopes equal to zero: (in this case only karasburg should have pvalue above 0.05)
		truth_vector <- append(truth_vector,		all(genewise_dataframe[which(genewise_dataframe$Row.names %in% pops_zero_trend_1Letter),]$p.value	> 0.05) 	)#&&
			# and absolute value of slope must also be below threshold for zero slopes
		if (is.numeric(threshold_for_slope)){
			truth_vector <- append(truth_vector,		all(abs(genewise_dataframe[which(genewise_dataframe$Row.names %in% pops_zero_trend_1Letter),]$slope)	< threshold_for_slope) 	)
		} #&&
			# all pops not in the vector mentioned above, must have slopes different from zero (ergo pvalue below 0.05)
		truth_vector <- append(truth_vector,		all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$p.value	< 0.05) 	)#&&

		#must be either:
			# all population with a slope different from 0, should have a positive slope above set threshold
		if (is.numeric(threshold_for_slope)){
			if( use_THRESH_for_slopeNotZero ){
				truth_vector <- append(truth_vector,	(	all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	> threshold_for_slope) ||
															all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	< -(threshold_for_slope))	)		
								)
			}else if( use_THRESH_for_slopeNotZero == F ){
				truth_vector <- append(truth_vector,	(	all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	> 0) ||
															all(genewise_dataframe[which(genewise_dataframe$Row.names %!in% pops_zero_trend_1Letter),]$slope	< 0)	)		
								)
			}
		}
		if(	length(truth_vector)>0 ){
			all(truth_vector)
		}
}




























# function for returning only slope from lm function
	#sub_dat = dframe subset to contain samples from a specific population and a specific gene(feature)
	# sub_dat_treatment= a vector of the corresponding treatment temperature
	# returns the slope.
	lm_slope <- function(sub_dat,sub_dat_treatment){ ## lm_slope(subset_dat$ToD_mean, subset_dat$Treatment)
		lmis <- lm(sub_dat ~ sub_dat_treatment)
		return(lmis$coefficients[2])
	}
	
# function for returning only intercept from lm function
	lm_intercept <- function(sub_dat,sub_dat_treatment){ ## lm_slope(subset_dat$ToD_mean, subset_dat$Treatment)
		lmis <- lm(sub_dat ~ sub_dat_treatment)
		return(lmis$coefficients[1])
	}

# function for calculating population specific slope in expression data
	# using above function (this function will be used in an apply setting)
	# pop= a better indicating the pop
	# gene_nr= the gene/feature number to look at the slope for
	# dataset dataframe with samples in columns and genes in rows
	# metadataset indicating which samples belong to which population, samples in rows, info in columns. (with col called "population")
	# meta_column_name= meta data column to calculate the slope based on. Must be numeric? default is "temperature"
	pop_spec_lm_slope <- function(pop,gene_nr, dataset, meta_dataset, meta_column_name="temperature"){ ## lm_slope(subset_dat$ToD_mean, subset_dat$Treatment)
		data_subset <- dataset[,which(meta_dataset$population == pop)]
		#gene_is <- rownames(data_subset[gene_nr,])
		treat_subset <- meta_dataset[which(meta_dataset$population == pop),which(colnames(meta_dataset) == meta_column_name)]
		lm_slope(t(data_subset[gene_nr,]),treat_subset)
	}	# for metabolites dataset must be transformed (t()) OBS!
	

# function for calculating population specific intercept using above function (this function will be used in an apply setting)
	pop_spec_lm_intercept <- function(pop,gene_nr, dataset, meta_dataset, meta_column_name="temperature"){ ## lm_slope(subset_dat$ToD_mean, subset_dat$Treatment)
		data_subset <- dataset[,which(meta_dataset$population == pop)]
		#gene_is <- rownames(data_subset[gene_nr,])
		treat_subset <- meta_dataset[which(meta_dataset$population == pop),		which(colnames(meta_dataset) == meta_column_name)]
		lm_intercept(t(data_subset[gene_nr,]),treat_subset)
	}
	
# function for calculating population specific mean (this function will be used in an apply setting)
	pop_spec_mean <- function(pop,gene_nr, dataset, meta_dataset){ 
		data_subset <- dataset[,which(meta_dataset$population == pop)]
		mean(t(data_subset[gene_nr,]),na.rm=T)
	}


# FROM GO SCRIPT - function for subsetting data for plotting, including jitter_move
	pop_subset_counts <- function(intensity_data, meta_data_frame, meta_dat_column_name="temperature", pop, pop_nr, gene_name, jitter_move_vec){ ## lm_slope(subset_dat$ToD_mean, subset_dat$Treatment)
			gene_nr <- grep(paste0(gene_name,"$"),rownames(intensity_data))
			data_subset <- intensity_data[gene_nr,which(meta_data_frame$population == pop)]
			#gene_is <- rownames(data_subset[gene_nr,])
			treat_subset <- meta_data_frame[which(meta_data_frame$population == pop),which(colnames(meta_data_frame) == meta_dat_column_name)]
			treat_subset <- treat_subset+jitter_move_vec[pop_nr]
			#message(treat_subset)
			data.frame(treat_subset,t(data_subset))
		}
	# FROM GO script
	plot_genes_with_pattern <- function(intensity_data, meta_data_frame, meta_dat_column_name="temperature", gene_nr, pop_vector, cols_vector, dataset_slopes_intercept,ylab_text="Expression (VST scaled counts)"){	
		jitter_move=seq(-.5,.5,length.out=4)
		for( pop in pop_vector ){
			rownames(dataset_slopes_intercept) <- dataset_slopes_intercept[,1]
			pop_number <- which(pop_vector == pop)
			gene_name <- rownames(dataset_slopes_intercept)[gene_nr]
			ylims <- c(min(intensity_data[which(rownames(intensity_data)==gene_name),],na.rm=T),max(intensity_data[which(rownames(intensity_data)==gene_name),],na.rm=T))
			subset_is <- pop_subset_counts(intensity_data=intensity_data, meta_data_frame=meta_data_frame, meta_dat_column_name=meta_dat_column_name, pop=pop, pop_nr=pop_number, gene_name=gene_name, jitter_move_vec=jitter_move)
			#print(subset_is)
			#print(ylims)
				 # looks like this: RNA (rownames), x_value (temp+jitter), y_value (expression count vstd)
			interceptee <- dataset_slopes_intercept[colnames(subset_is)[2],	grep(paste0("^",pop,".*","intercept","$"),colnames(dataset_slopes_intercept))]	#get intercept
			slopee <- dataset_slopes_intercept[colnames(subset_is)[2],	grep(paste0("^",pop,".*","slopes","$"),colnames(dataset_slopes_intercept))]	#get intercept	#get slope
			message(gene_name, ":    ",pop,": ",interceptee," + ",slopee)
			if (pop_number==1){
				par(new=F)
				plot(jitter(subset_is[,1],factor=.15),subset_is[,2],ylab=ylab_text, col=cols_vector[pop_number],
					xlab="Temperature acclimation", main=colnames(subset_is)[2], xaxt="n",xpd=NA,	xlim=c(0,29.5), ylim=ylims) #xlim=c(14.5,29.5))
				axis(1, at=c(0,15,17,19,21,23,25,27,29), labels=c(0,15,NA,19,NA,23,25,NA,29), xpd=NA)
				abline(interceptee,slopee, col=cols_vector[pop_number])
			}else{
				par(new=T)
				plot(jitter(subset_is[,1],factor=.15),subset_is[,2], col=cols_vector[pop_number],xlab=NA,ylab=NA,xaxt="n",yaxt="n",xlim=c(0,29.5), ylim=ylims)#xlim=c(14.5,29.5)) #axes=F
				abline(interceptee,slopee, col=cols_vector[pop_number])
			}
		}
	}
	


	#returns a named vector of significance letters or number is I like. 
	# it is made to be used across each row of the dataframe given to the function (using lapply)
gene_wise_model_methylation_means <- function(row_nr, dataframe=count_meth){
	each_gene_sub <-  gather(count_meth[row_nr,]) # make into long format
	each_gene_sub$Temperature <- sapply(strsplit(each_gene_sub$key,"([A-Z])"), "[[",2)
	each_gene_sub$Population <- sapply(strsplit(each_gene_sub$key,"([0-9])",perl=T), "[[",1)
	model <- aov(lm(each_gene_sub$value ~ each_gene_sub$Population))
	vector_letters <- cld(emmeans(model, "Population"),Letters=letters)$.group
	#vector_letters <- cld(emmeans(model, "Population"))$.group
	names(vector_letters) <- cld(emmeans(model, "Population"))$Population
	return(vector_letters)
}


	# genewise dataframe supplied by lapply function?
	# takes the significance letter display for a named vector (named by population/sample) and compares it to a named "trend vector" of significance letters.
	# returns true if all elements are similar
fit_ctmax_improved_filter_methylation <- function(genewise_vector_methylation_means, ctmax_trends=ctmax_trends_diff_vec){
		#genewise_vector_methylation_means=sig_dif_means_meth[[i]]
			#fix used subsets/variables
				# make population vector with slopes not different from zero into 1 letter vector(factor)
		veridity_vec=vector()
		genewise_vector_methylation_means <- gsub(" ", "", genewise_vector_methylation_means, fixed = TRUE)
		names_trend <- names(ctmax_trends)
		ctmax_trends <- gsub(" ", "", as.vector(ctmax_trends),fixed=T)
		names(ctmax_trends) <- names_trend	;	rm(names_trend)
		#ctmax_trends <- gsub(" ", "", ctmax_trends, fixed = TRUE)
		
		for( each_name in names(genewise_vector_methylation_means)){
			#each_name <- names(genewise_vector_methylation_means[1])
			meth_sub <- genewise_vector_methylation_means[get("each_name")]
			temp_nr <- grep(each_name,names(ctmax_trends))
			# got to here:
			# subset ctmax_trends to temp_nr and check if similar to each_name
			ctmax_trends_sub <- ctmax_trends[temp_nr]
			veridity_vec <- append(veridity_vec,	ctmax_trends_sub == meth_sub	)
			#if true, next iteration
			# if false .... no make vector of true false and check all in the end.
		}

		if(	length(veridity_vec)>0 ){
			all(veridity_vec) 
		}
		
}



fit_ccr_improved_filter_methylation <- function(genewise_vector_methylation_means, ctmax_trends=ccr_means_diff_vec){
		#genewise_vector_methylation_means=sig_dif_means_meth[[i]]
			#fix used subsets/variables
				# make population vector with slopes not different from zero into 1 letter vector(factor)
		veridity_vec=vector()
		genewise_vector_methylation_means <- gsub(" ", "", genewise_vector_methylation_means, fixed = TRUE)
		names_trend <- names(ctmax_trends)
		ctmax_trends <- gsub(" ", "", as.vector(ctmax_trends),fixed=T)
		names(ctmax_trends) <- names_trend	;	rm(names_trend)
		#ctmax_trends <- gsub(" ", "", ctmax_trends, fixed = TRUE)
		
		print(levels(as.factor(ctmax_trends)))
		ctmax_trends_names_1 <- names(ctmax_trends)[grep(levels(as.factor(ctmax_trends))[1],ctmax_trends)]
		message("Level 1: ", levels(as.factor(ctmax_trends))[1], ".  Populations: ",		ctmax_trends_names_1 )
		if( 	NROW(levels(as.factor(ctmax_trends))) > 1 ){
			if( 	NROW(levels(as.factor(ctmax_trends))) ==3 ){
				number_for_other_level <- which(	nchar((levels(as.factor(ctmax_trends)))) < 2	)
				number_for_other_level <- number_for_other_level[NROW(number_for_other_level)]
				ctmax_trends_names_2 <- names(ctmax_trends)[grep(levels(as.factor(ctmax_trends))[number_for_other_level],ctmax_trends)]
				message("Level 2: ", levels(as.factor(ctmax_trends))[number_for_other_level], ".  Populations: ",		ctmax_trends_names_2 )
			}
		}
		# I am trying to figure out a way to see if a vector contain the same letter (any letter) in all instances

			#check whether a, b or c is found in all instances in vector for populations that should have same letter 
		if( all(grepl("a",genewise_vector_methylation_means[ctmax_trends_names_1])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("b",genewise_vector_methylation_means[ctmax_trends_names_1])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("c",genewise_vector_methylation_means[ctmax_trends_names_1])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("d",genewise_vector_methylation_means[ctmax_trends_names_1])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else{veridity_vec <- append(veridity_vec,FALSE)}
			#Check the sqame for vector 2
		if( all(grepl("a",genewise_vector_methylation_means[ctmax_trends_names_2])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("b",genewise_vector_methylation_means[ctmax_trends_names_2])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("c",genewise_vector_methylation_means[ctmax_trends_names_2])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else if( all(grepl("d",genewise_vector_methylation_means[ctmax_trends_names_2])) ){	veridity_vec <- append(veridity_vec,TRUE)	
		}else{veridity_vec <- append(veridity_vec,FALSE)}

		if( NROW(ctmax_trends_names_1)<4 ){
			if( all(grepl("a",genewise_vector_methylation_means)) ){	veridity_vec <- append(veridity_vec,FALSE)	
			}else if( all(grepl("b",genewise_vector_methylation_means)) ){	veridity_vec <- append(veridity_vec,FALSE)	
			}else if( all(grepl("c",genewise_vector_methylation_means)) ){	veridity_vec <- append(veridity_vec,FALSE)	
			}else if( all(grepl("d",genewise_vector_methylation_means)) ){	veridity_vec <- append(veridity_vec,FALSE)	
			}
		}
		if( any(nchar(ctmax_trends) > 1)){
			# can only handle one with multiple letters assignmed at the moment.
			pop_is_now <- names(which(nchar(ctmax_trends) ==2))
			if( any(nchar(genewise_vector_methylation_means) > 1) ){
				if( table(nchar(genewise_vector_methylation_means) )[2] == 2 ){ veridity_vec <- append(veridity_vec,FALSE) }			
				
					# get this only if the population is the same as pop_is_now
				population_from_data_nchar <- names(which(nchar(genewise_vector_methylation_means)>1) )
				if( table(nchar(genewise_vector_methylation_means) )[2] == 1 && pop_is_now==population_from_data_nchar ){ veridity_vec <- append(veridity_vec,TRUE) }			
				if( table(nchar(genewise_vector_methylation_means) )[2] == 1 && pop_is_now!=population_from_data_nchar ){ veridity_vec <- append(veridity_vec,FALSE) }			
			}

		#for( each_name in names(genewise_vector_methylation_means)){
				#each_name <- names(genewise_vector_methylation_means[1])
			#meth_sub <- genewise_vector_methylation_means[get("each_name")]
			#temp_nr <- grep(each_name,names(ctmax_trends))
				# got to here:
				# subset ctmax_trends to temp_nr and check if similar to each_name
			#ctmax_trends_sub <- ctmax_trends[temp_nr]
			#veridity_vec <- append(veridity_vec,	ctmax_trends_sub == meth_sub	)
				#if true, next iteration
				# if false .... no make vector of true false and check all in the end.
		#}

		if(	length(veridity_vec)>0 ){
			all(veridity_vec) 
		}
		}
}



# function filters (methylated) genes such that pattern resembles that of heat tolerance
		# K slope is equal to zero (pattern here is: same significance group)
		# S, B and O are different from zero (pattern here is: different significance group from K, but similar direction)
	# K must be in significance group consisting of one charachter. 
	# SBO must be in another group, but both higher or lower than BO. (but groups may differ)
		# could be tightened, allowing S to only be category next to K, but BO can be any category equal to or higher/lower than S.
			# Would remove less than 10 genes (estimated by visual inspection)
fit_CTmax1_improved_filter_methylation <- function(genewise_vector_methylation_means, ctmax_trends=ctmax_zero_pops){
		# ctmax_trends=ctmax_trends_diff_vec
		# ctmax_trends - trend different from zero

		#genewise_vector_methylation_means=sig_dif_means_meth[[i]]
		pops_zero_trend_1Letter <- substr(ctmax_trends,1,1)
			# S større end K, men B og O i samme gruppe

				# make population vector with slopes not different from zero into 1 letter vector(factor)
		veridity_vec=vector()
		genewise_vector_methylation_means <- gsub(" ", "", genewise_vector_methylation_means, fixed = TRUE)

		names_diff_zero <- substr(	levels(ctmax_trends)[ which(levels(ctmax_trends) != ctmax_trends) ]	,1,1	)

		sign_letter_zeropops <- genewise_vector_methylation_means[which(names(genewise_vector_methylation_means) == pops_zero_trend_1Letter)]
		sign_letter_NONzeropops <- genewise_vector_methylation_means[which(names(genewise_vector_methylation_means) != pops_zero_trend_1Letter)]		

		if( NROW(unique(genewise_vector_methylation_means)) == 1){ 		veridity_vec <- append(veridity_vec,FALSE)		}

		if( NROW(unique(sign_letter_zeropops)) == 1){ # if there is just one significance group in zero-vector
			if( nchar(unique(sign_letter_zeropops)) == 1 ){ # There must only be one significance letter - not two+
				if( all( sign_letter_zeropops != sign_letter_NONzeropops ) ) { # all from not-zero group must be different from zero-group
					if( any(grepl(unique(sign_letter_zeropops), sign_letter_NONzeropops)) ){ # if the significance letter is found in non-zero group, eval false
						veridity_vec <- append(veridity_vec,FALSE)		
					}else{ # in these statements, it could be added that either all of second vector must be similar, or stampriet must be "smaller" than O and B
						#if zero letter is a, all others must be one of b, c or d
						if( unique(sign_letter_zeropops)=="a" ){
							if( all(	grepl( "b|c|d",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else{veridity_vec <- append(veridity_vec,FALSE)}
						}else if( unique(sign_letter_zeropops)=="b" ){
								#if zero letter is b, all others must all be one of either a or c
							if( all(	grepl( "a",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else if( all(	grepl( "c",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else{veridity_vec <- append(veridity_vec,FALSE)}
						}else if( unique(sign_letter_zeropops)=="c" ){
								#if zero letter is c, all others must be one of a or b
							if( all(	grepl( "a|b",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else{veridity_vec <- append(veridity_vec,FALSE)}
						}else if( unique(sign_letter_zeropops)=="d" ){
							#if zero letter is d, all others must be one of a, b or c
							if( all(	grepl( "a|b|c",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else{veridity_vec <- append(veridity_vec,FALSE)}
						}else{veridity_vec <- append(veridity_vec,FALSE)}
						#sign_letter_NONzeropops>4
					}
					
				}else{ 	veridity_vec <- append(veridity_vec,FALSE)		}
			}else{		veridity_vec <- append(veridity_vec,FALSE)		}
		}else{ # if multiple sign. groups in zero vector eval false
			veridity_vec <- append(veridity_vec,FALSE)	
		}
		veridity_vec <- append(veridity_vec,TRUE)
		if(	length(veridity_vec)>0 ){ #all needs to be true to pass filter
			all(veridity_vec) 
		}
	}

`%!in%` = Negate(`%in%`)


# function filters (methylated) genes such that pattern resembles that of cold tolerance
		# B and O are equal to zero (pattern here is: in same significance group)
		# K and S are different from zero (pattern here is: different significance group from BO, but similar direction)
	# B and O must be in same significance group. 
	# KS must be in another group, but both higher or lower than BO.
	# could be loosened to encapsulate that BO may be in two different, but both higher/lower than KS?
fit_CCR1_improved_filter_methylation <- function(genewise_vector_methylation_means, ctmax_trends=ccr_zero_pops){
		# ctmax_trends=ctmax_trends_diff_vec
		# ctmax_trends - trend different from zero

		#genewise_vector_methylation_means=sig_dif_means_meth[[i]]
		pops_zero_trend_1Letter <- substr(ctmax_trends,1,1)
			# S større end K, men B og O i samme gruppe

				# make population vector with slopes not different from zero into 1 letter vector(factor)
		veridity_vec=vector()
		genewise_vector_methylation_means <- gsub(" ", "", genewise_vector_methylation_means, fixed = TRUE)

		names_diff_zero <- substr(	levels(ctmax_trends)[ which(levels(ctmax_trends) != ctmax_trends) ]	,1,1	)

		sign_letter_zeropops <- genewise_vector_methylation_means[which(names(genewise_vector_methylation_means) %in% pops_zero_trend_1Letter)]
		sign_letter_NONzeropops <- genewise_vector_methylation_means[which(names(genewise_vector_methylation_means) %!in% pops_zero_trend_1Letter)]		

		if( NROW(unique(genewise_vector_methylation_means)) == 1){ 		veridity_vec <- append(veridity_vec,FALSE)		}

		if( NROW(unique(sign_letter_zeropops)) == 1){ # if there is just one significance group in zero-vector
			if( nchar(unique(sign_letter_zeropops)) == 1 ){ # There must only be one significance letter - not two+
				if( all( sign_letter_zeropops != sign_letter_NONzeropops ) ) { # all from not-zero group must be different from zero-group
					if( any(grepl(unique(sign_letter_zeropops), sign_letter_NONzeropops)) ){ # if the significance letter is found in non-zero group, eval false
						veridity_vec <- append(veridity_vec,FALSE)		
					}else{ # in these statements, it could be added that either all of second vector must be similar, or stampriet must be "smaller" than O and B
						
						if( unique(sign_letter_zeropops)=="a" ){
								#if zero letter is a, all others must be one of b, c or d
							if( all(	grepl( "b|c|d",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		}
						}else if( unique(sign_letter_zeropops)=="b" ){
								#if zero letter is b, all others must all be one of either a or c
							if( all(	grepl( "a",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		
							}else if( all(	grepl( "c",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		}
						}else if( unique(sign_letter_zeropops)=="c" ){
								#if zero letter is c, all others must be one of a or b
							if( all(	grepl( "a|b",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		}
						}else if( unique(sign_letter_zeropops)=="d" ){
								#if zero letter is d, all others must be one of a, b or c
							if( all(	grepl( "a|b|c",sign_letter_NONzeropops ) ) ){	veridity_vec <- append(veridity_vec,TRUE)		}
						}else{veridity_vec <- append(veridity_vec,FALSE)}
						#sign_letter_NONzeropops>4
					}
					
				}else{ 	veridity_vec <- append(veridity_vec,FALSE)		}
			}else{		veridity_vec <- append(veridity_vec,FALSE)		}
		}else{ # if multiple sign. groups in zero vector eval false
			veridity_vec <- append(veridity_vec,FALSE)	
		}
		veridity_vec <- append(veridity_vec,TRUE)	
		if(	length(veridity_vec)>0 ){ #all needs to be true to pass filter
			all(veridity_vec) 
		}
	}

