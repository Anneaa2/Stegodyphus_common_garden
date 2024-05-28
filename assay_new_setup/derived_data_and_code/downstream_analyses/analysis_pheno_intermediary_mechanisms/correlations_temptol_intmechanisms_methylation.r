#!/usr/bin/env Rscript

##############################################################################
############ 		Correlation heat tol vs. plasticity in expression
##############################################################################



################################################
### Look at patterns similar to  phenotypes  ###
################################################
# 
#Using functions from script to calculate pop specific means and return list of named vectors with significant difference indication
	print("Now the script employs functions to calculate population wise means for each gene. This may be slow.")
	
count_meth <- count_meth[-which(rowMeans(count_meth,na.rm=T) < .01),]
	### OBS REMOVING genes with mean (across all pops) below 0.01 wmlvl
sig_dif_means_meth <- lapply(1:nrow(count_meth), function(i) gene_wise_model_methylation_means(i) )
names(sig_dif_means_meth)<- rownames(count_meth)

	
################################################################################################
### Use functions to find genes with similar expression patterns to CTmax
################################################################################################

print("Now the script employs functions to find patterns in slopes similar to the pattern for CTmax.")
# look for this pattern
print(slope_vec_CTmax)
print(ctmax_trends_diff)
	# create named vector to find pattern
	ctmax_trends_diff_vec <- ctmax_trends_diff$.group
	names(ctmax_trends_diff_vec) <- as.character(ctmax_trends_diff$Population)
	names(ctmax_trends_diff_vec) <- substr(names(ctmax_trends_diff_vec),1,1)
print(ctmax_means_diff)
#print(paste0("Populations with slopes NOT DIFFERENT from zero: ",ctmax_zero_pops))


#results_list <- lapply(1:NROW(sig_dif_means_meth),function(i) fit_ctmax_improved_filter_methylation(sig_dif_means_meth[[i]]))
results_list <- lapply(1:NROW(sig_dif_means_meth),function(i) fit_CTmax1_improved_filter_methylation(sig_dif_means_meth[[i]]))
names(results_list) <- names(sig_dif_means_meth)
NROW(which(results_list == TRUE))
# [1] 152

#for( gene in names(sig_dif_means_meth)){#
#	gene_df <- sig_dif_means_meth[[ which( names(sig_dif_means_meth) == gene ) ]]
#}


# Plot them
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,6,4,2)]	# B, E/O,	S, K, 
pop_vector <- c("B","O","S","K")
popcols <- popcols[order(pop_vector)]
pop_vector <- pop_vector[order(pop_vector)]

veridity_names <- names(which(results_list == T))
count_meth_sub <- count_meth[ which(rownames(count_meth) %in% veridity_names), ]


nrow(count_meth_sub)/20
name_vector <- sprintf("%s_%02d", paste0("CTmax_methylation_patternsfit_dec22_slope_ALL"), 1:30)
counter=0
for( gene in rownames(count_meth_sub)){
	#gene=rownames(count_meth_sub)[1]
	gene_nr <- which(rownames(count_meth_sub) == gene)
	count_meth_sub_sub <- count_meth_sub[which(rownames(count_meth_sub) == gene),]
	print(count_meth_sub_sub)
	jitter_move=seq(-.5,.5,length.out=4)
	ylims <- c( min(count_meth_sub_sub,na.rm=T), max(count_meth_sub_sub,na.rm=T))
	temp_seq <- as.numeric( unique(sapply(strsplit( names(count_meth_sub_sub), "[A-Z]"),"[[",2)) )

	if( gene_nr %in% seq(1,nrow(count_meth_sub),20)){
		counter=counter+1
		png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/methylation_vs_temperature_tolerances/",name_vector[counter],".png"),
							width=10,height=6.5,units="in",res=200)
		par(mfrow=c(4,5), mar=c(3.6,3.5,2.1,0.4), oma=c(0,0,2,.5), mgp=c(2,.7,0))
		pop_seq <- unique(sapply(strsplit( names(count_meth_sub_sub), "[1-9]"),"[[",1))	
		for( popsis in pop_seq ){
			pop_nr <- which(pop_seq == popsis)
			# subset to pop
			pop_sub <- count_meth_sub_sub[ grep(popsis,names(count_meth_sub_sub)) ]
			if( popsis == pop_seq[1] ){	
				par(new=F)		;		ylabis="WMLvl"		;		xlabis="Temperature acclimation"	;		mainis=rownames(pop_sub)
			}else{		
				par(new=T)		;		ylabis=NA		;		xlabis=NA	;		mainis=NA
			}
			plot(temp_seq+jitter_move[pop_nr], pop_sub,		ylab=ylabis, col=popcols[pop_nr],	xlab=xlabis, main=mainis, xaxt="n",yaxt="n",xpd=NA,	xlim=c(14.5,29.5), ylim=ylims) #xlim=c(14.5,29.5))
			if( popsis == pop_seq[1] ){
				axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)	
				axis(2, ylim=ylims, xpd=NA)	}
			abline(h=rowMeans(pop_sub,na.rm=T), col=popcols[pop_nr])
			#print(rowMeans(pop_sub,na.rm=T))
		}
		#legend()
	}else{
		for( popsis in pop_seq ){
			pop_nr <- which(pop_seq == popsis)
			# subset to pop
			pop_sub <- count_meth_sub_sub[ grep(popsis,names(count_meth_sub_sub)) ]
			if( popsis == pop_seq[1] ){	
				par(new=F)		;		ylabis="WMLvl"		;		xlabis="Temperature acclimation"	;		mainis=rownames(pop_sub)
			}else{		
				par(new=T)		;		ylabis=NA		;		xlabis=NA	;		mainis=NA
			}
			plot(temp_seq+jitter_move[pop_nr], pop_sub,		ylab=ylabis, col=popcols[pop_nr],	xlab=xlabis, main=mainis, xaxt="n",yaxt="n", xpd=NA,	xlim=c(14.5,29.5), ylim=ylims) #xlim=c(14.5,29.5))
			if( popsis == pop_seq[1] ){	
				axis(2, ylim=ylims, xpd=NA)	
				axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)	}
			abline(h=rowMeans(pop_sub,na.rm=T), col=popcols[pop_nr])
			#print(rowMeans(pop_sub,na.rm=T))
		}
		if( gene_nr %in% (seq(1,nrow(count_meth_sub),20)-1) || gene_nr==nrow(count_meth_sub) ){		
			par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
			plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
			legend("top", legend=pop_seq, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
			dev.off()		}
	}
}







################################################
### CCRTemp like patterns
##################################################
	# We find patterns similar to differences in means for ccrTemp, since trends are similar.
		ccr_trends_diff
	# means are different:
		ccr_means_diff
	#working on:
ccr_means_diff_vec <- ccr_means_diff$.group
	names(ccr_means_diff_vec) <- as.character(ccr_means_diff$Population)
	names(ccr_means_diff_vec) <- substr(names(ccr_means_diff_vec),1,1)
print(ccr_means_diff)


#results_list <- lapply(1:NROW(sig_dif_means_meth),function(i) fit_ccr_improved_filter_methylation(sig_dif_means_meth[[i]]))
results_list <- lapply(1:NROW(sig_dif_means_meth),function(i) fit_CCR1_improved_filter_methylation(sig_dif_means_meth[[i]]))
names(results_list) <- names(sig_dif_means_meth)
NROW(which((results_list) == T))
# [1] 274

#for( gene in names(sig_dif_means_meth)){#
#	gene_df <- sig_dif_means_meth[[ which( names(sig_dif_means_meth) == gene ) ]]
#}


# Plot cols etc.
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,6,4,2)]	# B, E/O,	S, K, 
pop_vector <- c("B","O","S","K")
popcols <- popcols[order(pop_vector)]
pop_vector <- pop_vector[order(pop_vector)]

veridity_names <- names(which((results_list) == T))
count_meth_sub <- count_meth[ which(rownames(count_meth) %in% veridity_names), ]

# Plot loop
name_vector <- sprintf("%s_%02d", paste0("CCR_methylation_patternsfit_dec22_slope_ALL"), 1:50)
counter=0
for( gene in rownames(count_meth_sub)){
	#gene=rownames(count_meth_sub)[1]
	gene_nr <- which(rownames(count_meth_sub) == gene)
	count_meth_sub_sub <- count_meth_sub[which(rownames(count_meth_sub) == gene),]
	jitter_move=seq(-.5,.5,length.out=4)
	ylims <- c( min(count_meth_sub_sub,na.rm=T), max(count_meth_sub_sub,na.rm=T))
	temp_seq <- as.numeric( unique(sapply(strsplit( names(count_meth_sub_sub), "[A-Z]"),"[[",2)) )

	if( gene_nr %in% seq(1,nrow(count_meth_sub),20)){
		counter=counter+1
		png(paste0("assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/methylation_vs_temperature_tolerances/",name_vector[counter],".png"),
							width=10,height=6.5,units="in",res=200)
		par(mfrow=c(4,5), mar=c(3.6,3.5,2.1,0.4), oma=c(0,0,2,.5), mgp=c(2,.7,0))
		pop_seq <- unique(sapply(strsplit( names(count_meth_sub_sub), "[1-9]"),"[[",1))	
		for( popsis in pop_seq ){
			pop_nr <- which(pop_seq == popsis)
			# subset to pop
			pop_sub <- count_meth_sub_sub[ grep(popsis,names(count_meth_sub_sub)) ]
			if( popsis == pop_seq[1] ){	
				par(new=F)		;		ylabis="WMLvl"		;		xlabis="Temperature acclimation"	;		mainis=rownames(pop_sub)
			}else{		
				par(new=T)		;		ylabis=NA		;		xlabis=NA	;		mainis=NA
			}
			plot(temp_seq+jitter_move[pop_nr], pop_sub,		ylab=ylabis, col=popcols[pop_nr],	xlab=xlabis, main=mainis, xaxt="n",yaxt="n",xpd=NA,	xlim=c(14.5,29.5), ylim=ylims) #xlim=c(14.5,29.5))
			if( popsis == pop_seq[1] ){
				axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)	
				axis(2, ylim=ylims, xpd=NA)	}
			abline(h=rowMeans(pop_sub,na.rm=T), col=popcols[pop_nr])
		}
		#legend()
	}else{
		for( popsis in pop_seq ){
			pop_nr <- which(pop_seq == popsis)
			# subset to pop
			pop_sub <- count_meth_sub_sub[ grep(popsis,names(count_meth_sub_sub)) ]
			if( popsis == pop_seq[1] ){	
				par(new=F)		;		ylabis="WMLvl"		;		xlabis="Temperature acclimation"	;		mainis=rownames(pop_sub)
			}else{		
				par(new=T)		;		ylabis=NA		;		xlabis=NA	;		mainis=NA
			}
			plot(temp_seq+jitter_move[pop_nr], pop_sub,		ylab=ylabis, col=popcols[pop_nr],	xlab=xlabis, main=mainis, xaxt="n",yaxt="n", xpd=NA,	xlim=c(14.5,29.5), ylim=ylims) #xlim=c(14.5,29.5))
			if( popsis == pop_seq[1] ){	
				axis(2, ylim=ylims, xpd=NA)	
				axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)	}
			abline(h=rowMeans(pop_sub,na.rm=T), col=popcols[pop_nr])
		}
		if( gene_nr %in% (seq(1,nrow(count_meth_sub),20)-1) || gene_nr==nrow(count_meth_sub) ){		
			par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
			plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
			legend("top", legend=pop_seq, lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9)
			dev.off()		}
	}
}




