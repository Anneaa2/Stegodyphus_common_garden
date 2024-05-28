#!/usr/bin/env Rscript
# Script to get methylation level of genes that are involved in plasticity in response to mean temperature change

#conda activate assay_copy ; R
# cd /home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation/plastic_genes_inEXP_methlevel ; R
# source("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation/plastic_genes_inEXP_methlevel/plastic_genes_inEXP_methlevel.r")

#import libraries
#library(Hmisc) # for negation %nin% function
library(scales) # for alpha function
library(reshape2)
library(ggplot2)
library(viridis)
myvio<-function(mylist,from=0,to=1,borders=NULL,cols="lightblue2",xlab=NA,ylab=NA,...)
{
	n<-length(mylist)
	if(length(borders)==1)borders<-rep(borders,n)
	if(length(cols)==1)cols<-rep(cols,n)
	plot(0,type="n",xlim=c(0.5,n+0.5),ylim=c(from,to),axes=F,xlab=xlab,ylab=ylab,...)
	for(i in 1:n)
	{
		a<-density(mylist[[i]],from=from,to=to,na.rm=T)
		x<-c(a$x,rev(a$x))
		y<-c(a$y,rev(0-a$y))/max(a$y)/2.5+i
		polygon(y,x,border=borders[i],col=cols[i])
		points(i,median(mylist[[i]],na.rm=T),pch=21,bg="red",cex=1.5)
	}
	#axis(1,at=1:n,labels=names(mylist))
	#axis(2)
	box()
}
`%!in%` <- Negate(`%in%`)
	

## expression
		#for expression data (also fdr - why two different cutoffs) earlier: 0.05
	padj_threshold = 0.1
	#padj_threshold = 0.05
		# import expression data:
	exp_file = "/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt"
	exp_dat <- read.table(exp_file, header=T,       row.names=1, stringsAsFactors=T)
	exp_dat1 <- exp_dat[-which(rowSums(is.na(exp_dat[ , 2:ncol(exp_dat)])) > 3),]
	exp_dat2 <- exp_dat[which(rowSums(is.na(exp_dat[ , 2:ncol(exp_dat)])) > 3),]
	exp_dat_Temp <- exp_dat[	which( exp_dat$"padj_T" < padj_threshold	)	,	]
	exp_dat_Temp_pop <- exp_dat[ 	which( exp_dat$padj_T < padj_threshold & exp_dat$padj_P < padj_threshold ), ]
	exp_dat_Pop <- exp_dat[	which( exp_dat$"padj_P" < padj_threshold	)	,	]
	exp_dat_Int <- exp_dat[	which( exp_dat$"padj_I" < padj_threshold	)	,	]
	
#final[rowSums(is.na(final[ , 5:6])) == 0, ]
#rowSums(is.na(exp_dat[ , 2:ncol(exp_dat)])) == 9

	#get expression counts
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
	counts_exp_file <- "/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/data_for_downstream/norm_counts_vst_scaled_2022.csv"
	#counts_exp_file <- "/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/data_for_downstream/norm_counts_2022.csv"
	counts_exp <- read.csv(file=counts_exp_file, header=T,stringsAsFactors=T) 
	nrow(counts_exp) == nrow(exp_dat)
	
		## OBS IF FALSE
	#remove genes with low variation in expression (the ones with no pvalues calculated). They are in exp_dat2
	counts_exp1 <- counts_exp[which( counts_exp[,1] %in% exp_dat1[,1]),]
	nrow(counts_exp1) == nrow(exp_dat1)

	
		# get only plastic genes counts (and the other categories)
	plastic_counts_names <- intersect(exp_dat_Temp[,1],counts_exp[,1])
	pop_counts_names <- intersect(exp_dat_Pop[,1],counts_exp[,1])
	plastic_pop_counts_names <- intersect(exp_dat_Temp_pop[,1],counts_exp[,1])
	interaction_counts_names <- intersect(exp_dat_Int[,1],counts_exp[,1])
	#plastic_counts <- counts_exp[counts_exp[,1] %in% plastic_counts_names,]
	
	non_plastic_counts_names <- setdiff(counts_exp1[,1], exp_dat_Temp[,1])
	non_pop_counts_names <- setdiff(counts_exp1[,1], exp_dat_Pop[,1])
	non_plastic_pop_counts_names <- setdiff(counts_exp1[,1], exp_dat_Temp_pop[,1])
	non_interaction_counts_names <- setdiff(counts_exp1[,1], exp_dat_Int[,1])
	#non_plastic_counts <- counts_exp[counts_exp[,1] %in% non_plastic_counts_names,]
	
	

	# meta data for expression: 
	meta_exp <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/analyzing/RNAextradction_samples_results022618_dosimp.txt",sep="\t", row.names=1,header=T)
	#meta_exp <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/analyzing/RNAextradction_samples_results022618_dosimp.txt",sep="\t", row.names=1,header=T)
	
	meta_exp$pop_temp <- gsub("N|E","O",meta_exp$pop_temp)
	meta_exp$population <- gsub("N|E","O",meta_exp$population)

		# order expression files to be pop/temp group means.		
	rownames(counts_exp1) <- counts_exp1[,1]	; counts_exp1 <- counts_exp1[,-1]
	counts_mean <- aggregate(t(counts_exp1[1:nrow(counts_exp1),]), by=list(meta_exp$pop_temp), mean)
	rownames(counts_mean) <- counts_mean[,1]		; counts_mean <- counts_mean[,-1]	; counts_mean <- t(counts_mean)
	#plastic_counts_mean <- data.frame(t(plastic_counts_mean))
	



# methylation level (weighted methylation level, genewise per population)
	count_meth <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/wmlvl_allpops_d10D32_gene_cpg.txt",sep="\t", row.names=1,header=T)
	colnames(count_meth) <- gsub("N|E","O",colnames(count_meth))
	
	# common genes between methylation and expression (since methylation data has fewer genes - possibly because of 0 genes?)
	common_genes1 <- merge(count_meth,exp_dat, by.x=0,by.y=1)
	
	common_genes1$Plasticity <- TRUE	; 	common_genes1$Pop_Effect <- TRUE	;	common_genes1$Plasticity <- TRUE 	; common_genes1$Interaction <- TRUE
	common_genes1[ which( common_genes1$Row.names %in% plastic_counts_names), "Plasticity" ] <- "Plastic"
	common_genes1[ which( common_genes1$Row.names %in% non_plastic_counts_names), "Plasticity" ] <- "Non-Plastic"
	
	common_genes1[ which( common_genes1$Row.names %in% pop_counts_names), "Pop_Effect" ] <- "Population"
	common_genes1[ which( common_genes1$Row.names %in% non_pop_counts_names), "Pop_Effect" ] <- "NAN"
	
	common_genes1[ which( common_genes1$Row.names %in% plastic_pop_counts_names), "Plastic_Pop_effect" ] <- "Plastic_Pop"
	common_genes1[ which( common_genes1$Row.names %in% non_plastic_pop_counts_names), "Plastic_Pop_effect" ] <- "NAN"
	
	common_genes1[ which( common_genes1$Row.names %in% interaction_counts_names), "Interaction" ] <- "Interaction"
	common_genes1[ which( common_genes1$Row.names %in% non_interaction_counts_names), "Interaction" ] <- "NAN"
	
	common_genes <- common_genes1[, grep("^[A-Z]", colnames(common_genes1))]
	common_genes_long <- melt(common_genes, id.vars=c("Row.names", "Plasticity","Pop_Effect","Plastic_Pop_effect","Interaction"))
	
	
meth_list_plastic <- list(  common_genes_long[ common_genes_long[,"Plasticity"]=="Plastic", "value" ],
							 common_genes_long[ common_genes_long[,"Plasticity"]=="Non-Plastic", "value" ])
names(meth_list_plastic) <- c("Plastic","Non-plastic")
names(meth_list_plastic) <- paste0(names(meth_list_plastic),"ally\nexpressed genes","\n(",sapply(meth_list_plastic,NROW)/20," genes)")
color_pattern <- c(	viridis(1, begin = 0, end = .03), "white", 		
					viridis(4, begin = .23, end = .5), "white",		
					viridis(1, begin = .8, end = .85), "white",		
					viridis(1, begin = 1, end = 1))
meth_color <- color_pattern[8]
bord_col <- color_pattern[10]


#pdf(paste0("plastic_vs_nonPlastic_genes_methylation_VIOPLOT_FDR_",padj_threshold,".pdf"))
#myvio(meth_list_plastic, cols=meth_color, from=min(common_genes_long$value,na.rm=T), to=max(common_genes_long$value,na.rm=T),xlab="Gene status in gene expression in response to mean temperature change", ylab="Weighted Methylation Level")
#axis(1,at=1:length(meth_list_plastic),lwd.ticks=0,labels=paste0(names(meth_list_plastic),"\n(N: ",sapply(meth_list_plastic,NROW)/20," g per group/treat)") )
#axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
##text(0.375,-22,"\\\\",xpd=NA,cex=1.2,srt=90)
#dev.off()


vioplot_df1 <- data.frame(meth_list_plastic[1],names(meth_list_plastic[1]))		; names(vioplot_df1) <- c("WML", "Plasticity")
shapiro.test(sample(vioplot_df1[,1],size=5000))	# significantly different from normal
vioplot_df2 <- data.frame(meth_list_plastic[2],names(meth_list_plastic[2]))		; names(vioplot_df2) <- c("WML", "Plasticity")
shapiro.test(sample(vioplot_df2[,1],size=5000))	# significantly different from normal
median_vec <- c(median(vioplot_df1[,1], na.rm=T),	median(vioplot_df2[,1], na.rm=T))
mean_vec <- c(mean(vioplot_df1[,1], na.rm=T),	mean(vioplot_df2[,1], na.rm=T))
vioplot_df <- rbind(vioplot_df1, vioplot_df2)
wilcox.test(vioplot_df1[,1],vioplot_df2[,1])	# Man whitney U/wilcoxon rank sum test: Null hyp: These data are from same population (ie. are they similar?) (non parametric test)
	# They are statistically different (***).
	
library(coin)
# wilcox_test(vioplot_df1[,1],vioplot_df2[,1]) # cannot get this to work


plotti <- ggplot( aes(x=Plasticity, y=WML), data=vioplot_df ) + 
    theme_bw() 
plotti <- plotti +
	geom_violin(width=1.1, position = position_dodge(), show.legend=F, fill=meth_color, colour="grey48", lwd=.01) +	# add violins # , draw_quantiles=c(.25,.5,.75)
	annotate("point", x=c(1.01,2.01), y=median_vec, colour="black", size=3.5) +
	annotate("point", x=c(1.000,2.000), y=mean_vec, colour="grey50", size=3.5,pch="-") +
	theme(axis.text.x=element_text(family="Helvetica", size=9, colour="black")) + #c("black","transparent","black","transparent","black","black","transparent","black"))) +
	theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
	theme(plot.margin=unit(c(0,.3,0,.4), "lines"), axis.title.y=element_text(family="Helvetica", size=9, colour="black")) +
	ylab("Weighted methylation level") +
	theme(axis.ticks.x=element_blank(), panel.border=element_rect(colour="grey60"))
pdf(paste0("assay_new_setup/derived_data_and_code/methylation/plastic_genes_inEXP_methlevel/S24_plastic_vs_nonPlastic_genes_methylation_VIOPLOT_FDR_",padj_threshold,".pdf"), height=4, width=3.5)
plotti
dev.off()

# for png
plotti <- ggplot( aes(x=Plasticity, y=WML), data=vioplot_df ) + 
    theme_bw() 
plotti <- plotti +
	geom_violin(width=1.1, position = position_dodge(), show.legend=F, fill=meth_color, colour="grey48", lwd=.3) +	# add violins # , draw_quantiles=c(.25,.5,.75)
	annotate("point", x=c(1.0,2.0), y=median_vec, colour="black", size=1.5) +
	annotate("point", x=c(1.000,2.000), y=mean_vec, colour="grey50", size=5, pch="-") +
	theme(axis.text.x=element_text(family="Helvetica", size=9, colour="black")) + #c("black","transparent","black","transparent","black","black","transparent","black"))) +
	theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
	theme(plot.margin=unit(c(0,.3,0,.4), "lines"), axis.title.y=element_text(family="Helvetica", size=9, colour="black")) +
	ylab("Weighted methylation level") +
	theme(axis.ticks.x=element_blank(), panel.border=element_rect(colour="grey60"))
ggsave(paste0("assay_new_setup/derived_data_and_code/methylation/plastic_genes_inEXP_methlevel/S24_plastic_vs_nonPlastic_genes_methylation_VIOPLOT_FDR_",padj_threshold,".png"), plot=plotti, height=4, width=3.5, units="in",dpi=800)







	# cut up plastic genes based on methylation level and get their respective expression count (normalized and vst scaled)
meth_list_plastic_cut <- cut(meth_list_plastic$Plastic,c(-1,0.2,0.7,1),labels=c("[0,0.2]","(0.2,0.7]","(0.7,1]"))
mylist<-list( counts_mean[meth_list_plastic_cut=="[0,0.2]"],counts_mean[meth_list_plastic_cut=="(0.2,0.7]"],counts_mean[meth_list_plastic_cut=="(0.7,1]"])
names(mylist)<-c("Low","Medium","High")

pdf(paste0("plastic_genes_mehylationLvL_vs_EXP_VIOPLOT_FDR_",padj_threshold,".pdf"), width=7, height=5)
myvio(mylist,from=min(counts_mean,na.rm=T),to=max(counts_mean,na.rm=T),xlab="Weighted methylation level of Plastic genes",ylab="Normalized Expression Count")
axis(1,at=1:length(mylist),lwd=0,labels=paste0(names(mylist),"\nN: ",sapply(mylist,NROW)/20) )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
dev.off()
	# conclusion, level of expression not clearly different between highly and lowly methylated genes showing plasticity in gene expression

 # plot level of plasticity to check for signs.
 # common_genes1$stat_T -> statistics of the genes in relation to plasticity
mylist<-list( common_genes1$stat_T[meth_list_plastic_cut=="[0,0.2]"],common_genes1$stat_T[meth_list_plastic_cut=="(0.2,0.7]"],common_genes1$stat_T[meth_list_plastic_cut=="(0.7,1]"])
names(mylist)<-c("Low","Medium","High")

pdf(paste0("plastic_genes_mehylationLvL_vs_lvlOfPlasticity_stat_VIOPLOT_FDR_",padj_threshold,".pdf"), width=7, height=5)
myvio(mylist,from=min(common_genes1$stat_T,na.rm=T),to=max(common_genes1$stat_T,na.rm=T),xlab="Weighted methylation level of Plastic genes",ylab="Test-statistics (LRT-test)")
axis(1,at=1:length(mylist),lwd=0,labels=paste0(names(mylist),"\nN: ",sapply(mylist,NROW)/20) )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
dev.off()

mylist<-list( common_genes1$padj_T[meth_list_plastic_cut=="[0,0.2]"],common_genes1$padj_T[meth_list_plastic_cut=="(0.2,0.7]"],common_genes1$padj_T[meth_list_plastic_cut=="(0.7,1]"])
names(mylist)<-c("Low","Medium","High")

pdf(paste0("plastic_genes_mehylationLvL_vs_lvlOfPlasticity_fdr_VIOPLOT_FDR_",padj_threshold,".pdf"), width=7, height=5)
myvio(mylist,from=min(common_genes1$padj_T,na.rm=T),to=max(common_genes1$padj_T,na.rm=T),xlab="Weighted methylation level of Plastic genes",ylab="FDR (LRT-test)")
axis(1,at=1:length(mylist),lwd=0,labels=paste0(names(mylist),"\nN: ",sapply(mylist,NROW)/20) )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
dev.off()

#plot combi plot


pdf(paste0("plasticity_at_play_panel_FDR_",padj_threshold,".pdf"), width=7, height=5)
par(mfrow=c(1,2),cex=.8,cex.lab=.9, cex.axis=.8)
	# subplot 1
myvio(meth_list_plastic, from=min(common_genes_long$value,na.rm=T), to=max(common_genes_long$value,na.rm=T),xlab="Gene status in gene expression", ylab="Weighted Methylation Level")
axis(1,at=1:length(meth_list_plastic),lwd.ticks=0,labels=paste0(names(meth_list_plastic),"\n",sapply(meth_list_plastic,NROW)/20," genes/group") )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))

	#subplot 2
mylist<-list( common_genes1$pvalue_T[meth_list_plastic_cut=="[0,0.2]"],common_genes1$pvalue_T[meth_list_plastic_cut=="(0.2,0.7]"],common_genes1$pvalue_T[meth_list_plastic_cut=="(0.7,1]"])
names(mylist)<-c("Low","Medium","High")
myvio(mylist,from=min(common_genes1$padj_T,na.rm=T),to=max(common_genes1$padj_T,na.rm=T),xlab="Weighted methylation level of Plastic genes",ylab="Pvalue (LRT-test)")
axis(1,at=1:length(mylist),lwd=0,labels=paste0(names(mylist),"\nN: ",sapply(mylist,NROW)/20) )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
dev.off()


# remove lowly expressed genes or zero genes?


meth_list_all <- list(  common_genes_long[ common_genes_long[,"Plasticity"]=="Plastic", "value" ],
						common_genes_long[ common_genes_long[,"Plasticity"]=="Non-Plastic", "value" ],
						common_genes_long[ common_genes_long[,"Pop_Effect"]=="Population", "value" ],
						common_genes_long[ common_genes_long[,"Pop_Effect"]=="NAN", "value" ],
						common_genes_long[ common_genes_long[,"Plastic_Pop_effect"]=="Plastic_Pop", "value" ],
						common_genes_long[ common_genes_long[,"Interaction"]=="Interaction", "value" ])
names(meth_list_all) <- c("Plastic","Non-Plastic","Pop_effect","Non-Pop","Plastic & Pop","Interaction")

pdf(paste0("plastic_vs_allCategories_genes_methylation_VIOPLOT_FDR_",padj_threshold,".pdf"), width=7, height=5)
par(cex=.8,cex.lab=.9, cex.axis=.8)
myvio(meth_list_all, from=min(common_genes_long$value,na.rm=T), to=max(common_genes_long$value,na.rm=T),xlab="Gene status in gene expression in response to mean temperature change", ylab="Weighted Methylation Level")
axis(1,at=1:length(meth_list_all),lwd=0,labels=paste0(names(meth_list_all),"\nN: ",sapply(meth_list_all,NROW)/20) )
axis(2)#,at=c(-25,-10,0,10),labels=c("NA","-10","0","10"))
#text(0.375,-22,"\\\\",xpd=NA,cex=1.2,srt=90)
dev.off()

