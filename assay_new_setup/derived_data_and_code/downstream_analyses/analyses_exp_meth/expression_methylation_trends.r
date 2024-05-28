

 ## For SOV results and temperature
 ## /faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/only_temp_prec/genewise_SMPS_res/
	## temp_prec_partial_mantel_res_genes_smps_refpop_ALL.csv

# cd /faststorage/project/spider2/assay_study_anne/analyses_exp_meth/ ; 
# conda activate assay_copy; R
# setwd("/faststorage/project/spider2/assay_study_anne/analyses_exp_meth/") 
# Rscript
library(scales)
library(vioplot)
library(vegan)
library(reshape2)
library(agricolae)
outplace="assay_new_setup/derived_data_and_code/downstream_analyses/analyses_exp_meth/"


#popcols=c("steelblue1","blue","tan4","green3","snow4","red") # B, K, N, S, G, E/O
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")


# methylation level (weighted methylation level, genewise per population)
	count_meth <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/wmlvl_allpops_d10D32_gene_cpg.txt",sep="\t", row.names=1,header=T)
	# /home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation
	colnames(count_meth) <- gsub("N|E","O",colnames(count_meth))
# pick out which genes are significant with respect to population (since none were with respect to treatment temperature.
	meth_dat <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/methylation/smps/gene_region/DSS_output/DSS_result_addit_d10D32_gene_CPG_population_significant.txt",sep="\t", row.names=1,header=T)


# expression genes with effect of population: filter from this file:
	exp_dat <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt",sep="\t", row.names=1,header=T)
# file with counts of expreesion: (that can be compared between samples, but not wihtin)
	norm_counts_exp <- read.csv(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/data_for_downstream/norm_counts_2022.csv",row.names=1)
	vst_counts_exp <- read.csv(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/results/data_for_downstream/norm_counts_vst_scaled_2022.csv",row.names=1)
# meta data for expression: 
	meta_exp <- read.table(file="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/expression/analyzing/RNAextradction_samples_results022618_dosimp.txt",sep="\t", row.names=1,header=T)
	meta_exp$pop_temp <- gsub("N|E","O",meta_exp$pop_temp)
	meta_exp$population <- gsub("N|E","O",meta_exp$population)
	
# get genes significan for pop effect
	fdr_thresh <- 0.1
	fdr_thresh <- 0.05
	exp_dat_sub <- exp_dat[which(exp_dat$padj_P < fdr_thresh),]
	exp_dat_sub_temp <- exp_dat[which(exp_dat$padj_T < fdr_thresh),]
	meth_dat_sub <- meth_dat[which(meth_dat$fdrs < fdr_thresh),]

	common_genes <- merge(meth_dat_sub,exp_dat_sub, by.x=0,by.y=1)
	common_genes_all <- merge(count_meth,exp_dat, by.x=0,by.y=1)

nrow(exp_dat_sub)
nrow(meth_dat_sub)
nrow(common_genes)
nrow(common_genes_all)
# going with this:
' fdr= .05 for expression and methylation
[1] 12089
[1] 1136
[1] 572
[1] 34971
'
# fdr =.1 in exp as well 
'
nrow(common_genes)
[1] 13998
[1] 1537
[1] 844
[1] 34971
'

# overlap with interactions
#exp_dat_sub_int <- exp_dat[which(exp_dat$padj_I < .1),]
exp_dat_sub_int <- exp_dat[which(exp_dat$padj_I < fdr_thresh),]
common_genes_int <- merge(meth_dat_sub,exp_dat_sub_int, by.x=0,by.y=1)
nrow(exp_dat_sub_int)
nrow(common_genes_int)
' threshold .1
[1] 909
[1] 41
'' threshold .05
> nrow(exp_dat_sub_int)
[1] 711
> nrow(common_genes_int)
[1] 25
'

# Subset counts files to common genes:
	count_meth_sub <- merge(count_meth, common_genes, by.x=0, by.y=1, all.y=T)[,1:21]
	rownames(count_meth_sub) <- count_meth_sub[,1]	; count_meth_sub <- count_meth_sub[,-1]
	norm_counts_exp_sub <- merge(norm_counts_exp, common_genes, by.x=0, by.y=1, all.y=T)[,1:200]
	vst_counts_exp_sub <- merge(vst_counts_exp, common_genes, by.x=0, by.y=1, all.y=T)[,1:200]
	# similar for all genes
	count_meth_all_common <- merge(count_meth, common_genes_all, by.x=0, by.y=1, all.y=T)[,1:21]
	rownames(count_meth_all_common) <- count_meth_all_common[,1]	; count_meth_all_common <- count_meth_all_common[,-1]
	norm_counts_exp_all_common <- merge(norm_counts_exp, common_genes_all, by.x=0, by.y=1, all.y=T)[,1:200]
	vst_counts_exp_all_common <- merge(vst_counts_exp, common_genes_all, by.x=0, by.y=1, all.y=T)[,1:200]
	
# fix expression files to calculate means for each pop/temp group
	rownames(norm_counts_exp_sub) <- norm_counts_exp_sub[,1]	; norm_counts_exp_sub <- norm_counts_exp_sub[,-1]
	norm_counts_exp_sub_mean <- aggregate(t(norm_counts_exp_sub[1:nrow(norm_counts_exp_sub),]), by=list(meta_exp$pop_temp), mean)
	rownames(norm_counts_exp_sub_mean) <- norm_counts_exp_sub_mean[,1]		; norm_counts_exp_sub_mean <- norm_counts_exp_sub_mean[,-1]
	norm_counts_exp_sub_mean <- data.frame(t(norm_counts_exp_sub_mean))

	rownames(vst_counts_exp_sub) <- vst_counts_exp_sub[,1]	; vst_counts_exp_sub <- vst_counts_exp_sub[,-1]
	vst_counts_exp_sub_mean <- aggregate(t(vst_counts_exp_sub[1:nrow(vst_counts_exp_sub),]), by=list(meta_exp$pop_temp), mean)
	rownames(vst_counts_exp_sub_mean) <- vst_counts_exp_sub_mean[,1]		; vst_counts_exp_sub_mean <- vst_counts_exp_sub_mean[,-1]
	vst_counts_exp_sub_mean <- data.frame(t(vst_counts_exp_sub_mean))
	# for all genes:
	rownames(norm_counts_exp_all_common) <- norm_counts_exp_all_common[,1]	; norm_counts_exp_all_common <- norm_counts_exp_all_common[,-1]
	norm_counts_exp_all_common_mean <- aggregate(t(norm_counts_exp_all_common[1:nrow(norm_counts_exp_all_common),]), by=list(meta_exp$pop_temp), mean)
	rownames(norm_counts_exp_all_common_mean) <- norm_counts_exp_all_common_mean[,1]		; norm_counts_exp_all_common_mean <- norm_counts_exp_all_common_mean[,-1]
	norm_counts_exp_all_common_mean <- data.frame(t(norm_counts_exp_all_common_mean))

	rownames(vst_counts_exp_all_common) <- vst_counts_exp_all_common[,1]	; vst_counts_exp_all_common <- vst_counts_exp_all_common[,-1]
	vst_counts_exp_all_common_mean <- aggregate(t(vst_counts_exp_all_common[1:nrow(vst_counts_exp_all_common),]), by=list(meta_exp$pop_temp), mean)
	rownames(vst_counts_exp_all_common_mean) <- vst_counts_exp_all_common_mean[,1]		; vst_counts_exp_all_common_mean <- vst_counts_exp_all_common_mean[,-1]
	vst_counts_exp_all_common_mean <- data.frame(t(vst_counts_exp_all_common_mean))
# order methylation and expression data:
	norm_counts_exp_sub_mean <- norm_counts_exp_sub_mean[order(rownames(norm_counts_exp_sub_mean)),]
	vst_counts_exp_sub_mean <- vst_counts_exp_sub_mean[order(rownames(vst_counts_exp_sub_mean)),]
	count_meth_sub <- count_meth_sub[order(rownames(count_meth_sub)),]
	# all genes
	norm_counts_exp_all_common_mean <- norm_counts_exp_all_common_mean[order(rownames(norm_counts_exp_all_common_mean)),]
	vst_counts_exp_all_common_mean <- vst_counts_exp_all_common_mean[order(rownames(vst_counts_exp_all_common_mean)),]
	count_meth_all_common <- count_meth_all_common[order(rownames(count_meth_all_common)),]


#######################################################################################################
###		Does methylation group according to populations acroos genes with popeffect?	 ##############
#######################################################################################################

output_type = "png"

# Subset counts files to genes with pop effect:
	count_meth_fPCA <- merge(count_meth, meth_dat_sub, by.x=0, by.y=0, all.y=T)[,1:21]
	rownames(count_meth_fPCA) <- count_meth_fPCA[,1]	; count_meth_fPCA <- count_meth_fPCA[,-1]

	#norm_counts_exp_sub <- merge(norm_counts_exp, common_genes, by.x=0, by.y=1, all.y=T)[,1:200]
	#vst_counts_exp_sub <- merge(vst_counts_exp, common_genes, by.x=0, by.y=1, all.y=T)[,1:200]


pca_popMeth_genes <- prcomp(t(na.omit(count_meth_fPCA)))
summa <- summary(pca_popMeth_genes)
summa
'   Importance of components: (0.05)
                          PC1    PC2    PC3     PC4     PC5     PC6     PC7
Standard deviation     2.0161 1.2083 0.9188 0.36004 0.31754 0.31057 0.30048
Proportion of Variance 0.5464 0.1963 0.1135 0.01743 0.01356 0.01297 0.01214
Cumulative Proportion  0.5464 0.7427 0.8562 0.87363 0.88719 0.90016 0.91230
                           PC8     PC9    PC10    PC11    PC12    PC13    PC14
Standard deviation     0.28464 0.27071 0.25738 0.24636 0.24217 0.23639 0.22810
Proportion of Variance 0.01089 0.00985 0.00891 0.00816 0.00788 0.00751 0.00699
Cumulative Proportion  0.92319 0.93304 0.94195 0.95011 0.95799 0.96550 0.97250
                          PC15    PC16    PC17    PC18   PC19      PC20
Standard deviation     0.21745 0.20701 0.19884 0.19416 0.1928 9.448e-16
Proportion of Variance 0.00636 0.00576 0.00532 0.00507 0.0050 0.000e+00
Cumulative Proportion  0.97886 0.98462 0.98993 0.99500 1.0000 1.000e+00

'

pca_pop_order <- sapply(strsplit(rownames(pca_popMeth_genes$x),"1|2"),"[[",1)
pop_order <- vector()
popcols <- popcols[grep("N",names(popcols), invert=T)]
for( each in pca_pop_order ){	pop_order <- append(pop_order,grep(each,names(popcols)))	}
popcols1 <- popcols[pop_order]

#if (	output_type == "pdf"	)	{
#	pdf(paste0(outplace,"pca_methy_plots.pdf"),width=7,height=7) 		}else if ( 	output_type == "png"	){
#	png(paste0(outplace,"pca_methy_plots.png"),width=7,height=7,units="in",res=400) 	}

if (	output_type == "pdf"	)	{
	pdf(paste0(outplace,"pca_methy_plots_",fdr_thresh, "_mod.pdf"),width=7,height=7) 		
}else if ( 	output_type == "png"	){
	png(paste0(outplace,"pca_methy_plots_",fdr_thresh, "_mod.png"),width=7,height=7,units="in",res=400) 	}

par(mar=c(4.7,4.2,3.7,1.5), mfrow=c(2,2),cex=.7)
plottit <- ordiplot(pca_popMeth_genes,choices=c(1,2),type="n", 
			xlim=c(-3,3),ylim=c(-3,3), cex=0.8, main="Methylation, population",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"))
points(plottit, "sites",col=popcols1, pch=1)
ordiellipse(plottit,popcols1,conf=0.80, col=levels(as.factor(popcols1)))
legend("topright",legend=sort(names(popcols)),col=popcols[sort(names(popcols))],pch=1, 
			inset=c(-(.42),-.21), xpd=NA, ncol=2)

plottit <- ordiplot(pca_popMeth_genes,choices=c(2,3),type="n", 
			xlim=c(-3,3),ylim=c(-3,3), cex=0.8, main="Methylation, population",
			xlab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"))
points(plottit, "sites",col=popcols1, pch=1)#
ordiellipse(plottit,popcols1,conf=0.80, col=levels(as.factor(popcols1)))

screeplot(pca_popMeth_genes, type="l", main="Screeplot Methylation")
biplot(pca_popMeth_genes, main="Biplot methylation")
dev.off()



#######################################################################################################
###		Does EXPRESSION group according to populations acroos genes with popeffect?	 	##############
#######################################################################################################


#subset counts file to genes with pop effect
	vst_counts_pca <- merge(vst_counts_exp, exp_dat_sub, by.x=0, by.y=1, all.y=T)[,1:200]
	rownames(vst_counts_pca) <- vst_counts_pca[,1]	; vst_counts_pca <- vst_counts_pca[,-1]

	
pca_popExp_genes <- prcomp(t(na.omit(vst_counts_pca)))
summa <- summary(pca_popExp_genes)
summa
'Importance of components:
                           PC1     PC2      PC3      PC4      PC5     PC6
Standard deviation     36.8752 34.7727 27.50402 21.77716 20.10785 18.6252
Proportion of Variance  0.1639  0.1457  0.09116  0.05715  0.04873  0.0418
Cumulative Proportion   0.1639  0.3096  0.40074  0.45789  0.50662  0.5484
                           PC7     PC8      PC9     PC10     PC11     PC12
Standard deviation     17.2369 16.7965 13.43278 11.92988 11.29166 11.10551
Proportion of Variance  0.0358  0.0340  0.02174  0.01715  0.01537  0.01486
Cumulative Proportion   0.5842  0.6182  0.63997  0.65712  0.67249  0.68735
                          PC13    PC14    PC15    PC16    PC17    PC18   PC19
Standard deviation     9.91977 9.08043 8.93805 8.32391 8.14298 7.70135 7.4571
Proportion of Variance 0.01186 0.00994 0.00963 0.00835 0.00799 0.00715 0.0067
Cumulative Proportion  0.69921 0.70914 0.71877 0.72712 0.73511 0.74226 0.7490

'


meta_exp_dat <- data.frame(rownames(meta_exp),meta_exp$population) 	; colnames(meta_exp_dat) <- c("rna_nr","population")
	# change rna nr to be similar to the other df (RNA_0001) not RNA 1
meta_exp_dat$rna_nr_new <- sapply(strsplit(as.character(meta_exp_dat$rna_nr)," "), function(x)
    paste(x[1], sprintf("%04d",as.numeric(x[2])), sep="_")
)

pop_order <- vector()
for( each in meta_exp_dat$population ){	pop_order <- append(pop_order,grep(each,names(popcols)))	}
popcols2 <- popcols[pop_order]
#output_type = "png"



#if (	output_type == "pdf"	)	{ 	
#	pdf(paste0(outplace,"pca_EXPpop_plots.pdf"),width=7,height=7) 		}else if ( 	output_type == "png"	){	
#	png(paste0(outplace,"pca_EXPpop_plots.png"),width=7,height=7,units="in",res=400) 	}
#par(mar=c(4.7,4.2,3.7,1.5), mfrow=c(2,2), cex=.7)
#plottit <- ordiplot(pca_popExp_genes,choices=c(1,2),type="n", 
#			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression, population",
#			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"))
#points(plottit, "sites",col=popcols2, pch=1)#
#ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))
#legend("topright",legend=sort(names(popcols)),col=popcols[sort(names(popcols))],pch=1, 
#			inset=c(-(.42),-.21), xpd=NA, ncol=2)

#plottit <- ordiplot(pca_popExp_genes,choices=c(2,3),type="n", 
#			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression, population",
#			xlab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"))
#points(plottit, "sites",col=popcols2, pch=1)#
#ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))

#screeplot(pca_popExp_genes, type="l", main="Screeplot Expression, population")
#biplot(pca_popExp_genes, main="Biplot Expression")
#dev.off()


if (	output_type == "pdf"	)	{ 	
	pdf(paste0(outplace,"pca_EXPpop_plots_",fdr_thresh, "_mod.pdf"),width=3.5,height=7) 		}else if ( 	output_type == "png"	){	
	png(paste0(outplace,"pca_EXPpop_plots_",fdr_thresh, "_mod.png"),width=3.5,height=7,units="in",res=400) 	}

par(mar=c(4.7,4.2,3.7,1.5), mfrow=c(2,1), cex=.7, oma=c(3,0,0,0))
plottit <- ordiplot(pca_popExp_genes,choices=c(1,2),type="n", 
			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), xaxt="n", yaxt="n")
axis(1, xlim=c(-100,100))
axis(2, ylim=c(-100,100), las=2)
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))

if( length( levels(as.factor(popcols2)) ) == 5) { 
	mtext("b", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
if( length( levels(as.factor(popcols2)) ) == 4) { 
	mtext("a", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
plottit <- ordiplot(pca_popExp_genes,choices=c(1,3),type="n", 
			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"), xaxt="n", yaxt="n")
axis(1, xlim=c(-100,100))
axis(2, ylim=c(-100,100), las=2)
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))
legend("bottom",legend=sort(names(popcols)),col=popcols[sort(names(popcols))],pch=1, 
			inset=c(0,-.48), xpd=NA, ncol=2, horiz=T)
#screeplot(pca_popExp_genes, type="l", main="Screeplot Expression, population")
#biplot(pca_popExp_genes, main="Biplot Expression")
if( length( levels(as.factor(popcols2)) ) == 5) { 
		mtext("d", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
if( length( levels(as.factor(popcols2)) ) == 4) { 
	mtext("c", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
dev.off()



#######################################################################################################
###		How Does EXPRESSION group according to temperature acroos genes with temperature effect?	 	##############
#######################################################################################################


#subset counts file to genes with pop effect
	vst_counts_pca_temp <- merge(vst_counts_exp, exp_dat_sub_temp, by.x=0, by.y=1, all.y=T)[,1:200]
	rownames(vst_counts_pca_temp) <- vst_counts_pca_temp[,1]	; vst_counts_pca_temp <- vst_counts_pca_temp[,-1]

	
pca_tempExp_genes <- prcomp(t(na.omit(vst_counts_pca_temp)))
summa <- summary(pca_tempExp_genes)
summa
'Importance of components: (0.05)
                           PC1     PC2      PC3      PC4      PC5      PC6
Standard deviation     38.9251 31.2840 25.36480 20.54854 20.11738 17.19972
Proportion of Variance  0.2086  0.1348  0.08858  0.05813  0.05572  0.04073
Cumulative Proportion   0.2086  0.3433  0.43193  0.49007  0.54579  0.58652
                           PC7      PC8      PC9     PC10     PC11    PC12
Standard deviation     16.2153 14.70300 12.58887 10.70966 10.48077 9.96841
Proportion of Variance  0.0362  0.02976  0.02182  0.01579  0.01512 0.01368
Cumulative Proportion   0.6227  0.65248  0.67430  0.69009  0.70522 0.71890
                          PC13    PC14    PC15    PC16    PC17    PC18    PC19
Standard deviation     9.41719 8.73020 8.35598 7.74413 7.39211 7.07399 6.78978
Proportion of Variance 0.01221 0.01049 0.00961 0.00826 0.00752 0.00689 0.00635
Cumulative Proportion  0.73111 0.74160 0.75122 0.75947 0.76700 0.77389 0.78023

'
library(viridis)
tempcols <- plasma(5) ; names(tempcols) <- c("15","19","23","25","29")


meta_exp_dat_temp <- data.frame(rownames(meta_exp),meta_exp$Temperature) 	; colnames(meta_exp_dat_temp) <- c("rna_nr","temperature")
	# change rna nr to be similar to the other df (RNA_0001) not RNA 1
meta_exp_dat_temp$rna_nr_new <- sapply(strsplit(as.character(meta_exp_dat$rna_nr)," "), function(x)
    paste(x[1], sprintf("%04d",as.numeric(x[2])), sep="_")
)


pop_order <- vector()
for( each in meta_exp_dat_temp$temperature ){	pop_order <- append(pop_order,grep(each,names(tempcols)))	}
tempcols2 <- tempcols[pop_order]
#output_type = "png"

if (	output_type == "pdf"	)	{ 	
	pdf(paste0(outplace,"pca_EXPtemp_plots_",fdr_thresh, "_mod.pdf"),width=3.5,height=7) 		}else if ( 	output_type == "png"	){	
	png(paste0(outplace,"pca_EXPtemp_plots_",fdr_thresh, "_mod.png"),width=3.5,height=7,units="in",res=400) 	}

par(mar=c(4.7,4.2,3.7,1.5), mfrow=c(2,1), cex=.7, oma=c(3,0,0,0))
plottit <- ordiplot(pca_tempExp_genes,choices=c(1,2),type="n", 
			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), xaxt="n", yaxt="n")
axis(1, xlim=c(-100,100))
axis(2, ylim=c(-100,100), las=2)
points(plottit, "sites",col=tempcols2, pch=1)#
ordiellipse(plottit,tempcols2,conf=0.80, col=levels(as.factor(tempcols2)))

if( length( levels(as.factor(tempcols2)) ) == 5) { 
	mtext("b", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
if( length( levels(as.factor(tempcols2)) ) == 4) { 
	mtext("a", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}


plottit <- ordiplot(pca_tempExp_genes,choices=c(1,3),type="n", 
			xlim=c(-100,100),ylim=c(-100,100), cex=0.8, main="Expression",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"), xaxt="n", yaxt="n")
axis(1, xlim=c(-100,100))
axis(2, ylim=c(-100,100), las=2)
points(plottit, "sites",col=tempcols2, pch=1)#
ordiellipse(plottit,tempcols2,conf=0.80, col=levels(as.factor(tempcols2)))
legend("bottom",legend=sort(unique(names(tempcols2))),col=tempcols2[sort(unique(names(tempcols2)))],pch=1, 
			inset=c(0,-.48), xpd=NA, ncol=2, horiz=T)
#screeplot(pca_popExp_genes, type="l", main="Screeplot Expression, population")
#biplot(pca_popExp_genes, main="Biplot Expression")
if( length( levels(as.factor(tempcols2)) ) == 5) { 
		mtext("d", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
if( length( levels(as.factor(tempcols2)) ) == 4) { 
	mtext("c", outer=F, xpd = T, side = 2, las = 2, line = 1, padj = 0, at = par()$usr[4])
}
dev.off()





###########################################################################
###		Does higher methylation result in higher expression? ##############
###########################################################################

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

# correlate all rows, to see if higher methylation gives higher expression:
	vec_cor_norExp_meth <- sapply(1:nrow(count_meth_sub), function(i) cor(t(count_meth_sub[i,]), t(vst_counts_exp_sub_mean[i,])))
	Non_vst_vec_cor_norExp_meth <- sapply(1:nrow(count_meth_sub), function(i) cor(t(count_meth_sub[i,]), t(norm_counts_exp_sub_mean[i,])))
		# when plottet also looks normal ish (or sligtly bimodal)
# histogram of correlations, check if left skewed. # very subtle left skew.
'	pdf(file=paste0(outplace,"hist_corre_normExp_Meth_popRespGenes_",fdr_thresh, ".pdf"), width=5, height=4)
	par(cex=.7)
	hist(vec_cor_norExp_meth, main="Weighted methylation level\nvs. vst transformed normalized expression level (log2foldchange)\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=15,xlim=c(-1,1))
	abline(v=mean(vec_cor_norExp_meth,na.rm=T),col="green")
	dev.off()
'
	png(file=paste0(outplace,"hist_corre_normExp_Meth_popRespGenes_vst_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(vec_cor_norExp_meth, main="Weighted methylation level\nvs. vst transformed normalized expression level\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=15,xlim=c(-1,1))
	abline(v=mean(vec_cor_norExp_meth,na.rm=T),col="green")
	dev.off()
'	
# how many percent has positive correlation above .5 ?
	NROW(which(vec_cor_norExp_meth > .5))/NROW((vec_cor_norExp_meth))
	#	[1] 0.1574803
	vec_named_cor_norExp_meth <- vec_cor_norExp_meth[which(vec_cor_norExp_meth > .5)] 	; names(vec_named_cor_norExp_meth) <- rownames(count_meth_sub)[which(vec_cor_norExp_meth > .5)]

	norm_counts_exp_sub_mean_CORSUB <- vst_counts_exp_sub_mean[which(vec_cor_norExp_meth > .5),]
	count_meth_sub_CORSUB <- count_meth_sub[which(vec_cor_norExp_meth > .5),]
'
	# For ALL genes in common:
	# coorelate rows
	vec_cor_norExp_meth_ALL <- sapply(1:nrow(count_meth_all_common), function(i) cor(t(count_meth_all_common[i,]), t(vst_counts_exp_all_common_mean[i,])))
	Non_vst_vec_cor_norExp_meth_ALL <- sapply(1:nrow(count_meth_all_common), function(i) cor(t(count_meth_all_common[i,]), t(norm_counts_exp_all_common_mean[i,])))

	png(file=paste0(outplace,"hist_corre_normExp_Meth_ALLGenes_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(Non_vst_vec_cor_norExp_meth_ALL, main="Weighted methylation level\nvs. normalized expression level\nfor all genes",xlab="Pearsons correlation Coefficient", breaks=15,xlim=c(-1,1))
	abline(v=mean(Non_vst_vec_cor_norExp_meth_ALL,na.rm=T),col="green")
	dev.off()
	
	png(file=paste0(outplace,"hist_corre_normExp_Meth_ALLGenes_vst_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(vec_cor_norExp_meth_ALL, main="Weighted methylation level\nvs. vst transformed normalized expression level\nfor all genes",xlab="Pearsons correlation Coefficient", breaks=15,xlim=c(-1,1))
	abline(v=mean(vec_cor_norExp_meth_ALL,na.rm=T),col="green")
	dev.off()
###################################################################################
###		Does higher methylation result in more stable expression?	 ##############
###################################################################################

# Calculate stability of expression by calculating variance
	norm_counts_exp_sub_sd <- aggregate(t(norm_counts_exp_sub[1:nrow(norm_counts_exp_sub),]), by=list(meta_exp$pop_temp), sd)
	vst_counts_exp_sub_sd <- aggregate(t(vst_counts_exp_sub[1:nrow(vst_counts_exp_sub),]), by=list(meta_exp$pop_temp), sd)

	rownames(norm_counts_exp_sub_sd) <- norm_counts_exp_sub_sd[,1]		; norm_counts_exp_sub_sd <- norm_counts_exp_sub_sd[,-1]
	norm_counts_exp_sub_sd <- data.frame(t(norm_counts_exp_sub_sd))
	
	rownames(vst_counts_exp_sub_sd) <- vst_counts_exp_sub_sd[,1]		; vst_counts_exp_sub_sd <- vst_counts_exp_sub_sd[,-1]
	vst_counts_exp_sub_sd <- data.frame(t(vst_counts_exp_sub_sd))
	
	# for all genes
	norm_counts_exp_all_common_sd <- aggregate(t(norm_counts_exp_all_common[1:nrow(norm_counts_exp_all_common),]), by=list(meta_exp$pop_temp), sd)
	vst_counts_exp_all_common_sd <- aggregate(t(vst_counts_exp_all_common[1:nrow(vst_counts_exp_all_common),]), by=list(meta_exp$pop_temp), sd)

	rownames(norm_counts_exp_all_common_sd) <- norm_counts_exp_all_common_sd[,1]		; norm_counts_exp_all_common_sd <- norm_counts_exp_all_common_sd[,-1]
	norm_counts_exp_all_common_sd <- data.frame(t(norm_counts_exp_all_common_sd))
	
	rownames(vst_counts_exp_all_common_sd) <- vst_counts_exp_all_common_sd[,1]		; vst_counts_exp_all_common_sd <- vst_counts_exp_all_common_sd[,-1]
	vst_counts_exp_all_common_sd <- data.frame(t(vst_counts_exp_all_common_sd))

	
# order methylation and expression data:
	norm_counts_exp_sub_sd <- norm_counts_exp_sub_sd[order(rownames(norm_counts_exp_sub_sd)),]
	vst_counts_exp_sub_sd <- vst_counts_exp_sub_sd[order(rownames(vst_counts_exp_sub_sd)),]
	count_meth_sub <- count_meth_sub[order(rownames(count_meth_sub)),]
	# for all genes
	norm_counts_exp_all_common_sd <- norm_counts_exp_all_common_sd[order(rownames(norm_counts_exp_all_common_sd)),]
	vst_counts_exp_all_common_sd <- vst_counts_exp_all_common_sd[order(rownames(vst_counts_exp_all_common_sd)),]
	count_meth_all_common <- count_meth_all_common[order(rownames(count_meth_all_common)),]

# correlate all rows, to see if higher methylation gives higher expression:
	NON_vst_vec_cor_norExp_sd_meth <- sapply(1:nrow(count_meth_sub), function(i) cor(t(count_meth_sub[i,]), t(norm_counts_exp_sub_sd[i,])))
		#these correlations look almost similar to the variance stabilized ones.
	vec_cor_norExp_sd_meth <- sapply(1:nrow(count_meth_sub), function(i) cor(t(count_meth_sub[i,]), t(vst_counts_exp_sub_sd[i,])))
	# for all genes:
	NON_vst_vec_cor_norExp_sd_meth_ALL <- sapply(1:nrow(count_meth_all_common), function(i) cor(t(count_meth_all_common[i,]), t(norm_counts_exp_all_common_sd[i,])))
		#these correlations look almost similar to the variance stabilized ones.
	vec_cor_norExp_sd_meth_ALL <- sapply(1:nrow(count_meth_all_common), function(i) cor(t(count_meth_all_common[i,]), t(vst_counts_exp_all_common_sd[i,])))

# histogram of correlations, check if right skewed. lower variance in expression with higher methylation is the idea.
	#pdf(file=paste0(outplace,"hist_corre_normExp_variance_Meth_popRespGenes_",fdr_thresh, ".pdf"), width=5, height=4)
	#par(cex=.7)
	#hist(vec_cor_norExp_sd_meth, main="Weighted methylation level\nvs. SD of Normalized expression level (log2foldchange)\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=20,xlim=c(-1,1))
	#abline(v=mean(vec_cor_norExp_sd_meth,na.rm=T),col="green")
	#dev.off()
	png(file=paste0(outplace,"hist_corre_normExp_variance_Meth_popRespGenes_vst_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(vec_cor_norExp_sd_meth, main="Weighted methylation level\nvs. SD of Normalized vst transformed expression level\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=20,xlim=c(-1,1))
	abline(v=mean(vec_cor_norExp_sd_meth,na.rm=T),col="green")
	dev.off()
	# all genes:
	png(file=paste0(outplace,"hist_corre_normExp_variance_Meth_ALLGenes_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(NON_vst_vec_cor_norExp_sd_meth_ALL, main="Weighted methylation level\nvs. SD of Normalized expression level\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=20,xlim=c(-1,1))
	abline(v=mean(NON_vst_vec_cor_norExp_sd_meth_ALL,na.rm=T),col="green")
	dev.off()
	png(file=paste0(outplace,"hist_corre_normExp_variance_Meth_ALLGenes_vst_",fdr_thresh, ".png"), width=5, height=4,units="in",res=400)
	par(cex=.7)
	hist(vec_cor_norExp_sd_meth_ALL, main="Weighted methylation level\nvs. SD of Normalized vst transformed expression level\nfor genes with population response",xlab="Pearsons correlation Coefficient", breaks=20,xlim=c(-1,1))
	abline(v=mean(vec_cor_norExp_sd_meth_ALL,na.rm=T),col="green")
	dev.off()

# I actually want to look at all the 1500 genes with population effect in methylation, but I don't care if they show population effect in expression.
# I will plot methylation level on x, and SD of expression on y, and color all points based on population. Lets see.

# subset expression counts files to genes with population effect in methylation
	norm_counts_exp_popMETH <- merge(norm_counts_exp, meth_dat_sub, by.x=0, by.y=0, all.y=T)[,1:200]	;	dim(norm_counts_exp_popMETH)
	vst_counts_exp_popMETH <- merge(vst_counts_exp, meth_dat_sub, by.x=0, by.y=0, all.y=T)[,1:200]		;	dim(vst_counts_exp_popMETH)
	meth_count_popMETH <- merge(count_meth, meth_dat_sub, by.x=0, by.y=0, all.y=T)[,1:21]		;	dim(meth_count_popMETH)
	
# add metadata groups to expression files
	rownames(norm_counts_exp_popMETH) <- norm_counts_exp_popMETH[,1]	; norm_counts_exp_popMETH <- norm_counts_exp_popMETH[,-1]
	norm_counts_exp_popMETH_SD <- aggregate(t(norm_counts_exp_popMETH[1:nrow(norm_counts_exp_popMETH),]), by=list(meta_exp$pop_temp), sd)
	rownames(norm_counts_exp_popMETH_SD) <- norm_counts_exp_popMETH_SD[,1]		; norm_counts_exp_popMETH_SD <- norm_counts_exp_popMETH_SD[,-1]
	norm_counts_exp_popMETH_SD <- data.frame(t(norm_counts_exp_popMETH_SD))

	rownames(vst_counts_exp_popMETH) <- vst_counts_exp_popMETH[,1]	; vst_counts_exp_popMETH <- vst_counts_exp_popMETH[,-1]
	vst_counts_exp_popMETH_SD <- aggregate(t(vst_counts_exp_popMETH[1:nrow(vst_counts_exp_popMETH),]), by=list(meta_exp$pop_temp), sd)
	rownames(vst_counts_exp_popMETH_SD) <- vst_counts_exp_popMETH_SD[,1]		; vst_counts_exp_popMETH_SD <- vst_counts_exp_popMETH_SD[,-1]
	vst_counts_exp_popMETH_SD <- data.frame(t(vst_counts_exp_popMETH_SD))	

# order methylation and expression data:
	norm_counts_exp_popMETH_SD <- norm_counts_exp_popMETH_SD[order(rownames(norm_counts_exp_popMETH_SD)),]
	vst_counts_exp_popMETH_SD <- vst_counts_exp_popMETH_SD[order(rownames(vst_counts_exp_popMETH_SD)),]
	rownames(meth_count_popMETH) <- meth_count_popMETH[,1]	; meth_count_popMETH <- meth_count_popMETH[,-1]
	meth_count_popMETH <- meth_count_popMETH[order(rownames(meth_count_popMETH)),]

'	 dim(meth_count_popMETH)
[1] 1537   25
> dim(vst_counts_exp_popMETH)
[1] 1537  199
> dim(norm_counts_exp_popMETH)
[1] 1537  199
'
library(plyr)
	# Plot limits
xlimsare=c(min(meth_count_popMETH,na.rm=T),max(meth_count_popMETH,na.rm=T))
ylimsare=c(min(vst_counts_exp_popMETH_SD,na.rm=T), round_any( max(vst_counts_exp_popMETH_SD,na.rm=T), .6,f=ceiling) )
ylimsare1=c(min(norm_counts_exp_popMETH_SD,na.rm=T),max(norm_counts_exp_popMETH_SD,na.rm=T))
	
#Create violin plot with these limits on methylation:
#(low [0–0.2]; medium [0.2–0.7]; high [0.7–1])

 # convert to long format
 viodat_EXPR <- data.frame(rownames(vst_counts_exp_popMETH_SD),vst_counts_exp_popMETH_SD)
 viodat_EXPR1 <- data.frame(rownames(norm_counts_exp_popMETH_SD),norm_counts_exp_popMETH_SD)
 viodat_METHYL <- data.frame(rownames(meth_count_popMETH),meth_count_popMETH)
 viodat_EXPR_long <- melt(viodat_EXPR)
 viodat_EXPR_long1 <- melt(viodat_EXPR1)
 viodat_METHYLlong <- melt(viodat_METHYL)
 viodata <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long$value)
 viodata1 <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long1$value)
 colnames(viodata) <- c("Methy_level","EXP_value")
 colnames(viodata1) <- c("Methy_level","EXP_value")
 
 
 # create groups based on Liu et al 2019
 viodata$Methy_level[ which(viodata$Methy_level > .7) ] <- "High"
 viodata$Methy_level[ which(viodata$Methy_level <= .7 & viodata$Methy_level >= .2) ] <- "Medium"
 viodata$Methy_level[ which(viodata$Methy_level < .2) ] <- "Low"
 viodata$Methy_level <- as.factor(viodata$Methy_level)
 viodata$Methy_level <- factor(viodata$Methy_level, levels=c("Low","Medium","High"))

viodata1$Methy_level[ which(viodata1$Methy_level > .7) ] <- "High"
 viodata1$Methy_level[ which(viodata1$Methy_level <= .7 & viodata1$Methy_level >= .2) ] <- "Medium"
 viodata1$Methy_level[ which(viodata1$Methy_level < .2) ] <- "Low"
 viodata1$Methy_level <- as.factor(viodata1$Methy_level)
 viodata1$Methy_level <- factor(viodata1$Methy_level, levels=c("Low","Medium","High"))


# test if difference
lm(EXP_value~Methy_level,data=viodata1)
'Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Coefficients:
      (Intercept)  Methy_levelMedium    Methy_levelHigh
            975.7             -615.4             -704.2
'
anova(lm(EXP_value~Methy_level,data=viodata1))
'Analysis of Variance Table

Response: EXP_value
               Df     Sum Sq    Mean Sq F value    Pr(>F)
Methy_level     2 3.1815e+09 1590774168  29.861 1.108e-13 ***
Residuals   29748 1.5847e+12   53272444
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(lm(EXP_value~Methy_level,data=viodata1))
'
Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Residuals:
   Min     1Q Median     3Q    Max
  -976   -963   -341   -211 635300

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)         975.73      65.49  14.898  < 2e-16 ***
Methy_levelMedium  -615.36     102.01  -6.033 1.63e-09 ***
Methy_levelHigh    -704.22     102.32  -6.883 5.99e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7299 on 29748 degrees of freedom
  (989 observations deleted due to missingness)
Multiple R-squared:  0.002004,  Adjusted R-squared:  0.001936
F-statistic: 29.86 on 2 and 29748 DF,  p-value: 1.108e-13
'
HSD.test(lm(EXP_value~Methy_level,data=viodata1),"Methy_level")$group
       'EXP_value groups
Low     975.7289      a
Medium  360.3644      b
High    271.5102      b
'
# Draw the plot

#pdf(file=paste0(outplace,"stability_of_expression_assay_violin.pdf"), width=7, height=5)
#	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("orange",.7), rectCol=NA, lineCol=NA, colMed="gray99", xlab="Weighted methylation level", ylab="SD expression (log2)") 
#dev.off()
#png(file=paste0(outplace,"stability_of_expression_assay_violin.png"), width=7, height=5,units="in",res=400)
#	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("orange",.7), rectCol=NA, lineCol=NA, colMed="gray99", xlab="Weighted methylation level", ylab="SD expression (log2)") 
#dev.off()

### Double plots: dot plot and violin plot side by side

#pdf(file=paste0(outplace,"stability_of_expression_assay_dotViolin_vst_",fdr_thresh, ".pdf"), width=7, height=3)
#	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,4,1,0), cex=.8)
#	ylabis="SD of expression\n(normalized vst transformed counts)"
#	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab=ylabis,xpd=NA)
#	for(rownr in 1:nrow(meth_count_popMETH)){
#		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare)
#	}
#	text(-.1, 5.3, labels="a", xpd=NA, cex=1.2)
#	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
#			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare, yaxt="n", border=alpha("darkorchid4",.7)) 
#	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
#	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
#	text(0.2, 5.3, labels="b", xpd=NA, cex=1.2)
#	text(c(1,2,3),rep(4.65,3),labels=c("a","b","b"))
#dev.off()
png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_vst_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,4,1,0), cex=.8)
	ylabis="SD of expression\n(normalized vst transformed counts)"
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab=ylabis, xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare)
	}
	text(-.1, 5.3, labels="a", xpd=NA, cex=1.2)
	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(0.2, 5.3, labels="b", xpd=NA, cex=1.2)
	text(c(1,2,3),rep(4.65,3),labels=c("a","b","b"))
dev.off()


### Double plots: SAME just NOT VST counts

png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_NONvst_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,4,0,0), cex=.8)
	ylabis="SD of expression\n(normalized counts)"
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare1, xlab="Weighted methylation level",ylab=ylabis, xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(norm_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare1)
	}
	vioplot( EXP_value ~ Methy_level, data=viodata1,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare1, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(c(1,2,3),rep(max(viodata1[,2]),3),labels=c("a","b","b"))
dev.off()

##################################################################
## vioplots for ALL genes:
##########
norm_counts_exp_all_common_sd <- norm_counts_exp_all_common_sd[order(rownames(norm_counts_exp_all_common_sd)),]
vst_counts_exp_all_common_sd <- vst_counts_exp_all_common_sd[order(rownames(vst_counts_exp_all_common_sd)),]
count_meth_all_common <- count_meth_all_common[order(rownames(count_meth_all_common)),]


	# Plot limits
xlimsare=c(min(count_meth_all_common,na.rm=T),max(count_meth_all_common,na.rm=T))
ylimsare=c(min(vst_counts_exp_all_common_sd,na.rm=T), round_any( max(vst_counts_exp_all_common_sd,na.rm=T), .6,f=ceiling) )
ylimsare1=c(min(norm_counts_exp_all_common_sd,na.rm=T),max(norm_counts_exp_all_common_sd,na.rm=T))
	
#Create violin plot with these limits on methylation:
#(low [0–0.2]; medium [0.2–0.7]; high [0.7–1])

 # convert to long format
viodat_EXPR <- data.frame(rownames(vst_counts_exp_all_common_sd),vst_counts_exp_all_common_sd)
viodat_EXPR1 <- data.frame(rownames(norm_counts_exp_all_common_sd),norm_counts_exp_all_common_sd)
viodat_METHYL <- data.frame(rownames(count_meth_all_common),count_meth_all_common)
viodat_EXPR_long <- melt(viodat_EXPR)
viodat_EXPR_long1 <- melt(viodat_EXPR1)
viodat_METHYLlong <- melt(viodat_METHYL)
viodata <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long$value)
viodata1 <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long1$value)
colnames(viodata) <- c("Methy_level","EXP_value")
colnames(viodata1) <- c("Methy_level","EXP_value")
 
 
 # create groups based on Liu et al 2019
viodata$Methy_level[ which(viodata$Methy_level > .7) ] <- "High"
viodata$Methy_level[ which(viodata$Methy_level <= .7 & viodata$Methy_level >= .2) ] <- "Medium"
viodata$Methy_level[ which(viodata$Methy_level < .2) ] <- "Low"
viodata$Methy_level <- as.factor(viodata$Methy_level)
viodata$Methy_level <- factor(viodata$Methy_level, levels=c("Low","Medium","High"))

viodata1$Methy_level[ which(viodata1$Methy_level > .7) ] <- "High"
viodata1$Methy_level[ which(viodata1$Methy_level <= .7 & viodata1$Methy_level >= .2) ] <- "Medium"
viodata1$Methy_level[ which(viodata1$Methy_level < .2) ] <- "Low"
viodata1$Methy_level <- as.factor(viodata1$Methy_level)
viodata1$Methy_level <- factor(viodata1$Methy_level, levels=c("Low","Medium","High"))


# test if difference
lm(EXP_value~Methy_level,data=viodata1)
'Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Coefficients:
      (Intercept)  Methy_levelMedium    Methy_levelHigh  
            825.0             -399.6             -439.9 
'
anova(lm(EXP_value~Methy_level,data=viodata1))
'Analysis of Variance Table

Response: EXP_value
                Df     Sum Sq    Mean Sq F value    Pr(>F)    
Methy_level      2 2.7072e+10 1.3536e+10  164.63 < 2.2e-16 ***
Residuals   652072 5.3615e+13 8.2223e+07                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(lm(EXP_value~Methy_level,data=viodata1))
'
Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Residuals:
    Min      1Q  Median      3Q     Max 
   -825    -824    -804    -350 1496899 

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)         825.03      13.90   59.38   <2e-16 ***
Methy_levelMedium  -399.64      36.96  -10.81   <2e-16 ***
Methy_levelHigh    -439.93      26.83  -16.40   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 9068 on 652072 degrees of freedom
  (47345 observations deleted due to missingness)
Multiple R-squared:  0.0005047, Adjusted R-squared:  0.0005016 
F-statistic: 164.6 on 2 and 652072 DF,  p-value: < 2.2e-16
'
HSD.test(lm(EXP_value~Methy_level,data=viodata1),"Methy_level")$group
'       EXP_value groups
Low     825.0336      a
Medium  425.3951      b
High    385.1001      b
'
# Draw the plot
### Double plots: dot plot and violin plot side by side

png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_ALLgenes_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,4,1,0), cex=.8)
	ylabis="SD of expression\n(normalised vst transformed counts)"
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab=ylabis,xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare)
	}
	text(-.1, ylimsare[2]+(.09*(ylimsare[2]-ylimsare[1])), labels="a", xpd=NA, cex=1.2)
	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
		xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=seq(ylimsare[1], ylimsare[2],1),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(0.2, ylimsare[2]+(.09*(ylimsare[2]-ylimsare[1])), labels="b", xpd=NA, cex=1.2)
	text(c(1,2,3),rep(max(viodata[,2]),3),labels=c("a","b","b"))
dev.off()


### Double plots: SAME just NOT VST counts

png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_NONvst_ALLgenes_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,4,0,0), cex=.8)
	ylabis="SD of expression\n(normalised counts)"
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare1, xlab="Weighted methylation level",ylab=ylabis,xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(norm_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare1)
	}
	text(-.1, ylimsare[2]+(.09*(ylimsare[2]-ylimsare[1])), labels="a", xpd=NA, cex=1.2)
	vioplot( EXP_value ~ Methy_level, data=viodata1,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare1, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=seq(ylimsare[1], ylimsare[2],1),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(0.2, max(norm_counts_exp_popMETH_SD)+(.09*(max(norm_counts_exp_popMETH_SD)-0)), labels="b", xpd=NA, cex=1.2)
	text(c(1,2,3),rep(max(viodata1[,2]),3),labels=c("a","b","b"))
dev.off()

#############################################

##################################################################
## vioplots for genes subset to only include genes where expression yields a phenotypic effect
##########
pattern_fit_files <- list.files(path="/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/downstream_analyses/analysis_pheno_intermediary_mechanisms/", pattern="subset_fits_patterns.csv", full.names=T )
phenotype_like_cases <- lapply(pattern_fit_files[grep("Exp", pattern_fit_files)], read.csv)

### GOT TO HERE
# subset data to genes with similarity to phenotype
norm_counts_exp_all_common_sd <- norm_counts_exp_all_common_sd[order(rownames(norm_counts_exp_all_common_sd)),]
vst_counts_exp_all_common_sd <- vst_counts_exp_all_common_sd[order(rownames(vst_counts_exp_all_common_sd)),]
count_meth_all_common <- count_meth_all_common[order(rownames(count_meth_all_common)),]


	# Plot limits
xlimsare=c(min(count_meth_all_common,na.rm=T),max(count_meth_all_common,na.rm=T))
ylimsare=c(min(vst_counts_exp_all_common_sd,na.rm=T), round_any( max(vst_counts_exp_all_common_sd,na.rm=T), .6,f=ceiling) )
ylimsare1=c(min(norm_counts_exp_all_common_sd,na.rm=T),max(norm_counts_exp_all_common_sd,na.rm=T))
	
#Create violin plot with these limits on methylation:
#(low [0–0.2]; medium [0.2–0.7]; high [0.7–1])

 # convert to long format
viodat_EXPR <- data.frame(rownames(vst_counts_exp_all_common_sd),vst_counts_exp_all_common_sd)
viodat_EXPR1 <- data.frame(rownames(norm_counts_exp_all_common_sd),norm_counts_exp_all_common_sd)
viodat_METHYL <- data.frame(rownames(count_meth_all_common),count_meth_all_common)
viodat_EXPR_long <- melt(viodat_EXPR)
viodat_EXPR_long1 <- melt(viodat_EXPR1)
viodat_METHYLlong <- melt(viodat_METHYL)
viodata <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long$value)
viodata1 <- data.frame(viodat_METHYLlong$value, viodat_EXPR_long1$value)
colnames(viodata) <- c("Methy_level","EXP_value")
colnames(viodata1) <- c("Methy_level","EXP_value")
 
 
 # create groups based on Liu et al 2019
viodata$Methy_level[ which(viodata$Methy_level > .7) ] <- "High"
viodata$Methy_level[ which(viodata$Methy_level <= .7 & viodata$Methy_level >= .2) ] <- "Medium"
viodata$Methy_level[ which(viodata$Methy_level < .2) ] <- "Low"
viodata$Methy_level <- as.factor(viodata$Methy_level)
viodata$Methy_level <- factor(viodata$Methy_level, levels=c("Low","Medium","High"))

viodata1$Methy_level[ which(viodata1$Methy_level > .7) ] <- "High"
viodata1$Methy_level[ which(viodata1$Methy_level <= .7 & viodata1$Methy_level >= .2) ] <- "Medium"
viodata1$Methy_level[ which(viodata1$Methy_level < .2) ] <- "Low"
viodata1$Methy_level <- as.factor(viodata1$Methy_level)
viodata1$Methy_level <- factor(viodata1$Methy_level, levels=c("Low","Medium","High"))


# test if difference
lm(EXP_value~Methy_level,data=viodata1)
'Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Coefficients:
      (Intercept)  Methy_levelMedium    Methy_levelHigh  
            825.0             -399.6             -439.9 
'
anova(lm(EXP_value~Methy_level,data=viodata1))
'Analysis of Variance Table

Response: EXP_value
                Df     Sum Sq    Mean Sq F value    Pr(>F)    
Methy_level      2 2.7072e+10 1.3536e+10  164.63 < 2.2e-16 ***
Residuals   652072 5.3615e+13 8.2223e+07                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(lm(EXP_value~Methy_level,data=viodata1))
'
Call:
lm(formula = EXP_value ~ Methy_level, data = viodata1)

Residuals:
    Min      1Q  Median      3Q     Max 
   -825    -824    -804    -350 1496899 

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)         825.03      13.90   59.38   <2e-16 ***
Methy_levelMedium  -399.64      36.96  -10.81   <2e-16 ***
Methy_levelHigh    -439.93      26.83  -16.40   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 9068 on 652072 degrees of freedom
  (47345 observations deleted due to missingness)
Multiple R-squared:  0.0005047, Adjusted R-squared:  0.0005016 
F-statistic: 164.6 on 2 and 652072 DF,  p-value: < 2.2e-16
'
HSD.test(lm(EXP_value~Methy_level,data=viodata1),"Methy_level")$group
'       EXP_value groups
Low     825.0336      a
Medium  425.3951      b
High    385.1001      b
'
# Draw the plot
### Double plots: dot plot and violin plot side by side

pdf(file=paste0(outplace,"stability_of_expression_assay_dotViolin_ALLgenes_",fdr_thresh, ".pdf"), width=7, height=3)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,3.5,1,0), cex=.8)
		plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)",xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare)
	}
	text(-.1, 5.3, labels="a", xpd=NA, cex=1.2)
	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(0.2, 5.3, labels="b", xpd=NA, cex=1.2)
	text(c(1,2,3),rep(max(viodata[,2]),3),labels=c("a","b","b"))
dev.off()
png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_ALLgenes_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,3.5,1,0), cex=.8)
		plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)",xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare)
	}
	text(-.1, 5.3, labels="a", xpd=NA, cex=1.2)
	vioplot( EXP_value ~ Methy_level, data=viodata,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(0.2, 5.3, labels="b", xpd=NA, cex=1.2)
	text(c(1,2,3),rep(max(viodata[,2]),3),labels=c("a","b","b"))
dev.off()


### Double plots: SAME just NOT VST counts

pdf(file=paste0(outplace,"stability_of_expression_assay_dotViolin_NONVST_ALLgenes_",fdr_thresh, ".pdf"), width=7, height=3)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,3.5,0,0), cex=.8)
		plot(NA,NA, xlim=xlimsare,ylim=ylimsare1, xlab="Weighted methylation level",ylab="SD of expression (log2)",xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(norm_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare1)
	}
	vioplot( EXP_value ~ Methy_level, data=viodata1,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare1, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(c(1,2,3),rep(max(viodata1[,2]),3),labels=c("a","b","b"))
dev.off()
png(file=paste0(outplace,"stability_of_expression_assay_dotViolin_NONVST_ALLgenes_",fdr_thresh, ".png"), width=7, height=3,units="in",res=400)
	par(mfrow=c(1,2),mar=c(5,1,1,1),oma=c(0,3.5,0,0), cex=.8)
		plot(NA,NA, xlim=xlimsare,ylim=ylimsare1, xlab="Weighted methylation level",ylab="SD of expression (log2)",xpd=NA)
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(norm_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.2), xlim=xlimsare,ylim=ylimsare1)
	}
	vioplot( EXP_value ~ Methy_level, data=viodata1,  	col=alpha("darkorchid4",.7), rectCol=NA, lineCol=NA, colMed="gray99", pchMed="-",
			xlab="Weighted methylation level", ylab=NA, xpd=NA, ylim=ylimsare1, yaxt="n", border=alpha("darkorchid4",.7)) 
	axis(2, at=c(0,1,2,3,4),xpd=NA, labels=NA)
	axis(1, at=c(1,2,3), labels=c("Low","Medium","High"),xpd=NA)
	text(c(1,2,3),rep(max(viodata1[,2]),3),labels=c("a","b","b"))
dev.off()

#############################################

'
# dot plot alone:
pdf(file=paste0(outplace,"stability_of_expression_assay_dotPLOT.pdf"), width=7, height=5)
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)")
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.3), xlim=xlimsare,ylim=ylimsare)
	}
dev.off()
png(file=paste0(outplace,"stability_of_expression_assay_dotPLOT.png"), width=7, height=5,units="in",res=400)
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)")
	for(rownr in 1:nrow(meth_count_popMETH)){
		points(t(meth_count_popMETH[rownr,]),t(vst_counts_exp_popMETH_SD[rownr,]), col=alpha(popcols,.3), xlim=xlimsare,ylim=ylimsare)
	}
dev.off()

# plot regression line plots.

pdf(file=paste0(outplace,"stability_of_expression_assay_regLINE.pdf"), width=7, height=5)
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)")
	for(rownr in 1:nrow(meth_count_popMETH)){
		if( all(is.na(t(vst_counts_exp_popMETH_SD[rownr,]))) ){ next }
		test_is <- lm(t(vst_counts_exp_popMETH_SD[rownr,]) ~ t(meth_count_popMETH[rownr,]),na.action=na.omit)
		if( test_is$coefficients[2] >= .2){coloris="blue"}
		if( test_is$coefficients[2] < .2 & test_is$coefficients[2] > -.2 ){coloris="yellow"}
		if( test_is$coefficients[2] <= -.2){coloris="red"}
		abline(lm(t(vst_counts_exp_popMETH_SD[rownr,]) ~ t(meth_count_popMETH[rownr,]),na.action=na.omit), col=alpha(coloris,.1))
	}
dev.off()
png(file=paste0(outplace,"stability_of_expression_assay_regLINE.png"), width=7, height=5,units="in",res=400)
	plot(NA,NA, xlim=xlimsare,ylim=ylimsare, xlab="Weighted methylation level",ylab="SD of expression (log2)")
	for(rownr in 1:nrow(meth_count_popMETH)){
		if( all(is.na(t(vst_counts_exp_popMETH_SD[rownr,]))) ){ next }
		test_is <- lm(t(vst_counts_exp_popMETH_SD[rownr,]) ~ t(meth_count_popMETH[rownr,]),na.action=na.omit)
		if( test_is$coefficients[2] >= .2){coloris="blue"}
		if( test_is$coefficients[2] < .2 & test_is$coefficients[2] > -.2 ){coloris="yellow"}
		if( test_is$coefficients[2] <= -.2){coloris="red"}
		abline(lm(t(vst_counts_exp_popMETH_SD[rownr,]) ~ t(meth_count_popMETH[rownr,]),na.action=na.omit), col=alpha(coloris,.1))
	}
dev.off()

'

#######################################################################################################
###		Does EXPRESSION group according to populations acroos genes with Interaction effect?	 	##############
#######################################################################################################

nrow(exp_dat_sub_int)

#subset counts file to genes with int effect
	vst_counts_pca_int <- merge(vst_counts_exp, exp_dat_sub_int, by.x=0, by.y=1, all.y=T)[,1:200]
	rownames(vst_counts_pca_int) <- vst_counts_pca_int[,1]	; vst_counts_pca_int <- vst_counts_pca_int[,-1]

	
pca_popExp_genes_int <- prcomp(t(na.omit(vst_counts_pca_int)))
summa <- summary(pca_popExp_genes_int)
summa
'Importance of components:
                           PC1      PC2      PC3     PC4     PC5     PC6
Standard deviation     35.1091 15.69500 15.55729 12.4244 8.24586 8.00601
Proportion of Variance  0.4735  0.09462  0.09297  0.0593 0.02612 0.02462
Cumulative Proportion   0.4735  0.56811  0.66108  0.7204 0.74650 0.77112
                           PC7    PC8     PC9    PC10    PC11    PC12    PC13
Standard deviation     7.60654 6.7111 5.58779 5.33368 5.09090 4.72738 4.23600
Proportion of Variance 0.02223 0.0173 0.01199 0.01093 0.00996 0.00858 0.00689
Cumulative Proportion  0.79334 0.8106 0.82264 0.83357 0.84352 0.85211 0.85900
'


meta_exp_dat <- data.frame(rownames(meta_exp),meta_exp$population) 	; colnames(meta_exp_dat) <- c("rna_nr","population")
	# change rna nr to be similar to the other df (RNA_0001) not RNA 1
meta_exp_dat$rna_nr_new <- sapply(strsplit(as.character(meta_exp_dat$rna_nr)," "), function(x)
    paste(x[1], sprintf("%04d",as.numeric(x[2])), sep="_")
)

meta_exp_dat$population <- gsub("E", "O",meta_exp_dat$population)
pop_order <- vector()
for( each in meta_exp_dat$population ){	pop_order <- append(pop_order,grep(each,names(popcols)))	}
popcols2 <- popcols[pop_order]
#output_type="pdf"
#output_type="png"


if (	output_type == "pdf"	)	{ 	
	pdf(paste0(outplace,"pca_EXPint_plots_",fdr_thresh, ".pdf"),width=7,height=7) 		
}else if ( 	output_type == "png"	){
	png(paste0(outplace,"pca_EXPint_plots_",fdr_thresh, ".png"),width=7,height=7,units="in",res=400) 	
}
par(mar=c(4.7,4.2,3.7,1.5), mfrow=c(2,2),cex=.7)
plottit <- ordiplot(pca_popExp_genes_int,choices=c(1,2),type="n", 
			xlim=c(-80,80),ylim=c(-80,80), cex=0.8, main="Expression, Interaction",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"))
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))
legend("topright",legend=sort(names(popcols)),col=popcols[sort(names(popcols))],pch=1, 
			inset=c(-(.42),-.21), xpd=NA, ncol=2)
			
plottit <- ordiplot(pca_popExp_genes_int,choices=c(2,3),type="n", 
			xlim=c(-60,60),ylim=c(-60,60), cex=0.8, main="Expression, Interaction",
			xlab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"))
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))

screeplot(pca_popExp_genes_int, type="l", main="Screeplot Expression, Interaction")
biplot(pca_popExp_genes_int, main="Biplot Expression, Interaction")
dev.off()



