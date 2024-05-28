#!/home/anneaa/miniconda3/envs/assay_copy/bin/Rscript


###      OBS add to one scale / plot


# source("assay_new_setup/derived_data_and_code/figures_overall/overall_figure.r")
# srun --mem=10g --pty /bin/bash
## Overall figure with patterns of adaptation and plasticity

# row 1: 	Methylations			RNA						Micro
														# yaxis: c( 	"Temperature effect",	"population effect",	"interaction"	)
# row 2: 	Matabolites-ctmax		Matabolites-CCR		Matabolites-Tobias
														# yaxis: c( 	"Temperature effect",	"population effect",	"interaction"	)
# Row 3-6:				Effect of temp			effect of popq
# Row 3		Ctmax				
# Row 4		CCR			
# Row 5		Growth				
# Row 6		Survival
	
	# if you make new env, error will most likely happen. This solves it (from bash):  cp miniconda3/envs/YOUR_ENV/lib/liblapack.so miniconda3/envs/YOUR_ENV/lib/libRlapack.so
#conda activate assay_copy; R			
#setwd("assay_new_setup/derived_data_and_code/figures_overall")

#output_type="pdf" # "png" # "pdf"
output_type="pdf"

library(viridis)
library(stringr)
library(agricolae)
library(colorspace)
library(scales)
library(eulerr)
library(ggplot2)
library(multcomp)
library(emmeans)
library(multcompView)
library(plyr)

coltest=F  # logical
coltest_type= "desaturate"  # tritan	#deutan		#protan		#desaturate		#
	# desaturate(popcols)
if (coltest) {
	if (coltest_type == "desaturate") { 	popcols <- desaturate(popcols)	; coltest_type = "desaturate_" }
	if (coltest_type == "tritan") { 	popcols <- tritan(popcols)		;	coltest_type = "tritan_" }
	if (coltest_type == "protan") { 	popcols	<- protan(popcols)		; coltest_type = "protan_" }
	if (coltest_type == "deutan") { 	popcols	<- deutan(popcols)		; coltest_type = "deutan_" }
}else{	coltest_type = ""	}	



####################################
#------------- Methylations	-------#
####################################

fdr_threshold = 0.1
	# number of genes showing patterns
meth_file = "assay_new_setup/derived_data_and_code/methylation/smps/gene_region/DSS_output/DSS_result_all_d10D32_gene_CPG.txt"
meth_datis <- read.table(meth_file,	header=T,	stringsAsFactors=T)
meth_total_genes <- nrow(meth_datis)
meth_dat_pop_sign <- meth_datis[which(meth_datis$fdr_p < fdr_threshold), ]
meth_dat_temp_sign <- meth_datis[which(meth_datis$fdr_T < fdr_threshold), ]
meth_dat <- c(nrow(meth_dat_pop_sign),nrow(meth_dat_temp_sign),NROW(intersect( droplevels(meth_dat_pop_sign$region_name), droplevels(meth_dat_temp_sign$region_name) )),NA)
names(meth_dat) <- c("Population","Temperature","Pop+Temp","Interaction")
meth_dat
meth_dat <- meth_dat/meth_total_genes*100
meth_dat

###	######################
#	Expression data			####
###	######################
	#for expression data (also fdr - why two different cutoffs) earlier: 0.05
padj_threshold = 0.1
	# import methylation data
	# import expression data:
exp_file = "assay_new_setup/derived_data_and_code/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt"
exp_dat <- read.table(exp_file, header=T,       row.names=1, stringsAsFactors=T)
exp_total_genes <- nrow(exp_dat)

exp_dat_pop <-	exp_dat[	which( exp_dat[,which( colnames(exp_dat) == "padj_P" )] < padj_threshold	)	,	]
exp_dat_Temp <- exp_dat[	which( exp_dat[,which( colnames(exp_dat) == "padj_T" )] < padj_threshold	)	,	]

overlap_genes_exprP_methP <- NROW(intersect(exp_dat_pop[,1],meth_dat_pop_sign[,1]))
overlap_genes_exprP_methP_TOTAL <- NROW(intersect(exp_dat[,1],meth_datis[,1]))
message("Ekspression and methylation: Number of genes overlapping between EXSP-population and METH-population related genes: ", overlap_genes_exprP_methP)
'
Ekspression and methylation: Number of genes overlapping between EXSP-population and METH-population related genes: 844
'
#overlap_genes_exprT_methP <- NROW(intersect(exp_dat_Temp[,1],meth_dat_pop_sign[,1]))
#overlap_genes_exprT_expP_methP <- NROW(intersect(	intersect(exp_dat_Temp[,1],meth_dat_pop_sign[,1])	,	exp_dat_pop[,1]))

####
#	Expression data
###

overlap_genes_pop_temp <- NROW(intersect(exp_dat_pop[,1],exp_dat_Temp[,1]))
#message("Ekspression: Number of genes overlapping between population and temperaure related genes: ", overlap_genes_pop_temp)
'
Ekspression: Number of genes overlapping between population and temperaure related genes: 7669
'

datis <- c(	NROW(which( exp_dat$padj_P < padj_threshold )),
			NROW(which( exp_dat$padj_T < padj_threshold )),
			overlap_genes_pop_temp,
			NROW(which( exp_dat$padj_I < padj_threshold ))	)
datis
datis <- datis/exp_total_genes*100
datis
'
[1] 13998 11807  7669   909
'
###############################################
#------------- Overlap Methylation / Expression ----#
###############################################

meth_dat_plot <- meth_dat ; meth_dat_plot[ which(is.na(meth_dat_plot)) ] <- 0
combi_dat_frame <- rbind(rev(meth_dat_plot),rev(datis))
combi_dat_frame <- combi_dat_frame[,c(1,3,4)]
combi_overlap_frame <- rbind(rev(c(overlap_genes_exprP_methP,0,0,0)),rev(c(overlap_genes_exprP_methP,0,0,0)))
combi_overlap_frame <- combi_overlap_frame[,c(1,3,4)]
combi_overlap_frame <- combi_overlap_frame/c(meth_total_genes,exp_total_genes)*100


####################################
#------------- Metabolites ## 2  plot ## From tobias LCMS
####################################

NMR_AQ_total <- read.csv("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CTmax.txt")
NMR_AQ_total <- NMR_AQ_total[-c(1:4)]	; NMR_AQ_total <- ncol(NMR_AQ_total)
	# similar for AQ CCR:	"assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CCR.txt"
NMR_ORG_total <- read.csv("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CTmax.txt")
NMR_ORG_total <- NMR_ORG_total[-c(1:4)]	; NMR_ORG_total <- ncol(NMR_ORG_total)
	# similar for ORG CCR:	"assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CCR.txt"
LCMS_total <- read.csv("assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv")
LCMS_total <- LCMS_total[-c(1:2)]	; LCMS_total <- ncol(LCMS_total)
	# similar for CCR: 	"assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv"

'NMR_AQ_total
[1] 454
> NMR_ORG_total
[1] 240
> LCMS_total
[1] 95 '

#grep data for CTmax samples:
#metab_file = "assay_new_setup/derived_data_and_code/metabolites/data/Significance_table_ANOVA_metabolites_LCMS.csv"

metab_file = "assay_new_setup/derived_data_and_code/metabolites/data/Significance_table_ANOVA_metabolites_ALL.csv"
metab_dat <- read.csv(metab_file,	header=T, stringsAsFactors=T, row.names=1)
metab_dat <- metab_dat[,-ncol(metab_dat)]
metab_dat <- metab_dat[,-3]
#names(metab_dat)<-NULL
metab_dat <- metab_dat/c(NMR_ORG_total,NMR_AQ_total,NMR_ORG_total,NMR_AQ_total,LCMS_total,LCMS_total)
metab_dat <- metab_dat*100
metab_dat <- metab_dat[c(1,3,2,4,5,6),]
#metab_dat_nmr <- metab_dat[ grep("NMR",rownames(metab_dat)), ]

#sums_nmr_ctmax <- colSums(metab_dat_nmr[ grep("CTmax", rownames(metab_dat_nmr)), ])
#sums_nmr_ccr <- colSums(metab_dat_nmr[ grep("CCR", rownames(metab_dat_nmr)), ])
#metab_dat_nmr_org <- metab_dat_nmr[ grep("ORG",rownames(metab_dat_nmr)), ]
#metab_dat_lcms <- metab_dat[ grep("LCMS",rownames(metab_dat)), ]

#metab_data_new <- rbind(sums_nmr_ctmax, sums_nmr_ccr, metab_dat_lcms)
#rownames(metab_data_new)[c(1,2)] <- c("NMR_CTmax","NMR_CCR")
#rownames(metab_data_new) <- paste0("Metabolites\n", rownames(metab_data_new))
#rownames(metab_data_new) <- gsub("CCR","CCRTemp",rownames(metab_data_new))
#rownames(metab_data_new) <- gsub("_","-",rownames(metab_data_new))

#metab_data_new_add <- rbind( metab_dat_nmr_org, rep(0,3), rep(0,3))



####################################
#------------- Microbiome 1 plot---#
####################################
	# number of ASVs showing patterns (after cutoff of 25% prevalence in at least one pop/temp)

total_number_micro_file = "assay_new_setup/derived_data_and_code/microbiome/data_output/aov_betweenSbjs_all_ASVs.csv"
total_number_micro <- read.csv(total_number_micro_file,	header=T,	row.names=1, stringsAsFactors=T)
total_number_micro <- nrow(total_number_micro)
micro_file = "assay_new_setup/derived_data_and_code/microbiome/data_output/Count_sign_ASVs.csv"
micro_dat <- read.csv(micro_file,	header=T,	row.names=1, stringsAsFactors=T)
micro_dat <- micro_dat[-3]
micro_dat <- micro_dat/total_number_micro*100

#combi_dat <- rbind(rev(t(micro_dat)), as.matrix(rev(metab_data_new)), combi_dat_frame)
#rownames(combi_dat) <- c("Microbiome", rownames(metab_data_new), "DNA Methylation","Gene Expression")
combi_dat <- rbind(rev(t(micro_dat)), as.matrix(rev(metab_dat)), combi_dat_frame)
rownames(combi_dat) <- c("Microbiome", rownames(metab_dat), "DNA methylation","Gene expression")

'                intsig_vec tempsig_vec popsig_vec
Microbiome        1.282051     0.00000   7.692308
NMR_ORG_CTmax     0.000000    39.58333   0.000000
NMR_AQ_CTmax      0.000000    21.80617  17.400881
NMR_ORG_CCR       0.000000    26.66667   0.000000
NMR_AQ_CCR        0.000000    24.66960  30.176211
LCMS_CTmax       42.105263    47.36842  33.684211
LCMS_CCR          4.210526    22.10526  27.368421
DNA Methylation   0.000000     0.00000   4.395070
Gene Expression   2.433446    31.40076  37.227733'

#combi_dat_overlap <- rbind(c(0,0,0),as.matrix(rev(metab_data_new_add)),combi_overlap_frame)
#rownames(combi_dat_overlap) <- c("Micro", rownames(metab_data_new), "Methylation","Expression")
combi_dat_overlap <- rbind(c(0,0,0),c(0,0,0),c(0,0,0),c(0,0,0),c(0,0,0),c(0,0,0),c(0,0,0),combi_overlap_frame)
rownames(combi_dat_overlap) <- c("Micro", rownames(metab_dat), "methylation","expression")

'             [,1] [,2]     [,3]
Micro            0    0 0.000000
NMR_ORG_CTmax    0    0 0.000000
NMR_AQ_CTmax     0    0 0.000000
NMR_ORG_CCR      0    0 0.000000
NMR_AQ_CCR       0    0 0.000000
LCMS_CTmax       0    0 0.000000
LCMS_CCR         0    0 0.000000
Methylation      0    0 2.413428
Expression       0    0 2.244621'

############## Plotting ##############

#combi_dat_ORG <- combi_dat_overlap	; combi_dat_ORG[,3] <- 0
#combi_dat_overlap[,2] <- 0
	# for plotting it twice



N_for_legend_vector <- c( exp_total_genes, meth_total_genes, LCMS_total, LCMS_total, NMR_AQ_total, NMR_AQ_total, NMR_ORG_total, NMR_ORG_total, total_number_micro)
names(N_for_legend_vector) <- rev(rownames(combi_dat))


# trying new multiple groups plotting to add spacing between groups of grouped bars
combi_dat1 <- rbind(combi_dat[1,],rep(NA,3), combi_dat[2:7,], rep(NA,3),combi_dat[8,], rep(NA,3), combi_dat[9,])
combi_dat_overlap1 <- rbind(combi_dat_overlap[1,],rep(NA,3), combi_dat_overlap[2:7,], rep(NA,3),combi_dat_overlap[8,], rep(NA,3),combi_dat_overlap[9,])
#combi_dat_ORG1 <- rbind(combi_dat_ORG[1,],rep(NA,3), combi_dat_ORG[2:5,], rep(NA,3),combi_dat_ORG[6,], rep(NA,3),combi_dat_ORG[7,])
	# Add N to rownames:
rownames(combi_dat) <- paste0(rownames(combi_dat), " (N: ", rev(N_for_legend_vector), ")")
	# Add empty spaces for empty data series
rownames(combi_dat1) <- c(rownames(combi_dat)[1],"", rownames(combi_dat)[2:7], "",rownames(combi_dat)[8], "",rownames(combi_dat)[9])
rownames(combi_dat_overlap1) <- c(rownames(combi_dat_overlap)[1],"", rownames(combi_dat_overlap)[2:7], "",rownames(combi_dat_overlap)[8], "",rownames(combi_dat_overlap)[9])
#rownames(combi_dat_ORG1) <- c(rownames(combi_dat_ORG)[1],"", rownames(combi_dat_ORG)[2:5], "",rownames(combi_dat_ORG)[6], "",rownames(combi_dat_ORG)[7])
combi_dat_overlap1[combi_dat_overlap1 == 0] <- NA
#combi_dat_ORG1[combi_dat_ORG1 == 0] <- NA




color_pattern <- c(	viridis(1, begin = 0, end = .03), "white", 		
					viridis(6, begin = .23, end = .5), "white",		
					viridis(1, begin = .8, end = .85), "white",		
					viridis(1, begin = 1, end = 1))

#color_pattern <- desaturate( color_pattern )
new_xmax <- ceiling(max(combi_dat1,na.rm=T)/10)*10
space_is <- 2.5
legend_col_vector <- c(1:8,2,9:12)
legend_name_vector <- c(rownames(combi_dat1)[1:8],"Metabolites:",rownames(combi_dat1)[9:12])
legend_name_vector <- str_pad(legend_name_vector, width = max(nchar(legend_name_vector)), side="right", pad=" ")	# space pad names
legend_name_vector <- gsub("_","-",legend_name_vector)
legend_name_vector <- gsub("ORG","Org",legend_name_vector)
legend_name_vector <- gsub("AQ","Aq",legend_name_vector)
annotation_col_vec <- color_pattern
annotation_col_vec[which(annotation_col_vec != "white")] <- "black"	# vector for the pop and temp term
annotation_col_vec[which(annotation_col_vec == "white")] <- "NA"
annotation_col_vec1 <- annotation_col_vec ;	color_pattern1 <- color_pattern	# vector for the interaction term
annotation_col_vec1[10] <- NA ;	color_pattern1[10] <- "white"


width_is = 6.8		; height_is = 5.1
if (	output_type == "pdf"	)	{ 	
    pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"overview_fig_mar23_onefig1.pdf"),width=width_is,height=height_is) 		# a4 page is 8.3x11.7
}else if ( 	output_type == "png"	){
    png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"overview_fig_mar23_onefig1.png"),width=width_is,height=height_is,units="in",res=600) 	
}
par(mfrow = c(1,1), mar=c(2.3,2.9,.3,7.1), family = "Helvetica")
barppp <- barplot(combi_dat1, beside=T, horiz=T, xpd=NA, axes=F, mgp=c(3,1.7,0),	# mend 2. in mgp to move axis labels on y
	col=color_pattern, xlim=c(0,new_xmax), space=c(0,space_is),names.arg=c("Interaction","Temperature","Population"))
barplot(combi_dat_overlap1, beside=T, horiz=T, col="black", add=T, axisnames = F, space=c(0,space_is), 
		density=20, angle=45, axes=F) # add hashing
calc_string_widht_legend <- max(strwidth(legend_name_vector))
#par(lty=3, lwd=2) # Make dots wider		#barplot(combi_dat_ORG1, beside=T, horiz=T, col="grey", xpd=NA, add=T, axisnames =F, border=NA, space=c(0,space_is),		#		density=11, angle=20, axes=F) # Add dots		#par(lty=1, lwd=1)
#barplot(combi_dat_ORG1, beside=T, horiz=T, col="grey", xpd=NA, add=T, axisnames =F, border=NA, space=c(0,space_is),		#		density=15, angle=45, axes=F) # Add white lines
#barplot(combi_dat_ORG1, beside=T, horiz=T, xpd=NA, add=T, axisnames=F, space=c(0,space_is),		#		border="black", col=alpha("white",0), axes=F) # Add border
axis(1, at=seq(0,new_xmax,by=10), labels=seq(0,new_xmax,by=10),xpd=NA, mgp=c(3,.25,0),tcl=-.25, cex.axis=.8, line=-.5)
title(xlab="Percent of genes / metabolites / ASVs", line=1)
#text(x=combi_dat1[8,1]+1.5, y=barppp[8,1], labels="NA", cex=.8, xpd=NA, adj=0)	#	Add legend as lines of text alligned to graph
abline(h=barppp[nrow(barppp),1:2]+((1+space_is)/2), lty=2, col="grey")
par(cex=.8)
#legend("topright", inset=c(-.3,-.017), ncol=1, xpd=NA, fill=rev(color_pattern[legend_col_vector]), text.col="black", 
#		legend=rev(legend_name_vector), y.intersp=.78, border=NA, bg="white", box.col="white")
#legend("topright", inset=c(-.3,-.017), ncol=1, xpd=NA, fill=alpha("white",0), text.col=alpha("white",0), legend=rev(legend_name_vector), #, inset=c(-.05,-.017)
#		y.intersp=.78, border= rev(c("black", NA, rep("black",6), rep(NA,2), "black", NA, "black")), bty='n')	# add box
		# If I should decrease space between the blank legend elements, I think the easiest is plotting several.
# Add gene expression
inset_y <- -.011		; inset_move <- 0.037		; 	inset_x <- -.37
legend("topright", inset=c(inset_x,inset_y), ncol=1, xpd=NA, fill=rev(color_pattern[legend_col_vector])[1], text.col="black", 
		legend=rev(legend_name_vector)[1], y.intersp=.78, border="black", bg=NA, box.col=NA, text.width=calc_string_widht_legend)
# Add gene expression
legend("topright", inset=c(inset_x,inset_y+inset_move), ncol=1, xpd=NA, fill=rev(color_pattern[legend_col_vector])[3], text.col="black", 
		legend=rev(legend_name_vector)[3], y.intersp=.78, border="black", bg=NA, box.col=NA, text.width=calc_string_widht_legend)
# Add metabolites
legend("topright", inset=c(inset_x,inset_y+(inset_move*2)), ncol=1, xpd=NA, fill=rev(color_pattern[legend_col_vector])[5:11], text.col="black", 
		legend=rev(legend_name_vector)[5:11], y.intersp=.78, border=c("white",rep("black",6)), bg=NA, box.col=NA, text.width=calc_string_widht_legend)
# add microbiome:
legend("topright", inset=c(inset_x,inset_y+(inset_move*7.5)), ncol=1, xpd=NA, fill=rev(color_pattern[legend_col_vector])[13], text.col="black", 
		legend=rev(legend_name_vector)[13], y.intersp=.78, border="black", bg=NA, box.col=NA, text.width=calc_string_widht_legend)
par(cex=1)
rect(xleft=rep(-2.5,12), xright=rep(-1.5,12),ybottom=barppp[,3]-.5, ytop=barppp[,3]+.5, col=color_pattern, border=annotation_col_vec, xpd=NA)
rect(xleft=rep(-2.5,12), xright=rep(-1.5,12),ybottom=barppp[,2]-.5, ytop=barppp[,2]+.5, col=color_pattern, border=annotation_col_vec, xpd=NA)
rect(xleft=rep(-2.5,12), xright=rep(-1.5,12),ybottom=barppp[,1]-.5, ytop=barppp[,1]+.5, col=color_pattern1, border=annotation_col_vec1, xpd=NA) #annotation_col_vec1
dev.off()










