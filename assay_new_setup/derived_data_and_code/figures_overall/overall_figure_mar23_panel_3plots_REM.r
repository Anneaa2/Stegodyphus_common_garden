#!/home/anneaa/miniconda3/envs/assay_copy/bin/Rscript

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

fig_types=	"individuel_panels"# "overview" # "individuel_panels"
#output_type="pdf" # "png" # "pdf"
output_type="pdf"

library(agricolae)
library(colorspace)
library(scales)
library(eulerr)
library(ggplot2)
library(multcomp)
library(emmeans)
library(multcompView)
library(plyr)

#popcols=c("steelblue1","blue","tan4","green3","snow4","red") # B, K, N, S, G, E/O
popcols=c("lightskyblue1","blue2","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS modify lightblue (steelblue1)
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","snow4","red2") # B, K, N, S, G, E/O    #### OBS new modified popcols
popcols <- popcols[c(1,2,6,4)]	# B, K, E/O,	S
coltest=F  # logical
coltest_type= "desaturate"  # tritan	#deutan		#protan		#desaturate		#
proportion_down_from_max= 0.015 # for moving the significance groups down from the maxY limit. for phenotype plots.
proportion_down_for_xlabels_text = 0.15 # for moving the significance groups down from the maxY limit. for phenotype plots.
proportion_to_increase_ylim= 0.085 # for increasing the ymax of all phenotype plots
main_line_out=1

#par( mfrow=c(4,3))
col_vec =  c("red", "orange", "yellow","green",
				"turquoise","blue","purple","magenta4","maroon","tomato3")
				
opaque_col1 <- function(coloris, number=4){ ### ADD NUMBER TO THIS FUNCTION!
	if (number==3){alphas=rev(c(255,135,60))}
	if (number==4){alphas=rev(c(255,135,60,30))}
	colo <- rgb(t(col2rgb(rep(coloris,number))), max = 255, alpha = alphas)
	return(colo)
}

# make fake data for data-less plots
custom_dat = as.data.frame(t(c(10,10,10,10)))
colnames(custom_dat) <- rev(c("Population","Temperature", "Pop+Temp","Interaction"))

if (coltest) {
	if (coltest_type == "desaturate") { 	popcols		<- desaturate(popcols)	; coltest_type = "desaturate_" }
	if (coltest_type == "tritan") { 	popcols		<- tritan(popcols)		;	coltest_type = "tritan_" }
	if (coltest_type == "protan") { 	popcols		<- protan(popcols)		; coltest_type = "protan_" }
	if (coltest_type == "deutan") { 	popcols		<- deutan(popcols)		; coltest_type = "deutan_" }
}else{	coltest_type = ""	}	



#------------- Overall figure info

if (	output_type == "pdf"	)	{ 	
    pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"overview_fig_mar23_panel3plots1.pdf"),width=7,height=10.5) 		# a4 page is 8.3x11.7
}else if ( 	output_type == "png"	){
    png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"overview_fig_mar23_panel3plots1.png"),width=7,height=10.5,units="in",res=400) 	
}
par(mfrow = c(1,3))
#lay_mat <- matrix(c(		1,1,1,		2,2,2,		3,3,3,			4,4,4,		
 #                           5,5,5,		6,6,6,		7,7,7,			8,8,8,			
  #                      9,9,9,9,		10,10,10,10,		11,11,11,11,		
   #                     12,12,12,		13,13,13,		14,14,14,		15,15,15),
    #                nrow=4,ncol=12,byrow=T)

#layout(lay_mat)


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
####
#	Expression data
###	
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
overlap_genes_exprT_methP <- NROW(intersect(exp_dat_Temp[,1],meth_dat_pop_sign[,1]))
overlap_genes_exprT_expP_methP <- NROW(intersect(	intersect(exp_dat_Temp[,1],meth_dat_pop_sign[,1])	,	exp_dat_pop[,1]))

####
#	Expression data
###

overlap_genes_pop_temp <- NROW(intersect(exp_dat_pop[,1],exp_dat_Temp[,1]))
message("Ekspression: Number of genes overlapping between population and temperaure related genes: ", overlap_genes_pop_temp)
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
#------------- Expression / RNA plotting  ----#
###############################################

meth_dat_plot <- meth_dat; meth_dat_plot[ which(meth_dat_plot == 0) ] <- NA
combi_dat_frame <- rbind(rev(meth_dat_plot),rev(datis))
combi_dat_frame <- combi_dat_frame[,c(1,3,4)]
combi_dat_frame <- rbind(c(0,0,0),c(0,0,0),combi_dat_frame)
	#empty rows
combi_overlap_frame <- rbind(rev(c(overlap_genes_exprP_methP,0,0,0)),rev(c(overlap_genes_exprP_methP,0,0,0)))
combi_overlap_frame <- combi_overlap_frame[,c(1,3,4)]
combi_overlap_frame <- combi_overlap_frame/overlap_genes_exprP_methP_TOTAL*100
combi_overlap_frame <- rbind(c(0,0,0),c(0,0,0),combi_overlap_frame)
	#empty rows

new_xmax <- ceiling(max(combi_dat_frame,na.rm=T)/10)*10
barppp <- barplot(combi_dat_frame, beside=T, horiz=T, xpd=NA, axes=F, xlab="Number of genes", xlim=c(0,new_xmax))
barplot(combi_overlap_frame, beside=T, horiz=T, col="black", density=17, angle=135, xpd=NA, add=T, xaxt="n")
#axis(1, at=c(0,seq(2000,14000,by=2000)), labels=c(0,2000,NA,6000,NA,10000,NA,14000),xpd=NA)
axis(1, at=seq(0,new_xmax,by=10), labels=seq(0,new_xmax,by=10),xpd=NA)
mtext("a",side=3,line=0.5, at=max(barppp) * .05,cex=1.2)
title(main="Expression & Methylation", xpd=NA, line=2.9)
#dev.off()


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
names(metab_dat_sub)<-NULL
metab_dat <- metab_dat/c(NMR_ORG_total,NMR_AQ_total,NMR_ORG_total,NMR_AQ_total,LCMS_total,LCMS_total)
metab_dat <- metab_dat*100

metab_dat_nmr <- metab_dat[ grep("NMR",rownames(metab_dat)), ]
sums_nmr_ctmax <- colSums(metab_dat_nmr[ grep("CTmax", rownames(metab_dat_nmr)), ])
sums_nmr_ccr <- colSums(metab_dat_nmr[ grep("CCR", rownames(metab_dat_nmr)), ])
metab_dat_nmr_org <- metab_dat_nmr[ grep("ORG",rownames(metab_dat_nmr)), ]
metab_dat_lcms <- metab_dat[ grep("LCMS",rownames(metab_dat)), ]


metab_data_new <- rbind(sums_nmr_ctmax, sums_nmr_ccr, metab_dat_lcms)
rownames(metab_data_new)[c(1,2)] <- c("NMR_CTmax","NMR_CCR")
metab_data_new_add <- rbind( metab_dat_nmr_org, rep(0,3), rep(0,3))

new_xmax <- ceiling(max(metab_data_new)/10)*10
plotaa <- barplot(as.matrix(rev(metab_data_new)), beside=T, horiz=T, xpd=NA, axes=F, 
	xlim=c(0,new_xmax), xlab="Number of metabolites/peaks")
axis(1, at=seq(0,new_xmax,10), labels=seq(0,new_xmax,10))
barplot(as.matrix(rev(metab_data_new_add)), beside=T, horiz=T, col="black", density=17, angle=135/2, xpd=NA, add=T, xaxt="n", xlim=c(0,round(max(metab_data_new),digits=-2)))
mtext("b",side=3,line=0.5, at=max(plotaa) * .05,cex=1.2)
title(main="Metabolism", xpd=NA, line=2.9)
#dev.off()


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
micro_dat <- rbind(c(0,0,0),c(0,0,0),c(0,0,0),micro_dat)
	# EMPTY ROWS

ploooot <- barplot(as.matrix(rev(micro_dat)),      beside=T, horiz=T, xlab="% of ASVs", xlim=c(0,ceiling(max(micro_dat))))
mtext("c",side=3,line=0.5, at=max(ploooot) * .05,cex=1.2)
title(main="Microbiome", xpd=NA, line=2.9)
dev.off()
