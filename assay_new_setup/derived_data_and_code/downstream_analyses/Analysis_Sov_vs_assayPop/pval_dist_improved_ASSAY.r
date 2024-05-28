#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript

#srun --mem=8g --pty /bin/bash
#conda activate sources_of_var; R
#source("figures_for_paper/pval_dist_improved.r")

library(viridis)
library(VennDiagram)
library(ggplot2)
library(eulerr)
library(cowplot)
library(gridExtra)
library(magick)

#####libraries and data aquisition:
#setwd("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/")
fileout="pdf" #pdf or png
analysis_type="mantel" # mantel or mrm
breaks=50
source("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/make_data_lists_ASSAY.r")
	# makes the lists of data
percentile_is = .99
source("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/count_genes_in_percentile_nulldist_ASSAY.r")	#modify to fit mrm	
	# counts and gets the genes in the percentile of null dist

# list_name_overlap - gives the gene name
	# Make venn/euler diagrams?
#Prepare data: 	change names to be identical.
	# if smps genes names has a "_no" appended, it means it was run using simple mantel, since there were no snp variation to correct for. Removed in list_name_overlap_even
list_name_overlap_even <- list()
for (	i in 1:NROW(list_name_overlap) ){
	list_name_overlap_even[[i]] 	<- 	gsub("gene_snp_|gene_smp_|relative_|absolute_|_no","",x=list_name_overlap[[i]])
}
names(list_name_overlap_even) <- names(list_name_overlap)
	# for no good reason
dataset_touse <- list_name_overlap_even 	

	# get the number of genes (and names of them) that correlate strongly in Assay and in Sov respectively
assay_numbers <- grep("ASSAY" ,	names(dataset_touse)	)
sov_numbers <- grep("SOV" ,	names(dataset_touse)	)
assay_genes <- unique(sort(unlist(dataset_touse[assay_numbers])))
sov_genes <- unique(sort(unlist(dataset_touse[sov_numbers])))
NROW(assay_genes)
NROW(sov_genes)

	# Check overlap between Assay and Sov strongly correlating genes (across all env axes)
NROW(intersect(assay_genes, sov_genes))
	# percent of sov genes, also found in assay
NROW(intersect(assay_genes, sov_genes))/NROW(sov_genes)*100
	# if percentile:	overlap percentage: overlap nr:
	#		.9				90					19957
	#		.98				71					6166	
	#		.99				53					2789
	#		.999			5					25
	#		.9999			0					0





















































ovelap_dataframe 	<- data.frame(	 paste0("PC",axisnumberis),	snpN,		refpopN,					smpN,			overlap_snp_refpop,		overlap_snp_smp,		ovelap_smp_refpop,
										ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)


# make overlap table:
for ( axisnumberis in 1:5 ){
	list_elm_togrep 	<- grep( paste0("PC",axisnumberis),	names(dataset_touse)	)
	greplist <- names(dataset_touse)[list_elm_togrep]
	
		#smp <- list_elm_togrep[	grep("smp", greplist)	]
		
		taxa_rel <-  list_elm_togrep[	grep("genuslevel_relative_", greplist)	]		;		taxa_abs <-  list_elm_togrep[	grep("genuslevel_absolute_", greplist)	]
		ASV_rel <- list_elm_togrep[	grep("singleASV_relative_", greplist)	]			;		ASV_abs <-  list_elm_togrep[	grep("singleASV_absolute_", greplist)	]
		
		snpN <-	NROW(dataset_touse[[snp]])		;	refpopN <- NROW(dataset_touse[[refpop]])		;		smpN <- NROW(dataset_touse[[smp]])	
		taxa_relN <- 	NROW(dataset_touse[[taxa_rel]])		;		taxa_absN <-  	NROW(dataset_touse[[taxa_abs]])
		ASV_relN <- 	NROW(dataset_touse[[ASV_rel]])			;		ASV_absN <-  	NROW(dataset_touse[[ASV_abs]])
		
		ovelap_smp_refpop 		<- NROW(intersect(dataset_touse[[smp]],  dataset_touse[[refpop]]))
		overlap_snp_smp			<- NROW(intersect(dataset_touse[[snp]],  dataset_touse[[smp]]))
		#overlap_snp_refpop		<- NROW(intersect(dataset_touse[[snp]],  dataset_touse[[refpop]]))
		overlap_taxa_abs_rel	<- NROW(intersect(dataset_touse[[taxa_rel]],  dataset_touse[[taxa_abs]]))
		overlap_ASV_abs_rel		<- NROW(intersect(dataset_touse[[ASV_rel]],  dataset_touse[[ASV_abs]]))
		
		if ( paste0("PC",axisnumberis) == "PC1"	){
			ovelap_dataframe 	<- data.frame(	 paste0("PC",axisnumberis),	snpN,		refpopN,					smpN,			overlap_snp_refpop,		overlap_snp_smp,		ovelap_smp_refpop,
										ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)
		}else{
			ovelap_dataframe	<- rbind(ovelap_dataframe, data.frame(	 paste0("PC",axisnumberis),	snpN,		refpopN,					smpN,			overlap_snp_refpop,		overlap_snp_smp,		ovelap_smp_refpop,
																ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)	)
		}
}

colnames(ovelap_dataframe) <- c("env_parameter",	"snp_N",	"smp_refpop_N",		"smp_cisref_N",		"overlap_snp_refpop_N",		"overlap_snp_cisref_N",		"overlap_smp_refpop_cisref_N",
								"ASV_rel_N",	"ASV_abs_N",	"ASV_overlap_abs_rel_N",	"Genera_rel_N",	"Genera_abs_N",	"Genera_overlap_rel_abs_N")
write.csv(file= paste0("./figures_for_paper/COR_Coef_overlapDat_PCaxes_",analysis_type,dataset_is,".csv"),ovelap_dataframe)


# make white space figure:	for padding/space between plots. Found a better way.
	#	pdf(file="figures_for_paper/sum_up_plots/white_space.pdf", bg="white")	;plot.new()	;dev.off()

# PLOT euler venn 

colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[4],colorsare1[5],colorsare1[5])
widthis=7
#emptywidth= (widthis/5)* 0.09  	# making room for empty plots when combining the panel (they will be 10% of one of the five plots. Adjust in this line.
#newwidth= (widthis/5)- 2*(emptywidth) #redundant method now.
newwidth=widthis/5
margin_padding=.0
set.seed(10)
list_plots <- list()
for ( axisnumberis in 1:5 ){
		if (analysis_type=="mrm") {
			euler_vector <- c(	"A"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "snp_N") ]	,
								"B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_cisref_N") ],
								"A&B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_cisref_N") ]
								)
			colorvec_is=colorsare[1:2]
		}else{
			euler_vector <- c(	"A"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "snp_N") ]	,
								"B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_refpop_N") ],
								"C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_cisref_N") ],
								"A&B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_refpop_N") ],
								"A&C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_cisref_N") ],
								"B&C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_smp_refpop_cisref_N") ],
								"A&B&C"= NROW(intersect( intersect(dataset_touse[[snp]],  dataset_touse[[refpop]]) , dataset_touse[[smp]]))
								)
			colorvec_is=colorsare[1:3]
		}
		euler_diag <- euler(euler_vector, input="union")
		
		plottits1 <- plot(	euler_diag,	quantities = list(alpha=1, cex=.7), 				fill=colorvec_is,
				labels=NA,				alpha=.6	,adjust_labels=T ,	edges=list(alpha=.01, col="white")	)	#	main=paste0("PC ",axisnumberis),
		plottits <- arrangeGrob(plottits1, ncol=1,			padding=unit(margin_padding,"in"),top="",bottom="",right="",left="")
		list_plots	<- append(list_plots, list(plottits))
		
	
	if( COEF_ONLY==T ){
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_diagramm_jan22_coef_PC",axisnumberis,".pdf"), 
			plot=plottits,	height=heights_ofplots[3],			width=newwidth,units="in")
	#	ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_diagramm_jan22_coef_PC",axisnumberis,".png"), 
	#		plot=plottits, 	height=heights_ofplots[3],			width=newwidth,units="in",dpi=800)
	}
}

















# FROM NOW USING THE LISTS OF PVALUES AND COR COEF
	######## Decision was to run on coefficients only.
COEF_ONLY=T #; 	SIGN_ONLY=F ;		COEF_AND_SIG=F

##### fUNCTIONs FOR PLOTTING GRAPH:

source("assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop//plot_graph_func_3plots_ASSAY.r")	




# I do not believe I need to plot.
# Just get correct lists and check overlap between SOV and Assay treatments











##### Now to actual plotting:

#fileout="pdf"
#fileout="png"
heights_ofplots <- c(5,3,2)  # in combined figure ( coef, barplot, venn)
expected_theoretical_number=T # if I want to use theoretically expected number of genes/ ASVs
colorsare1 <- c(		plasma(1, begin=.05,end=.05,alpha=.95),		plasma(1, begin=.7,end=.7),		plasma(1, begin=.7,end=.7, alpha=.4), plasma(1, begin=.35,end=.35),		plasma(1, begin=.83,end=.83))
col_lines <- plasma(1, begin=.99,end=.99)		## COL2 colors
Density_is <- c(0,0,0,0,30,0,30)	;		angle_is <- c(0,0,0,0,45,0,45)	
if( COEF_ONLY==T ) {dataset_is <- paste0("_percentile_",percentile_is) }
#if( SIGN_ONLY==T ) {dataset_is <- paste0("_SignOnly")	}
#if( COEF_AND_SIG==T ) {dataset_is <- paste0("_SignAnd_percentile_",percentile_is)}

if (analysis_type == "mrm" ){ #it is mrm
	mainis=c("MRM on distance in\ngenewise SNPs" , "MRM on distance in\ngenewise WML" , "MRM on distance in\nMicrobiome")
	heightis=heights_ofplots[2] 	;	widthis=7
	if (fileout=="pdf"){
		pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
			"Figure_s14_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,".pdf"),
			height=heightis,width=widthis)
	}
	if (fileout=="png"){
		png(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
			"Figure_s14_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,".png"),
			height=heightis,width=widthis,units="in",res=600)
	}
	colorsare <- c(colorsare1[1:2],colorsare1[4],colorsare1[5],colorsare1[4],colorsare1[5])
	factor_for_micro_leg_inset = 5.147059; Density_is	<- Density_is[c(1,2,4,6,5,7)] ; angle_is <- angle_is[c(1,2,4,6,5,7)]
	annotateis <- c("a","b","c")	;	leg_inset=.40	;	leg_inset_micro=leg_inset/factor_for_micro_leg_inset	#.068		#
	par(mfrow=c(1,3),oma=c(7,.5,0,1.2),mar=c(0.1,4.1,4.1,1.1),mgp=c(2.8,1,0))
} else { # it is mantel
	mainis=c("Genewise SNPs",		"Genewise WML",	"Microbiome Relative Abundance")	;		heightis=heights_ofplots[2]	;	widthis=7
	if (fileout=="pdf"){
		pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
			"Figure_5bcd_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,".pdf"),
			height=heightis,width=widthis)
	}
	#if (fileout=="png"){
	#	png(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
	#		"Figure_5bcd_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,".png"),
	#		height=heightis,width=widthis,units="in",res=800)
	#}
	colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[4],colorsare1[5],colorsare1[5])
	factor_for_micro_leg_inset = 5.147059
	annotateis <- c("b","c","d")	;	leg_inset=.43	;	leg_inset_micro=leg_inset/factor_for_micro_leg_inset	#.068		#
	par(mfrow=c(1,3),oma=c(7,.5,0,.5),mar=c(0.1,4.1,3.6,1.1),mgp=c(2.8,1,0));		adjustis=-.3 # adjustis adjusts nothing as we speak.?
}
if ( analysis_type == "mantel" ) { plot_graph_now(); dev.off() }else{	plot_graph_now_mrm(); dev.off()	}


' Data for theoretical expected numbers of genes/ASV/Genera
dataset_exp_above_perctile_theory
       snp refpop   smps singleASV_relative singleASV_absolute
PC1 2.9393  3.183 3.1789              0.006             0.0057
PC2 2.9393  3.183 3.1789              0.006             0.0057
PC3 2.9393  3.183 3.1789              0.006             0.0057
PC4 2.9393  3.183 3.1789              0.006             0.0057
PC5 2.9393  3.183 3.1789              0.006             0.0057
    genuslevel_relative genuslevel_absolute
PC1              0.0043              0.0039
PC2              0.0043              0.0039
PC3              0.0043              0.0039
PC4              0.0043              0.0039
PC5              0.0043              0.0039
>dataset_N_correlations
      snp refpop  smps singleASV_relative singleASV_absolute
PC1 29393  31830 31789                 60                 57
PC2 29393  31830 31789                 60                 57
PC3 29393  31830 31789                 60                 57
PC4 29393  31830 31789                 60                 57
PC5 29393  31830 31789                 60                 57
    genuslevel_relative genuslevel_absolute
PC1                  43                  39
PC2                  43                  39
PC3                  43                  39
PC4                  43                  39
PC5                  43                  39


#### MRM:

> dataset_exp_above_perctile_theory 
      snp   smp singleasvs_relative singleasvs_absolute genuslevel_relative
PC1 2.9393 3.183               0.006              0.0057              0.0043
PC2 2.9393 3.183               0.006              0.0057              0.0043
PC3 2.9393 3.183               0.006              0.0057              0.0043
PC4 2.9393 3.183               0.006              0.0057              0.0043
PC5 2.9393 3.183               0.006              0.0057              0.0043
    genuslevel_absolute
PC1              0.0039
PC2              0.0039
PC3              0.0039
PC4              0.0039
PC5              0.0039
>dataset_N_correlations
      snp   smp singleasvs_relative singleasvs_absolute genuslevel_relative
PC1 29393 31830                  60                  57                  43
PC2 29393 31830                  60                  57                  43
PC3 29393 31830                  60                  57                  43
PC4 29393 31830                  60                  57                  43
PC5 29393 31830                  60                  57                  43
    genuslevel_absolute
PC1                  39
PC2                  39
PC3                  39
PC4                  39
PC5                  39

'

	
# Test color blindness  ## function imported in earlier script (source( functionfor...))

colname_is1 <- "_Fin" #"ORIGI" # COL1   COL2
colorblind_test()
'colorsare <- c("red2","blue","lightskyblue","green4","green4","lawngreen","lawngreen")#,"green3","lawngreen")#"green3")
colorsare <- rgb(t(col2rgb(colorsare)),alpha=120,maxColorValue=255)
'

source("figures_for_paper/plot_graph_func_3plots.r")	


# plot only microbiome:
#fileout="png"
#fileout="pdf"
expected_theoretical_number=T # if I want to use theoretically expected number of genes/ ASVs
mainis=c("Genewise SNPs",		"Genewise WML",	"Microbiome ASVs/Genera")	;		heightis=heights_ofplots[2]	;	widthis=3
if( COEF_ONLY==T ) {dataset_is <- paste0("_percentile_",percentile_is) }
#if( SIGN_ONLY==T ) {dataset_is <- paste0("_SignOnly")	}
#if( COEF_AND_SIG==T ) {dataset_is <- paste0("_SignAnd_percentile_",percentile_is)}
if (fileout=="pdf"){
		pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
			"Figure_s13_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_onlyMicro.pdf"),
			height=heightis,width=widthis)
	}
	if (fileout=="png"){
		png(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/sum_up_plots/",
			"Figure_s13_","COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_onlyMicro.png"),
			height=heightis,width=widthis,units="in",res=600)
	}
colorsare1 <- c(		plasma(1, begin=.05,end=.05,alpha=.95),		plasma(1, begin=.7,end=.7),		plasma(1, begin=.7,end=.7, alpha=.4), 	plasma(1, begin=.35,end=.35),		plasma(1, begin=.83,end=.83))
col_lines <- plasma(1, begin=.99,end=.99)		## COL2 colors
colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[4],colorsare1[5],colorsare1[5])
colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[5],colorsare1[4],colorsare1[5])
Density_is <- c(0,0,0,0,0,30,30)	;		angle_is <- c(0,0,0,0,0,45,45)	; 	factor_for_micro_leg_inset =  19.44444
annotateis <- c("b","c","d")	;	leg_inset=.46	;	leg_inset_micro=leg_inset/factor_for_micro_leg_inset	#.018
par(mfrow=c(1,1),oma=c(7.5,.5,0,2),mar=c(0.1,4.1,4.1,1.1),cex=.6)

plot_graph_now_microonly();		dev.off()




# Plotting the alternative distribution (the dist of cor coef from correlations between a randomized environmental axis).

colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[4],colorsare1[5],colorsare1[5])
# alternative distribution (corrected), real distribution (corrected), NO means, REDUCED NUMBER OF GRAPHS SHOWN

if (fileout=="pdf"){
	pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/","Figure_5a_","corcoef_distributions","_",analysis_type,"_RANDeach_reduced.pdf"),
			height=heights_ofplots[1],width=7)
}
#if (fileout=="png"){
#	png(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/","Figure_5a_","corcoef_distributions","_",analysis_type,"_RANDeach_reduced.png"),	height=heights_ofplots[1],width=7,units="in",res=1000)
#}
par(mfrow=c(3,5),oma=c(4,3.5,4,3),mar=c(1,1,1,1),xpd=NA)
for( datadist in names(listofdists)[grep("RealCor_COR", names(listofdists))] ){
	message(datadist)
	secdatadist_name <- sapply(strsplit(datadist,"RealCor_COR_"),"[[",2)
	names_togrepfrom <- names(listofdists)[	grep(secdatadist_name, names(listofdists))]
	header_name <- secdatadist_name
	for( secdatadist in names_togrepfrom[	grep("RandCor_COR",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
		message(secdatadist)
		if( grepl("absolute",secdatadist) & grepl("singleASV",secdatadist) ){ next } # skips microbiome ASV absolute
		if( grepl("genus",secdatadist) || grepl("smps",secdatadist) ){ next } #skips microbiome genuslevel and smps cis ref
		if( grepl("absolute",secdatadist) || grepl("relative",secdatadist) ){
			breaknumber=13	; maxisY <- 20 		}else{			breaknumber=150	; maxisY <- 800		} 
		h <- hist((listofdists[[secdatadist]]),breaks=breaknumber,plot=F)
		hs <- hist((listofdists[[datadist]]),breaks=breaknumber,plot=F)
		# get density function
			multiplier <- h$counts / h$density
			mydensity <- density(listofdists[[secdatadist]])
			mydensity$y <- mydensity$y * multiplier[1]
		colnameis <- sapply(strsplit(header_name,"_"),"[[",2)	;	rownameis <- sapply(strsplit(header_name,"_"),"[[",1)
		if( grepl("relative",secdatadist) ){		rownameis="Microbiome ASVs\nRel. abundance\nPop structure\ncorrection"	; color_new=colorsare[4]
			xaxisname="Correlation Coefficient" ; labelslogical=T		;		yaxisticks=c(0,5,10,15,20)	;	labelstosee <- c("0",NA,"10",NA,"20")
			coloris <- "black"		;		trancolorperc = 60		; 	trancolorperc_grey = 70 
		}else{		xaxisname=NA	;	coloris <- "grey69" 	;		trancolorperc = 80		; 	trancolorperc_grey = 90
			labelslogical=F		;	yaxisticks=c(0,200,400,600,800)		;	labelstosee <- c("0",NA,"400",NA,"800")	}
		if( grepl("refpop",secdatadist) ){			rownameis="WML\nPop structure\ncorrection"; color_new=colorsare[2]	}
		if( grepl("snp",secdatadist) ){		rownameis="SNPs\nPop structure\ncorrection"	; color_new=colorsare[1]		}
		plot(h, 	col=t_col("grey44",percent=trancolorperc_grey),	border=t_col("grey44",percent=trancolorperc_grey),		ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n", 	yaxt="n",	main=NA,	xlab=NA,		ylab=NA , 	cex.lab=1.3,	font.lab=1)
		lines(mydensity,col="grey60")
		plot(hs, 	col=t_col(color_new,percent=trancolorperc),			border=t_col(color_new,percent=trancolorperc),		ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n",	yaxt="n",	main=NA,	xlab=NA,		ylab=NA,	add=T)
		if ( grepl("PC5",secdatadist) ){				axis(4,	at=yaxisticks, labels=labelstosee)				}
		if ( grepl("PC5",secdatadist) & grepl("refpop",secdatadist) ) {title(ylab="Frequency", font.lab=1,cex.lab=1,outer=F, line= -10.5)	}
		if ( grepl("PC3",secdatadist) & grepl("relative",secdatadist) ) {mtext("Correlation Coefficient", side=1,	line=3,	font=1,cex=.7, adj=.3)	}
		if ( grepl("PC1",secdatadist) ){	title(ylab=rownameis, font.lab=1,cex.lab=1, line=0, adj=0)	}#mtext(rownameis, side=2, line=2, cex=.7,	las=2, 	font=2, padj=1.2) 	} # padj moves the text downward 
		axis(1,	at=c(-1,-.5,0,.5,1), labels=labelslogical, col=coloris)
		if( grepl("snp",secdatadist) ){	title(main=	sapply(strsplit(header_name,"_"),"[[",2)	, line=1,	cex.main=1)	}
		if( grepl("snp",secdatadist) & grepl("PC3",secdatadist)	){
			legend(-3.6,1400, legend=c(	"Expected null distribution (Permuted environmental axes)",
										"Actual distribution (Actual environmental axes)"),	
				fill=c(t_col("grey44",percent=trancolorperc_grey-15),t_col(color_new,percent=trancolorperc-15)),
				ncol=1, bty="n", cex=1,	border=c(t_col("grey44",percent=trancolorperc_grey-15),t_col(color_new,percent=trancolorperc-15)))
			#mtext("a",side=3,line=3.5,adj=-3.3,cex=1.2,outer=F)				
}	}	}	
dev.off()



'
# alternative distribution (NOT corrected), real distribution (corrected), NO means REDUCED NUMBER OF GRAPHS SHOWN

pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/corcoef_distributions","_",analysis_type,"_RANDeach_NOCORRECT_reduced.pdf"),
			height=heights_ofplots[1],width=7)
par(mfrow=c(3,5),oma=c(4,9,5,5),mar=c(1,1,1,1),xpd=NA)
for( datadist in names(listofdists)[grep("RealCor_COR", names(listofdists))] ){
	message("\n",datadist)
	secdatadist_name <- sapply(strsplit(datadist,"RealCor_COR_"),"[[",2)
	names_togrepfrom <- names(listofdists)[	grep(secdatadist_name, names(listofdists))]
	header_name <- secdatadist_name# sapply(strsplit(secdatadist_name,"Cor_"),"[[",2)
	for( secdatadist in names_togrepfrom[	grep("RandCor_NO",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
		message(secdatadist)
		if( grepl("absolute",secdatadist) & grepl("singleASV",secdatadist) ){ next } # skips microbiome ASV absolute
		if( grepl("genus",secdatadist) || grepl("smps",secdatadist) ){ next } #skips microbiome genuslevel and smps cis ref
		#if( grepl("PC4",secdatadist) || grepl("PC5",secdatadist) ){ next } #skips PC4 and PC5
		if( grepl("absolute",secdatadist) || grepl("relative",secdatadist) ){
			breaknumber=13	; maxisY <- 20 
		}else{
			breaknumber=150	; maxisY <- 800
		} 
		h <- hist((listofdists[[secdatadist]]),breaks=breaknumber,plot=F)
		hs <- hist((listofdists[[datadist]]),breaks=breaknumber,plot=F)		#;		maxisYy <- max(c(h$counts,hs$counts)) ; 	print(maxisYy) 
		# get density function
			multiplier <- h$counts / h$density
			mydensity <- density(listofdists[[secdatadist]])
			mydensity$y <- mydensity$y * multiplier[1]
		colnameis <- sapply(strsplit(header_name,"_"),"[[",2)	;	rownameis <- sapply(strsplit(header_name,"_"),"[[",1)
		message(rownameis) 
		#if( grepl("PC1",secdatadist) ){		yaxisname="Frequency"	}else{ 		yaxisname=NA	 }  #if first columns 
		if( grepl("relative",secdatadist) ){		rownameis="Microbiome ASVs\nRel. abundance\nPop structure\ncorrection"
			xaxisname="Correlation Coefficient" ; labelslogical=T		;		yaxisticks=c(0,5,10,15,20)	;	labelstosee <- c("0",NA,"10",NA,"20")
			coloris <- "black"		;		trancolorperc = 60		; 	trancolorperc_grey = 70
		}else{		xaxisname=NA	;	coloris <- "grey69" 	;		trancolorperc = 80		; 	trancolorperc_grey = 90
			labelslogical=F		;	yaxisticks=c(0,200,400,600,800)		;	labelstosee <- c("0",NA,"400",NA,"800")	}
		if( grepl("refpop",secdatadist) ){			rownameis="WML\nPop structure\ncorrection"	}
		if( grepl("snp",secdatadist) ){		rownameis="SNPs\nPop structure\ncorrection"			}
		plot(h, 	col=t_col("grey44",percent=trancolorperc_grey),	border=t_col("grey44",percent=trancolorperc_grey),		ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n", 	yaxt="n",	main=NA,	xlab=NA,		ylab=NA )
		lines(mydensity,col="grey60")
		plot(hs, 	col=t_col("red",percent=trancolorperc),			border=t_col("red",percent=trancolorperc),				ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n",	yaxt="n",	main=NA,	xlab=NA,		ylab=NA, 	cex.lab=1.3,	font.lab=1,	add=T)
		if ( grepl("PC5",secdatadist) ){	
			axis(4,	at=yaxisticks, labels=labelstosee)	
			}#else{	axis(4,	at=yaxisticks, labels=F, col="grey69")		}
		if ( grepl("PC5",secdatadist) & grepl("refpop",secdatadist) ) {title(ylab="Frequency", font.lab=1,cex.lab=1,outer=F, line= -12)	}
		if ( grepl("PC3",secdatadist) & grepl("relative",secdatadist) ) {mtext("Correlation Coefficient", side=1,	line=3,	font=1,cex=.7, adj=.3)	}
		if ( grepl("PC1",secdatadist) ){	mtext(rownameis, side=2, line=2, cex=.7,	las=1, 	font=2, padj=1.2) 	} # padj moves the text downward 
		axis(1,	at=c(-1,-.5,0,.5,1), labels=labelslogical, col=coloris)
		if( grepl("snp",secdatadist) ){	title(main=	sapply(strsplit(header_name,"_"),"[[",2)	, line=1,	cex.main=1)	}
		if( grepl("snp",secdatadist) & grepl("PC3",secdatadist)	){
			legend(-3.6,1400, legend=c("Expected null distribution NOT CORRECTED (Permuted environmental axes)","Real distribution CORRECTED (Real environmental axes)"),	fill=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)),
				ncol=1, bty="n", cex=1,	border=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)))	}
		### Get the other microbiomes plotted
}	}	
dev.off()
'
'
#	Test out
# alternative distribution (NOT corrected), real distribution (corrected), NO means, ALL PLOTS

# secdatadist = c("snp","refpop","smps","singleASV_absolute","singleASV_relative","genus_absolute","genus_absolute") 
pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/corcoef_distributions","_",analysis_type,"_RANDeach_NOCORRECT_ALLplots.pdf"),
			height=11,width=7)
par(mfrow=c(7,5),oma=c(4,10.5,5,4),mar=c(.8,.8,.8,.8),xpd=NA)
for( datadist in names(listofdists)[grep("RealCor_COR", names(listofdists))] ){
	message("\n",datadist)
	secdatadist_name <- sapply(strsplit(datadist,"RealCor_COR_"),"[[",2)
	names_togrepfrom <- names(listofdists)[	grep(secdatadist_name, names(listofdists))]
	header_name <- secdatadist_name# sapply(strsplit(secdatadist_name,"Cor_"),"[[",2)
	for( secdatadist in names_togrepfrom[	grep("RandCor_NO",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
		message(secdatadist)
		if( grepl("absolute",secdatadist) || grepl("relative",secdatadist) ){
			breaknumber=13	; maxisY <- 20 
		}else{
			breaknumber=150	; maxisY <- 800
		} 
		h <- hist((listofdists[[secdatadist]]),breaks=breaknumber,plot=F)
		hs <- hist((listofdists[[datadist]]),breaks=breaknumber,plot=F)		#;		maxisYy <- max(c(h$counts,hs$counts)) ; 	print(maxisYy) 
		# get density function
			multiplier <- h$counts / h$density
			mydensity <- density(listofdists[[secdatadist]])
			mydensity$y <- mydensity$y * multiplier[1]
		colnameis <- sapply(strsplit(header_name,"_"),"[[",2)	;	rownameis <- sapply(strsplit(header_name,"_"),"[[",1)
		message(rownameis) 
		#if( grepl("PC1",secdatadist) ){		yaxisname="Frequency"	}else{ 		yaxisname=NA	 }  #if first columns 
		if( grepl("absolute",secdatadist) & grepl("genus",secdatadist) ){		rownameis="Microbiome Genus\nAbs. abundance\nPop structure\ncorrection"
			xaxisname="Correlation Coefficient" ; labelslogical=T		;		yaxisticks=c(0,5,10,15,20)	;	labelstosee <- c("0",NA,"10",NA,"20")
			coloris <- "black"		;		trancolorperc = 60		; 	trancolorperc_grey = 70
		} else if ( grepl("relative",secdatadist) || grepl("absolute",secdatadist) ) {		trancolorperc = 60		; 	trancolorperc_grey = 70		;	yaxisticks=c(0,5,10,15,20) ;	labelstosee <- c("0",NA,"10",NA,"20")
		}else{		xaxisname=NA	;	coloris <- "grey69" 	;		trancolorperc = 85		; 	trancolorperc_grey = 90
			labelslogical=F		;	yaxisticks=c(0,200,400,600,800)		;	labelstosee <- c("0",NA,"400",NA,"800")	}
		if( grepl("relative",secdatadist) & grepl("genus",secdatadist) ){		rownameis="Microbiome Genus\nRel. abundance\nPop structure\ncorrection" 	}
		if( grepl("absolute",secdatadist) & grepl("ASV",secdatadist) ){		rownameis="Microbiome ASVs\nAbs. abundance\nPop structure\ncorrection" 	}
		if( grepl("relative",secdatadist) & grepl("ASV",secdatadist) ){		rownameis="Microbiome ASVs\nRel. abundance\nPop structure\ncorrection" 	}
		if( grepl("smps",secdatadist) ){		rownameis="WML\nCis Var.\ncorrection"		}
		if( grepl("refpop",secdatadist) ){			rownameis="WML\nPop structure\ncorrection"	}
		if( grepl("snp",secdatadist) ){		rownameis="SNPs\nPop structure\ncorrection"			}
		plot(h, col=t_col("grey44",percent=trancolorperc_grey),	border=t_col("grey44",percent=trancolorperc_grey),	ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n", 	yaxt="n",	main=NA,	xlab=NA,		ylab=NA)
		lines(mydensity,col="grey60")
		plot(hs, col=t_col("red",percent=trancolorperc),		border=t_col("red",percent=trancolorperc),			ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n",	yaxt="n",	main=NA,	xlab=NA,		ylab=NA,	cex.lab=1.3,	font.lab=1,	add=T)
		if ( grepl("PC5",secdatadist) ){	
			axis(4,	at=yaxisticks, labels=labelstosee)	
			}#else{	axis(4,	at=yaxisticks, labels=F, col="grey69")		}
		if ( grepl("PC5",secdatadist) & grepl("relative",secdatadist) & grepl("ASV",secdatadist) ) {title(ylab="Frequency", font.lab=1,cex.lab=1.3,outer=F, line= -11)	}
		if ( grepl("PC3",secdatadist) & grepl("absolute",secdatadist) & grepl("genus",secdatadist) ) {mtext("Correlation Coefficient", side=1,	line=3,	font=1,cex=.9, adj=.3)	}
		if ( grepl("PC1",secdatadist) ){	mtext(rownameis, side=2, line=1.5, cex=.8,	las=1, 	font=2, padj=1.2) 	} # padj moves the text downward 
		axis(1,	at=c(-1,-.5,0,.5,1), labels=labelslogical, col=coloris)
		if( grepl("snp",secdatadist) ){	title(main=	sapply(strsplit(header_name,"_"),"[[",2)	, line=1,	cex.main=1.3)	}
		if( grepl("snp",secdatadist) & grepl("PC3",secdatadist)	){
			legend(-5.6,1410, legend=c("Expected null distribution NOT CORRECTED (Permuted environmental axes)","Real distribution CORRECTED (Real environmental axes)"),	fill=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)),
				ncol=1, bty="n", cex=1.3,	border=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)))	
		}
}	}
dev.off()
'


# alternative distribution (corrected), real distribution (corrected), NO means, ALL PLOTS
if (fileout=="pdf"){
	pdf(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/","Figure_s11_","corcoef_distributions","_",analysis_type,"_RANDeach_ALLplots.pdf"),
			height=10,width=7)
}
if (fileout=="png"){
	png(file=paste0("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/figures_for_paper/pval_distributions/","Figure_s11_","corcoef_distributions","_",analysis_type,"_RANDeach_ALLplots.png"),	height=10,width=7,units="in",res=600)
}
par(mfrow=c(7,5),oma=c(4,10.5,5,4),mar=c(.8,.8,.8,.8),xpd=NA)
for( datadist in names(listofdists)[grep("RealCor_COR", names(listofdists))] ){
	message("\n",datadist)
	secdatadist_name <- sapply(strsplit(datadist,"RealCor_COR_"),"[[",2)
	names_togrepfrom <- names(listofdists)[	grep(secdatadist_name, names(listofdists))]
	header_name <- secdatadist_name# sapply(strsplit(secdatadist_name,"Cor_"),"[[",2)
	for( secdatadist in names_togrepfrom[	grep("RandCor_COR",names(listofdists)[grep(secdatadist_name, names(listofdists))])	]	) {
		message(secdatadist)
		if( grepl("absolute",secdatadist) || grepl("relative",secdatadist) ){
			breaknumber=13	; maxisY <- 20 
		}else{
			breaknumber=150	; maxisY <- 800
		} 
		h <- hist((listofdists[[secdatadist]]),breaks=breaknumber,plot=F)
		hs <- hist((listofdists[[datadist]]),breaks=breaknumber,plot=F)		#;		maxisYy <- max(c(h$counts,hs$counts)) ; 	print(maxisYy) 
		# get density function
			multiplier <- h$counts / h$density
			mydensity <- density(listofdists[[secdatadist]])
			mydensity$y <- mydensity$y * multiplier[1]
		colnameis <- sapply(strsplit(header_name,"_"),"[[",2)	;	rownameis <- sapply(strsplit(header_name,"_"),"[[",1)
		message(rownameis) 
		#if( grepl("PC1",secdatadist) ){		yaxisname="Frequency"	}else{ 		yaxisname=NA	 }  #if first columns 
		if( grepl("absolute",secdatadist) & grepl("genus",secdatadist) ){		rownameis="Microbiome Genus\nAbs. abundance\nPop structure\ncorrection"
			xaxisname="Correlation Coefficient" ; labelslogical=T		;		yaxisticks=c(0,5,10,15,20)	;	labelstosee <- c("0",NA,"10",NA,"20")
			coloris <- "black"		;		trancolorperc = 60		; 	trancolorperc_grey = 70
		} else if ( grepl("relative",secdatadist) || grepl("absolute",secdatadist) ) {		trancolorperc = 60		; 	trancolorperc_grey = 70		;	yaxisticks=c(0,5,10,15,20) ;	labelstosee <- c("0",NA,"10",NA,"20")
		}else{		xaxisname=NA	;	coloris <- "grey69" 	;		trancolorperc = 85		; 	trancolorperc_grey = 90
			labelslogical=F		;	yaxisticks=c(0,200,400,600,800)		;	labelstosee <- c("0",NA,"400",NA,"800")	}
		if( grepl("relative",secdatadist) & grepl("genus",secdatadist) ){		rownameis="Microbiome Genus\nRel. abundance\nPop structure\ncorrection" 	}
		if( grepl("absolute",secdatadist) & grepl("ASV",secdatadist) ){		rownameis="Microbiome ASVs\nAbs. abundance\nPop structure\ncorrection" 	}
		if( grepl("relative",secdatadist) & grepl("ASV",secdatadist) ){		rownameis="Microbiome ASVs\nRel. abundance\nPop structure\ncorrection" 	}
		if( grepl("smps",secdatadist) ){		rownameis="WML\nCis Var.\ncorrection"		}
		if( grepl("refpop",secdatadist) ){			rownameis="WML\nPop structure\ncorrection"	}
		if( grepl("snp",secdatadist) ){		rownameis="SNPs\nPop structure\ncorrection"			}
		plot(h, col=t_col("grey44",percent=trancolorperc_grey),	border=t_col("grey44",percent=trancolorperc_grey),	ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n", 	yaxt="n",	main=NA,	xlab=NA,		ylab=NA)
		lines(mydensity,col="grey60")
		plot(hs, col=t_col("red",percent=trancolorperc),		border=t_col("red",percent=trancolorperc),			ylim=c(0,maxisY),	xlim=c(-1,1),	xaxt="n",	yaxt="n",	main=NA,	xlab=NA,		ylab=NA,	cex.lab=1.3,	font.lab=1,	add=T)
		if ( grepl("PC5",secdatadist) ){	
			axis(4,	at=yaxisticks, labels=labelstosee)	
			}#else{	axis(4,	at=yaxisticks, labels=F, col="grey69")		}
		if ( grepl("PC5",secdatadist) & grepl("relative",secdatadist) & grepl("ASV",secdatadist) ) {title(ylab="Frequency", font.lab=1,cex.lab=1.3,outer=F, line= -10)	}
		if ( grepl("PC3",secdatadist) & grepl("absolute",secdatadist) & grepl("genus",secdatadist) ) {mtext("Correlation Coefficient", side=1,	line=3,	font=1,cex=.9, adj=.3)	}
		if ( grepl("PC1",secdatadist) ){	mtext(rownameis, side=2, line=1.5, cex=.8,	las=1, 	font=2, padj=1.2) 	} # padj moves the text downward 
		axis(1,	at=c(-1,-.5,0,.5,1), labels=labelslogical, col=coloris)
		if( grepl("snp",secdatadist) ){	title(main=	sapply(strsplit(header_name,"_"),"[[",2)	, line=1,	cex.main=1.3)	}
		if( grepl("snp",secdatadist) & grepl("PC3",secdatadist)	){
			legend(-5.6,1410, legend=c("Expected null distribution (Permuted environmental axes)","Real distribution (Real environmental axes)"),	fill=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)),
				ncol=1, bty="n", cex=1.3,	border=c(t_col("grey44",percent=trancolorperc_grey-15),t_col("red",percent=trancolorperc-15)))	
				# both were CORRECTED for population structure
		}
}	}
dev.off()














# Make venn/euler diagrams

#Prepare data: 	change names to be identical.
	# if smps genes names has a "_no" appended, it means it was run using simple mantel, since there were no snp variation to correct for. Removed in list_name_overlap_even

list_name_overlap_even <- list()
for (	i in 1:NROW(list_name_overlap) ){
	list_name_overlap_even[[i]] 	<- 	gsub("gene_snp_|gene_smp_|relative_|absolute_|_no","",x=list_name_overlap[[i]])
}
names(list_name_overlap_even) <- names(list_name_overlap)

list_name_overlap_even_coef_pval <- list()
for (	i in 1:NROW(list_name_overlap_coef_pval) ){
	list_name_overlap_even_coef_pval[[i]] 	<- 	gsub("gene_snp_|gene_smp_|relative_|absolute_|_no","",x=list_name_overlap_coef_pval[[i]])
}
names(list_name_overlap_even_coef_pval) <- names(list_name_overlap_coef_pval)


if( COEF_ONLY==T ) {	dataset_touse <- list_name_overlap_even 	}
#if( SIGN_ONLY==T ) {	message("ERROR: Venn diagrams: no dataset for significance only (may exists, but not registered in code)")}
#if( COEF_AND_SIG==T ) {		dataset_touse <- list_name_overlap_even_coef_pval 	}	# list_name_overlap_even (one of them is only coef, the other pval and coef.


# make overlap table:
for ( axisnumberis in 1:5 ){
	list_elm_togrep 	<- grep( paste0("PC",axisnumberis),	names(dataset_touse)	)
	greplist <- names(dataset_touse)[list_elm_togrep]
	
	if( analysis_type=="mrm" ){
		snp <- list_elm_togrep[	grep("snp", greplist)	]		;		smp <- list_elm_togrep[	grep("smp", greplist)	]
		taxa_rel <-  list_elm_togrep[	grep("genuslevel_relative_", greplist)	]		;		taxa_abs <-  list_elm_togrep[	grep("genuslevel_absolute_", greplist)	]
		ASV_rel <- list_elm_togrep[	grep("singleasvs_relative_", greplist)	]			;		ASV_abs <-  list_elm_togrep[	grep("singleasvs_absolute_", greplist)	]
		
		snpN <-	NROW(dataset_touse[[snp]])		;	smpN <- NROW(dataset_touse[[smp]])	
		taxa_relN <- 	NROW(dataset_touse[[taxa_rel]])		;		taxa_absN <-  	NROW(dataset_touse[[taxa_abs]])
		ASV_relN <- 	NROW(dataset_touse[[ASV_rel]])			;		ASV_absN <-  	NROW(dataset_touse[[ASV_abs]])
		overlap_snp_smp			<- NROW(intersect(dataset_touse[[snp]],  dataset_touse[[smp]]))
		overlap_taxa_abs_rel	<- NROW(intersect(dataset_touse[[taxa_rel]],  dataset_touse[[taxa_abs]]))
		overlap_ASV_abs_rel		<- NROW(intersect(dataset_touse[[ASV_rel]],  dataset_touse[[ASV_abs]]))
		
		if ( paste0("PC",axisnumberis) == "PC1"	){
			ovelap_dataframe 	<- data.frame(	 paste0("PC",axisnumberis),	snpN,		smpN,			overlap_snp_smp,
			ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)
		}else{
			ovelap_dataframe	<- rbind(ovelap_dataframe, data.frame(	 paste0("PC",axisnumberis),	snpN,		smpN,			overlap_snp_smp,		
																ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)	)
		}
	}else{
		snp <- list_elm_togrep[	grep("snp", greplist)	]		;		refpop <- list_elm_togrep[	grep("refpop", greplist)	]		;			smp <- list_elm_togrep[	grep("smp", greplist)	]
		taxa_rel <-  list_elm_togrep[	grep("genuslevel_relative_", greplist)	]		;		taxa_abs <-  list_elm_togrep[	grep("genuslevel_absolute_", greplist)	]
		ASV_rel <- list_elm_togrep[	grep("singleASV_relative_", greplist)	]			;		ASV_abs <-  list_elm_togrep[	grep("singleASV_absolute_", greplist)	]
		
		snpN <-	NROW(dataset_touse[[snp]])		;	refpopN <- NROW(dataset_touse[[refpop]])		;		smpN <- NROW(dataset_touse[[smp]])	
		taxa_relN <- 	NROW(dataset_touse[[taxa_rel]])		;		taxa_absN <-  	NROW(dataset_touse[[taxa_abs]])
		ASV_relN <- 	NROW(dataset_touse[[ASV_rel]])			;		ASV_absN <-  	NROW(dataset_touse[[ASV_abs]])
		ovelap_smp_refpop 		<- NROW(intersect(dataset_touse[[smp]],  dataset_touse[[refpop]]))
		overlap_snp_smp			<- NROW(intersect(dataset_touse[[snp]],  dataset_touse[[smp]]))
		overlap_snp_refpop		<- NROW(intersect(dataset_touse[[snp]],  dataset_touse[[refpop]]))
		overlap_taxa_abs_rel	<- NROW(intersect(dataset_touse[[taxa_rel]],  dataset_touse[[taxa_abs]]))
		overlap_ASV_abs_rel		<- NROW(intersect(dataset_touse[[ASV_rel]],  dataset_touse[[ASV_abs]]))
		
		if ( paste0("PC",axisnumberis) == "PC1"	){
			ovelap_dataframe 	<- data.frame(	 paste0("PC",axisnumberis),	snpN,		refpopN,					smpN,			overlap_snp_refpop,		overlap_snp_smp,		ovelap_smp_refpop,
										ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)
		}else{
			ovelap_dataframe	<- rbind(ovelap_dataframe, data.frame(	 paste0("PC",axisnumberis),	snpN,		refpopN,					smpN,			overlap_snp_refpop,		overlap_snp_smp,		ovelap_smp_refpop,
																ASV_relN,	ASV_absN,	overlap_ASV_abs_rel,	taxa_relN,		taxa_absN,				overlap_taxa_abs_rel)	)
		}
	}
}
if (analysis_type=="mrm") {
	colnames(ovelap_dataframe) <- c("env_parameter",	"snp_N",	"smp_cisref_N",		"overlap_snp_cisref_N",		
								"ASV_rel_N",	"ASV_abs_N",	"ASV_overlap_abs_rel_N",	"Genera_rel_N",	"Genera_abs_N",	"Genera_overlap_rel_abs_N")
}else{
	colnames(ovelap_dataframe) <- c("env_parameter",	"snp_N",	"smp_refpop_N",		"smp_cisref_N",		"overlap_snp_refpop_N",		"overlap_snp_cisref_N",		"overlap_smp_refpop_cisref_N",
								"ASV_rel_N",	"ASV_abs_N",	"ASV_overlap_abs_rel_N",	"Genera_rel_N",	"Genera_abs_N",	"Genera_overlap_rel_abs_N")
}
write.csv(file= paste0("./figures_for_paper/COR_Coef_overlapDat_PCaxes_",analysis_type,dataset_is,".csv"),ovelap_dataframe)


# make white space figure:	for padding/space between plots. Found a better way.
	#	pdf(file="figures_for_paper/sum_up_plots/white_space.pdf", bg="white")	;plot.new()	;dev.off()

# PLOT euler venn 

colorsare <- c(colorsare1[1:3],colorsare1[4],colorsare1[4],colorsare1[5],colorsare1[5])
widthis=7
#emptywidth= (widthis/5)* 0.09  	# making room for empty plots when combining the panel (they will be 10% of one of the five plots. Adjust in this line.
#newwidth= (widthis/5)- 2*(emptywidth) #redundant method now.
newwidth=widthis/5
margin_padding=.0
set.seed(10)
list_plots <- list()
for ( axisnumberis in 1:5 ){
		if (analysis_type=="mrm") {
			euler_vector <- c(	"A"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "snp_N") ]	,
								"B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_cisref_N") ],
								"A&B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_cisref_N") ]
								)
			colorvec_is=colorsare[1:2]
		}else{
			euler_vector <- c(	"A"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "snp_N") ]	,
								"B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_refpop_N") ],
								"C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "smp_cisref_N") ],
								"A&B"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_refpop_N") ],
								"A&C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_snp_cisref_N") ],
								"B&C"= ovelap_dataframe[axisnumberis, which(colnames(ovelap_dataframe) == "overlap_smp_refpop_cisref_N") ],
								"A&B&C"= NROW(intersect( intersect(dataset_touse[[snp]],  dataset_touse[[refpop]]) , dataset_touse[[smp]]))
								)
			colorvec_is=colorsare[1:3]
		}
		euler_diag <- euler(euler_vector, input="union")
		
		plottits1 <- plot(	euler_diag,	quantities = list(alpha=1, cex=.7), 				fill=colorvec_is,
				labels=NA,				alpha=.6	,adjust_labels=T ,	edges=list(alpha=.01, col="white")	)	#	main=paste0("PC ",axisnumberis),
		plottits <- arrangeGrob(plottits1, ncol=1,			padding=unit(margin_padding,"in"),top="",bottom="",right="",left="")
		list_plots	<- append(list_plots, list(plottits))
		
	
	if( COEF_ONLY==T ){
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_diagramm_jan22_coef_PC",axisnumberis,".pdf"), 
			plot=plottits,	height=heights_ofplots[3],			width=newwidth,units="in")
	#	ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_diagramm_jan22_coef_PC",axisnumberis,".png"), 
	#		plot=plottits, 	height=heights_ofplots[3],			width=newwidth,units="in",dpi=800)
	}
}



### combine panels

if (analysis_type=="mantel"){
			#combine sum-up plot and venn-plot
	panel_top <- rasterGrob(image_read_pdf(paste0("figures_for_paper/pval_distributions/Figure_5a_corcoef_distributions_mantel_RANDeach_reduced.pdf")),interpolate=T)				# read in plots from picture file
			# heigh 5
	#empty_panel <- rasterGrob(image_read_pdf("figures_for_paper/sum_up_plots/white_space.pdf")		,interpolate=T)
	if( COEF_ONLY==T ) {	
		panel_mid <- rasterGrob(image_read_pdf(paste0("figures_for_paper/sum_up_plots/Figure_5bcd_COR_Coef_plots_sunimp2_mantel_percentile_0.9999.pdf")),interpolate=T)				# read in plots from picture file
			# height 3.5
		}
		
	panel_bottom	<- 	arrangeGrob(	list_plots[[1]],	list_plots[[2]],	list_plots[[3]],	list_plots[[4]],	list_plots[[5]],
												ncol=5,	heights=unit(heights_ofplots[3], c("in")), widths=unit(rep(widthis/5,5), c("in")),
												padding=unit(0,"line"),top="",bottom="",right="",left="")
	combined_panel 	<- arrangeGrob(panel_top,	panel_mid,		panel_bottom, nrow=3,	heights=unit(heights_ofplots, c("in"))	)	

	dev.new(height=sum(heights_ofplots), width=widthis, units = "in", res = 1000)	
	height_annotate_middle =	((heights_ofplots[2]+heights_ofplots[3])/sum(heights_ofplots))+.02
	height_annotate_low =	heights_ofplots[3]/sum(heights_ofplots)+.01
	ggdraw() + 
		draw_plot(combined_panel) + 
		draw_plot_label(c("a", "b", "c", "d", "e"), c(0, 0, widthis/3/widthis,	widthis/3*2/widthis,	0), c(1, height_annotate_middle, height_annotate_middle, height_annotate_middle,	height_annotate_low), size = 15) +
		draw_plot_label(c("PC1", "PC2", "PC3", "PC4", "PC5"), c(seq(from=0,to=widthis,length.out=100))[c(10,29,49,69,89)]/widthis, .03, size = 8, hjust = 0, vjust=2)

	if( COEF_ONLY==T ) {	
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/Figure_5_",analysis_type,"_dist_PCA_venn_coefOnly_EULER.pdf"), height=sum(heights_ofplots) ,width=widthis ,units="in")
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/Figure_5_",analysis_type,"_dist_PCA_venn_coefOnly_EULER.png"), height=sum(heights_ofplots) ,width=widthis ,units="in",dpi=600)
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/Figure_5.pdf"), height=sum(heights_ofplots) ,width=widthis ,units="in")
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/Figure_5.png"), height=sum(heights_ofplots) ,width=widthis ,units="in",dpi=800)
		}
	dev.off()
}else{
	panel_bottom	<- 	arrangeGrob(	list_plots[[1]],	list_plots[[2]],	list_plots[[3]],	list_plots[[4]],	list_plots[[5]],
											ncol=5,	heights=unit(heights_ofplots[3], c("in")), widths=unit(rep(widthis/5,5), c("in")),
											padding=unit(0,"line"),top="",bottom="",right="",left="")
											
	ggdraw() + 
		draw_plot(panel_bottom) + 
		draw_plot_label(c("PC1", "PC2", "PC3", "PC4", "PC5"), c(seq(from=0,to=widthis,length.out=100))[c(10,29,49,69,89)]/widthis, .2, size = 8, hjust = 0, vjust=2)
		
	if( COEF_ONLY==T ) {	
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_coefOnly_EULER.pdf"), height=heights_ofplots[3] ,width=widthis ,units="in")
		ggsave(filename=paste0("figures_for_paper/sum_up_plots/",analysis_type,"_PCA_venn_coefOnly_EULER.png"), height=heights_ofplots[3] ,width=widthis ,units="in",dpi=600)
		}
}
