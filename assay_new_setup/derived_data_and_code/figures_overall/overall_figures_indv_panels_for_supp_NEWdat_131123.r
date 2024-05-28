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

#output_type="pdf" # "png" # "pdf"
output_type="png"

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


####################################
#------------- Growth--------------#
####################################

if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_Suppl_pheno1.pdf"),width=7,height=4) 	
}else if( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_Suppl_pheno1.png"),width=7,height=4,units="in",res=400) 	
}
par( mfrow=c(1,2), mar=c(1,3,1,1), oma=c(3.5,.5,1.5,0), cex=.9)		
par()$cex

Grow <- read.table("assay_new_setup/derived_data_and_code/phenotypes/data/growth_surv/growth_and_survival.txt", header = TRUE, sep = "", dec = '.',stringsAsFactors=F,as.is=c(F,F,F,T,T,	T,T,T,T,T,	T))

xlab_name=expression(paste("Acclimation Temperature (", degree,"C)"))
ylimsis = c(min(Grow$growth,na.rm=T),	max(Grow$growth,na.rm=T) + (max(Grow$growth,na.rm=T)-min(Grow$growth,na.rm=T))	*proportion_to_increase_ylim	)	#c(-.2,4.2)
jitter_move=seq(-.5,.5,length.out=4)
for ( popul in levels(Grow$Population) ){
	counter=grep(popul,levels(Grow$Population))
	#print(counter) ;	print(popul)
	plot_rows <- which(Grow$Population == popul)
	subset_dat <- Grow[plot_rows,]
	if ( popul == levels(Grow$Population)[1]	){		
		par(new=F)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$growth, xpd=NA,
			col=popcols[counter], ylab=NA, xlab=xlab_name, ylim=ylimsis, xlim=c(14.5,29.5), xaxt="n", las=2)
		title(ylab="Growth factor", line=2, xpd=NA)
		legend("topright",legend=levels(Grow$Population),col=popcols, horiz=T, xpd=NA, pch=1, bty="n", inset=c(-1.05,-.17), text.width=5.5)
		axis(1, at=c(15,17,19,21,23,25,27,29), labels=NA, xpd=NA, col.ticks="grey80", col="white", tcl=-.35)
		axis(1, at=c(15,19,23,25,29), labels=c(15,19,23,25,29), xpd=NA, col.ticks="black")
		print(par()$cex)
	}else{	
		par(new=T)
		# message("entering if else")
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$growth, 	xpd=NA,
			col=popcols[counter], ylab=NA, xlab=NA, ylim=ylimsis, xlim=c(14.5,29.5),yaxt="n", xaxt="n")
		print(par()$cex)
	}
	loess_fit <- loess( subset_dat$growth ~ subset_dat$Treatment)
	i_ord <- order(subset_dat$Treatment)
	abline(lm(subset_dat$growth ~ subset_dat$Treatment),col=popcols[counter])
} 
mtext("a",side=3,line=0, at=13,cex=1.2)
print(par()$cex)
####################################
#------------- Survival------------#
####################################
surv_deaths_dat <- cbind(Grow$spiders_3,Grow$spiders_1-Grow$spiders_3)
colnames(surv_deaths_dat) <- c("N_alive","N_dead")
surv_model_T <- glm(surv_deaths_dat ~ as.factor(Grow$Treatment) + 0,  family=binomial(link="logit"))
#summary(surv_model_T)

xlab_name=expression(paste("Acclimation Temperature (", degree,"C)"))
ylimsis=c(min(Grow$survival,na.rm=T),	max(Grow$survival,na.rm=T) + (max(Grow$survival,na.rm=T)-min(Grow$survival,na.rm=T))*proportion_to_increase_ylim)
for ( popul in levels(Grow$Population) ){
	counter=grep(popul,levels(Grow$Population))
	plot_rows <- which(Grow$Population == popul)
	subset_dat <- Grow[plot_rows,]
	if ( popul == levels(Grow$Population)[1]	) {		
		par( new=F)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$survival, xpd=NA,
			col=popcols[counter], ylab=NA, xlab=xlab_name, ylim=ylimsis,xlim=c(14.5,29.5),xaxt="n", las=2)
		title(ylab="Survival", line=2.5, xpd=NA)
		axis(1, at=c(15,17,19,21,23,25,27,29), labels=NA, xpd=NA, col.ticks="grey80", col="white", tcl=-.35)
		axis(1, at=c(15,19,23,25,29), labels=c(15,19,23,25,29), xpd=NA, col.ticks="black")
	}else{	
		par(new=T)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$survival, xpd=NA,
			col=popcols[counter], ylab=NA, xlab=NA, ylim=ylimsis, xlim=c(14.5,29.5),yaxt="n", xaxt="n", las=2)
	}
	loess_fit <- loess( subset_dat$survival ~ subset_dat$Treatment)
	i_ord <- order(subset_dat$Treatment)
	lines(subset_dat$Treatment[i_ord], predict(loess_fit)[i_ord], col =	popcols[counter])
} 
mtext("b",side=3,line=0, at=13,cex=1.2)

dev.off()


####################################
#------------- CTmax---------------#
####################################

if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_2_phenoTemp1.pdf"),width=7,height=4) 		
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_2_phenoTemp1.png"),width=7,height=4,units="in",res=400) 	
}
par( mfrow=c(1,2),mar=c(1,3,1,1),oma=c(3.5,.5,1.5,0))

ToD_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/ToD/ToD_nest-temp_mean.csv", header=T,stringsAsFactors=T) 

mainname="CTmax"; xlab_name=expression(paste("Acclimation Temperature (", degree,"C)"))
ylimsis <- c(min(ToD_final$ToD_mean,na.rm=T), max(ToD_final$ToD_mean,na.rm=T) + (max(ToD_final$ToD_mean,na.rm=T)-min(ToD_final$ToD_mean,na.rm=T))*proportion_to_increase_ylim)
jitter_move=seq(-.5,.5,length.out=4)
for ( popul in levels(ToD_final$Population) ){
	counter=grep(popul,levels(ToD_final$Population))
	print(popul)
	plot_rows <- which(ToD_final$Population == popul)
	subset_dat <- ToD_final[plot_rows,]
	if ( popul == levels(ToD_final$Population)[1]	) {		
		par( new=F, cex=.9)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$ToD_mean, xpd=NA, 	col=popcols[counter], ylab=NA,xlab=xlab_name,		ylim=ylimsis,	xlim=c(14.5,29.5),xaxt="n",  las=2)#, ylim=c(44,54))
		title(ylab=expression(paste("CTmax (",degree,"C)")), line=2, xpd=NA)
		legend("topright",legend=levels(ToD_final$Population),col=popcols, horiz=T, xpd=NA, pch=1, bty="n", inset=c(-1.05,-.17), text.width=5.5)
		axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)
	}else{	
		par( new=T, cex=.9)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$ToD_mean, 	xpd=NA,	col=popcols[counter], ylab=NA, xlab=NA, ylim=ylimsis, xlim=c(14.5,29.5),yaxt="n", xaxt="n")#, ylim=c(44,54))
	}
	loess_fit <- loess( subset_dat$ToD_mean ~ subset_dat$Treatment)
	i_ord <- order(subset_dat$Treatment)
	abline(lm(subset_dat$ToD_mean ~ subset_dat$Treatment),col=popcols[counter])
	print(lm(subset_dat$ToD_mean ~ subset_dat$Treatment))
} 
print(summary(lm(ToD_mean ~  Treatment * Population, data=ToD_final)))
mtext("a",side=3,line=0, at=13,cex=1.2)
'Effect size (slope):

[1] "Betta"
Call:
lm(formula = subset_dat$ToD_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
            48.34632               0.09675

[1] "Karasburg"
Call:
lm(formula = subset_dat$ToD_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
            49.94567               0.01103

[1] "Otavi"
Call:
lm(formula = subset_dat$ToD_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
            48.85100               0.07302

[1] "Stampriet"
Call:
lm(formula = subset_dat$ToD_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
            49.51064               0.04184

'

####################################
#------------- CCR-----------------#
####################################
CCR_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/CCR/CCR_nest-temp_mean.csv", header=T,stringsAsFactors=T) 

# Effect of temperature
	#ctmax
#par(mfcol=c(2,2))

mainname="CCR"	; xlab_name=expression(paste("Acclimation Temperature (", degree,"C)"))
ylimsis <- c(min(CCR_final$CCR_mean,na.rm=T), max(CCR_final$CCR_mean,na.rm=T) + (max(CCR_final$CCR_mean,na.rm=T)-min(CCR_final$CCR_mean,na.rm=T))*proportion_to_increase_ylim)
jitter_move=seq(-.5,.5,length.out=4)
for ( popul in levels(CCR_final$Population) ){
	counter=grep(popul,levels(CCR_final$Population))
	print(popul)
	plot_rows <- which(CCR_final$Population == popul)
	subset_dat <- CCR_final[plot_rows,]
	if ( popul == levels(CCR_final$Population)[1]	) {		
		par( new=F)
		par( cex=.9)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$CCR_mean, xpd=NA, 	col=popcols[counter], ylab=NA,xlab=xlab_name,		ylim=ylimsis,	xlim=c(14.5,29.5),xaxt="n", las=2)#, ylim=c(44,54))
		title(ylab=expression(paste("CCR (",degree,"C)")), line=2, xpd=NA)
		axis(1, at=c(15,17,19,21,23,25,27,29), labels=c(15,NA,19,NA,23,25,NA,29), xpd=NA)
	}else{	
		par( new=T, cex=.9)
		plot(jitter(subset_dat$Treatment+jitter_move[counter], factor=.15), subset_dat$CCR_mean, 	xpd=NA,	col=popcols[counter], ylab=NA, xlab=NA, ylim=ylimsis, xlim=c(14.5,29.5),yaxt="n", xaxt="n")#, ylim=c(44,54))
	}
	loess_fit <- loess( subset_dat$CCR_mean ~ subset_dat$Treatment)
	i_ord <- order(subset_dat$Treatment)
	abline(lm(subset_dat$CCR_mean ~ subset_dat$Treatment),col=popcols[counter])
	print(lm(subset_dat$CCR_mean ~ subset_dat$Treatment))
} 
mtext("b",side=3,line=0, at=13,cex=1.2)
dev.off()
' Effect sizes:
[1] "Betta"
Call:
lm(formula = subset_dat$CCR_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
           11.658871             -0.003623

[1] "Karasburg"
Call:
lm(formula = subset_dat$CCR_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
              7.3243                0.1668

[1] "Otavi"
Call:
lm(formula = subset_dat$CCR_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
           11.945270              0.003831

[1] "Stampriet"
Call:
lm(formula = subset_dat$CCR_mean ~ subset_dat$Treatment)

Coefficients:
         (Intercept)  subset_dat$Treatment
              8.8665                0.1139
'



####################################
#------------- Methylations	-------#
####################################

fdr_threshold = 0.1
fdr_threshold = 0.05
	# number of genes showing patterns
meth_file = "assay_new_setup/derived_data_and_code/methylation/smps/gene_region/DSS_output/DSS_result_all_d10D32_gene_CPG.txt"
meth_datis <- read.table(meth_file,	header=T,	stringsAsFactors=T)
meth_dat_pop_sign <- meth_datis[which(meth_datis$fdr_p < fdr_threshold), ]
meth_dat_temp_sign <- meth_datis[which(meth_datis$fdr_T < fdr_threshold), ]
meth_dat <- c(nrow(meth_dat_pop_sign),nrow(meth_dat_temp_sign),NROW(intersect( droplevels(meth_dat_pop_sign$region_name), droplevels(meth_dat_temp_sign$region_name) )),NA)
names(meth_dat) <- c("Population","Temperature","Pop+Temp","Interaction")
#write.table(file = "/faststorage/project/spider2/assay_study_anne/methylation/smps/gene_region//DSS_result_counteffects_d10D32_gene.txt", t(meth_dat))
meth_dat
' Population Temperature    Pop+Temp Interaction
       1537           0           0          NA
0.05	
 Population Temperature    Pop+Temp Interaction 
       1136           0           0          NA 
'
	
	#for expression data (also fdr - why two different cutoffs) earlier: 0.05
padj_threshold = 0.1
padj_threshold = fdr_threshold
	# import methylation data
	# import expression data:
exp_file = "assay_new_setup/derived_data_and_code/expression/results/likelyhood_ratio_tests/ptpt_effects_table.txt"
exp_dat <- read.table(exp_file, header=T,       row.names=1, stringsAsFactors=T)
exp_dat_pop <-	exp_dat[	which( exp_dat[,which( colnames(exp_dat) == "padj_P" )] < padj_threshold	)	,	]
exp_dat_Temp <- exp_dat[	which( exp_dat[,which( colnames(exp_dat) == "padj_T" )] < padj_threshold	)	,	]

overlap_genes_exprP_methP <- NROW(intersect(exp_dat_pop[,1],meth_dat_pop_sign[,1]))
message("Ekspression and methylation: Number of genes overlapping between EXSP-population and METH-population related genes: ", overlap_genes_exprP_methP)
'
Ekspression and methylation: Number of genes overlapping between EXSP-population and METH-population related genes: 844
0.05
Ekspression and methylation: Number of genes overlapping between EXSP-population and METH-population related genes: 572
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
'
[1] 13998 11807  7669   909
'
###############################################
#------------- Expression / RNA plotting  ----#
###############################################
	
	
	# number of genes showing patterns
if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_3_Meth_Exp1_",fdr_threshold,".pdf"),width=7,height=3) 		
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_3_Meth_Exp1_",fdr_threshold,".png"),width=7,height=3,units="in",res=400) 	
}
par( mfrow=c(1,2), mar=c(1,1.5,0,1), oma=c(3,3.5,4.3,1)	, cex=.8, cex.axis=.9)

barppp <- barplot(rev(datis),        beside=T, horiz=T, col=opaque_col1(col_vec[2]), 
		xlab=paste0("Number of genes (N: ",nrow(exp_dat),")"), xpd=NA, axes=F)
title(main="Expression", xpd=NA, line=2.9)
barplot(rev(c(overlap_genes_exprP_methP,0,0,0)),        beside=T, horiz=T, col="black", names.arg=NA,
		density=c(rep(17,3)), angle=c(rep(135,3)), xpd=NA, add=T, xaxt="n")
axis(1, at=c(0,seq(2000,12000,by=2000)), labels=c(0,NA,4000,NA,8000,NA,12000),xpd=NA)
legend("top",legend="N Genes also showing\npopulation effect in methylation", fill=opaque_col1(col_vec[2])[4], 
		bty="n", cex=.9, inset=c(0,-.25), xpd=NA)
legend("top",legend="N Genes also showing\npopulation effect in methylation", density=17,angle=135, fill="black", 
		bty="n", cex=.9, inset=c(0,-.25), xpd=NA)		
text(	x=-max(datis,na.rm=T)*.09,
		y = barppp,
		labels=colnames(custom_dat), xpd=NA,
		srt=45, adj=c(.9,1.5))
mtext("a",side=3,line=0.5, at=max(barppp[4,]) * .05,cex=1.2)

####################################
#	Methylation plotting 	 ###
####################################
# number of genes showing patterns

plotaa <- barplot(rev(meth_dat),        beside=T, horiz=T, col=opaque_col1(col_vec[1]), 
		xlab=paste0("Number of genes (N: ",nrow(meth_datis),")"), names.arg=NA,
		ylab= NA, xpd=NA,  xlim=c(0,1200), xaxt="n")
title(main="Methylation", xpd=NA, line=2.9)
#axis(side=1, at=seq(0,1600,by=200),	labels=c(0,NA,400,NA,800,NA,1200,NA,1600))
text(	x=100,	y = plotaa[1,],	labels="NA", , xpd=NA)
# add overlap between ekpression P and methylation P
barplot(rev(c(overlap_genes_exprP_methP,0,0,0)),        beside=T, horiz=T, col="black", names.arg=NA,
		density=c(rep(17,3)), angle=c(rep(135,3)), xpd=NA, add=T, xlim=c(0,1200), xaxt="n")
axis(1, at=c(0,400,800,1200), labels=c(0,400,800,1200))
legend("top",legend="N Genes also showing\npopulation effect in expression", fill=opaque_col1(col_vec[1])[4], 
		bty="n", cex=.9, inset=c(0,-.25), xpd=NA)
legend("top",legend="N Genes also showing\npopulation effect in expression", density=17,angle=135, fill="black", 
		bty="n", cex=.9, inset=c(0,-.25), xpd=NA)
		# Add arrow indicating the number of genes showing population response (in expression) that by chance, would also show population response in methylation.
arrow_x <- 	meth_dat[1] / nrow(meth_datis) * 	datis[1]		# ( N diff METH genes / N All genes meth ) * N diff exp G Pop
arrows(x0=arrow_x,		y0=plotaa[3]+(plotaa[4]-plotaa[3])/3.35,		# increase division factor to move arrow down
		x1=arrow_x, 	y1=plotaa[3]+(plotaa[4]-plotaa[3])/1.95,		lwd=1.5, length=.045, angle=45, col="gray39")
'text(	x=-max(meth_dat,na.rm=T)*.09,
		y = plotaa,
		labels=colnames(custom_dat), xpd=NA,
		srt=45, adj=c(.9,1.5))'

mtext("b",side=3,line=0.5, at=max(plotaa[4,]) * .05,cex=1.2)
dev.off()


	# overlap diagram:
euler_vector <- c(	"A"= meth_dat[[1]]	,	# METH POP
					"B"= datis[1]	,	# EXP POP
					"C"= datis[2]	,	# EXP TEMP
					"A&B"= overlap_genes_exprP_methP,
					"A&C"= overlap_genes_exprT_methP,
					"B&C"= overlap_genes_pop_temp,
					"A&B&C"= overlap_genes_exprT_expP_methP					)	
euler_diag <- euler(euler_vector, input="union")
eul_plot <- plot(	euler_diag,	quantities = list(alpha=1, cex=.7), fill=c(	opaque_col1(col_vec[1])[2],	opaque_col1(col_vec[2])[3],	opaque_col1(col_vec[2])[2]),
				labels=list(labels=c("Methylation\nPopulation effect", "Expression\nPopulation effect", "Expression\nTemperature effect")),	cex=.7,			
				alpha=.8, adjust_labels=T,	edges=list(alpha=.01, col="white")	)	#	main=paste0("PC ",axisnumberis),
#center_dat <- eul_plot$data$centers
ggsave(plot=eul_plot,	filename=paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"_EULER_exp_meth",fdr_threshold,".pdf"), 
		height=6 ,width=6 ,units="in")
ggsave(plot=eul_plot,	filename=paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"_EULER_exp_meth",fdr_threshold,".png"), 
		height=6 ,width=6 ,units="in",dpi=400)

####################################
#------------- Microbiome 1 plot---#
####################################
	# number of ASVs showing patterns (after cutoff of 25% prevalence in at least one pop/temp)

if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_4_Micro_ASV.pdf"),width=4,height=4) 		
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_4_Micro_ASV.png"),width=4,height=4,units="in",res=400) 	
}		
par(mar=c(5.1,4.2,4.1,2.1))

micro_file = "assay_new_setup/derived_data_and_code/microbiome/data_output/Count_sign_ASVs.csv"
micro_dat <- read.csv(micro_file,	header=T,	row.names=1, stringsAsFactors=T)

ploooot <- barplot(rev(t(micro_dat)),      beside=T, horiz=T, col=opaque_col1(col_vec[3]), main= "Microbiome",	xlab="Number of ASVs (N: 78)")
text(	x=-max(micro_dat) * .09,
		y = ploooot, 
		labels=colnames(custom_dat), xpd=NA,
		srt=45, adj=c(.9,1.5))
		
dev.off()

####################################
#------------- Metabolites ## 2  plot ## From tobias LCMS
####################################

if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_5_Metabo_mod.pdf"),width=7,height=6) 		
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_5_Metabo_mod.png"),width=7,height=6,units="in",res=400) 	
}		
par(mfrow=c(2,2), mar=c(4,2,4,1),oma=c(1,3,0,0))	

#grep data for CTmax samples:
#metab_file = "assay_new_setup/derived_data_and_code/metabolites/data/Significance_table_ANOVA_metabolites_LCMS.csv"
metab_file = "assay_new_setup/derived_data_and_code/metabolites/data/Significance_table_ANOVA_metabolites_ALL.csv"
metab_dat <- read.csv(metab_file,	header=T, stringsAsFactors=T, row.names=1)
metab_dat <- metab_dat[,-ncol(metab_dat)]
names(metab_dat)<-NULL
metab_dat_sub <- metab_dat[ grep("CTmax",rownames(metab_dat)), ]
row_to_grep <- grep("LCMS",rownames(metab_dat_sub))
# xlimis <- c( 0 , round_any(max(metab_dat_sub[row_to_grep,]),10,f=ceiling) )
xlimis <- c(0,45)
plotaa <- barplot(rev(as.matrix(metab_dat_sub[row_to_grep,])),xlim=xlimis,   beside=F, horiz=T, col=c(opaque_col1(col_vec[4])), main= "LCMS - CTmax",	
			xlab="Number of metabolites (N: 95)", names.arg=NULL, xpd=NA, axes=F)
labels_are=seq(0,45,10)	; axis(1, at=labels_are, labels=T, tcl=-.5)
labels_are=seq(0,45,5)	; labels_are[seq(1,NROW(labels_are),2)] <- NA	;	axis(1, at=labels_are, labels=NA, tcl=-.25)
text(	x=-max(metab_dat_sub[row_to_grep,]) * .09,
		y = plotaa, 
		labels=colnames(custom_dat), xpd=NA,
		srt=45, adj=c(.9,1.5))
mtext("a",side=3,line=0.5, at=-max(metab_dat_sub[row_to_grep,]) * .05,cex=1.2)

####################################
 #grep data for CRR samples:
metab_dat_sub <- metab_dat[ grep("CCR",rownames(metab_dat)), ]
row_to_grep <- grep("LCMS",rownames(metab_dat_sub))
# xlimis <- c( 0 , round_any(max(metab_dat_sub[row_to_grep,]),10,f=ceiling) )

barplot(rev(as.matrix(metab_dat_sub[row_to_grep,])), xlim=xlimis,  beside=F, horiz=T, col=c(opaque_col1(col_vec[5])), main= "LCMS - CCR",	
		xlab="Number of metabolites (N: 95)", names.arg=NULL, xpd=NA, axes=F)
labels_are=seq(0,45,10)	; axis(1, at=labels_are, labels=T, tcl=-.5)
labels_are=seq(0,45,5)	; labels_are[seq(1,NROW(labels_are),2)] <- NA	;	axis(1, at=labels_are, labels=NA, tcl=-.25)
mtext("b",side=3,line=0.5, at=-max(metab_dat_sub[row_to_grep,]) * .05,cex=1.2)


####################################
#------------- Metabolites ## 2 plots ##
####################################

	# number of genes showing patterns
#metab_file = "assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/Significance_table_ANOVA_metabolites_NMR.csv"
#metab_dat <- read.csv(metab_file,	header=T, stringsAsFactors=T, row.names=1)

metab_dat_sub <- metab_dat[grep("NMR",rownames(metab_dat)), ] 

'org_CCR <- 112
org_CCR <- 112
org_CCR <- 112
org_CCR <- 112
'
####################################
 #CTmax:
row_to_grep <- grep("CTmax",rownames(metab_dat_sub))
metab_dat_CTmax <- metab_dat_sub[row_to_grep,]
metab_dat_CTmax <- t(apply(metab_dat_CTmax,1,rev))
# xlimis <- c( 0 , round_any(max(colSums(metab_dat_CTmax)),10,f=ceiling) )
xlimis <- c(0,200)
print(xlimis)
colnames(metab_dat_CTmax)<-NULL
plotaa <- barplot(as.matrix(metab_dat_CTmax),  xlim=xlimis, beside=F, horiz=T, col="white", main= "NMR - CTmax",	
			xlab=paste0("Number of peaks"), names.arg=NULL, xpd=NA)
for (i in 1:ncol(metab_dat_CTmax)){
    xx = metab_dat_CTmax
    xx[,-i] <- NA
	names(xx)<-NULL
	if(0 %in% xx[1,]){
		print("yes")
		shading_lines=0
		}else{
		shading_lines=15
		}
	if( i != 3 ){ #ncol(metab_dat_CTmax[,-1])
		par(lwd=2) # add broad shading
		barplot(as.matrix(xx),xlim=xlimis,col=c(opaque_col1(col_vec[6])[i]), density=c(NA,NA),add=T, axes=F, beside=F,horiz=T,ann=F,border=F,names.arg=NULL) 
		barplot(as.matrix(xx),xlim=xlimis,col="white", density=c(shading_lines,0),add=T, axes=F, beside=F,horiz=T,ann=F,border=F,names.arg=NULL) 
		par(lwd=1)	 # add normal box
		barplot(as.matrix(xx),xlim=xlimis,col=NA,add=T, axes=F, beside=F,horiz=T,ann=F, names.arg=NULL) 
	}else{
		par(lwd=2) # add broad shading
		legend_text= paste0(lapply(strsplit(rownames(metab_dat_CTmax),"_"), "[[", 2), list(" (N: "," (N: "),list(240,454),list(")",")"))
		barplot(as.matrix(xx),xlim=xlimis,col=c(opaque_col1(col_vec[6])[i]), density=c(NA,NA),add=T, axes=F, beside=F,horiz=T,ann=F,border=F, legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text,bty = "n", x = "top", inset=-.15, ncol = 2,border=F,text.col="white", cex=.9)) 
		barplot(as.matrix(xx),xlim=xlimis,col="white", density=c(shading_lines,0),add=T, axes=F, beside=F,horiz=T,ann=F,border=F, legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text,bty = "n", x = "top", inset=-.15, ncol = 2,border=F,text.col="white", cex=.9)) 
		par(lwd=.5) 
		barplot(as.matrix(xx),xlim=xlimis,col=NA,add=T, axes=F, beside=F,horiz=T,ann=F,legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text,bty = "n", x = "top", inset=-.15, ncol = 2,border=T,col=NULL, cex=.9)) 
	}
} 
text(	x=-max(colSums(metab_dat_CTmax)) * .06,
		y = plotaa,
		labels=colnames(custom_dat), xpd=NA,
		srt=45, adj=c(.9,1.5))
mtext("c",side=3,line=0.5, at=-max(colSums(metab_dat_CTmax)) * .04,cex=1.2)

####################################
 #CCR:
row_to_grep <- grep("CCR",rownames(metab_dat_sub))
metab_dat_CCR <- metab_dat_sub[row_to_grep,]
metab_dat_CCR <- t(apply(metab_dat_CCR,1,rev))
# xlimis <- c( 0 , round_any(max(colSums(metab_dat_CCR)),50,f=ceiling) )
colnames(metab_dat_CCR)<-NULL
barplot(as.matrix(metab_dat_CCR),  xlim=xlimis, beside=F, horiz=T, col="white", main= "NMR - CCR",	
				xlab=paste0("Number of peaks"),names.arg=NULL, xpd=NA)
for (i in 1:ncol(metab_dat_CCR)){
    xx = metab_dat_CCR
    xx[,-i] <- NA
	names(xx)<-NULL
	if(0 %in% xx[1,]){
		print("yes")
		shading_lines=0
		}else{
		shading_lines=15
		}
	if( i != 3 ){ #ncol(metab_dat_CCR[,-1])
		par(lwd=2) # add broad shading
		barplot(as.matrix(xx),xlim=xlimis,col=c(opaque_col1(col_vec[7])[i]), density=c(NA,NA),add=T, axes=F, beside=F,horiz=T,ann=F,border=F,names.arg=NULL) 
		barplot(as.matrix(xx),xlim=xlimis,col="white", density=c(shading_lines,0),add=T, axes=F, beside=F,horiz=T,ann=F,border=F,names.arg=NULL) 
		par(lwd=1)	 # add normal box
		barplot(as.matrix(xx),xlim=xlimis,col=NA,add=T, axes=F, beside=F,horiz=T,ann=F,names.arg=NULL) 
	}else{
		par(lwd=2) # add broad shading
		legend_text= paste0(lapply(strsplit(rownames(metab_dat_CCR),"_"), "[[", 2), list(" (N: "," (N: "),list(240,454),list(")",")"))
		barplot(as.matrix(xx),xlim=xlimis,col=c(opaque_col1(col_vec[7])[i]), density=c(NA,NA),add=T, axes=F, beside=F,horiz=T,ann=F,border=F, legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text,bty = "n", x = "top", inset=-.15, ncol = 2,border=F,text.col="white", cex=.9)) 
		barplot(as.matrix(xx),xlim=xlimis,col="white", density=c(shading_lines,0),add=T, axes=F, beside=F,horiz=T,ann=F,border=F, legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text,bty = "n", x = "top", inset=-.15, ncol = 2,border=F,text.col="white", cex=.9)) 
		par(lwd=1) 
		barplot(as.matrix(xx),xlim=xlimis,col=NA,add=T, axes=F, beside=F,horiz=T,ann=F,legend=T,names.arg=NULL,
			args.legend = list(xpd=NA, legend=legend_text, bty = "n", x = "top", inset=-.15, ncol = 2,border=T,col=NULL, cex=.9)) 
	}
	print(par()$mar)
}
mtext("d",side=3,line=0.5, at=-max(colSums(metab_dat_CCR)) * .05,cex=1.2)

dev.off()
print(metab_dat)
