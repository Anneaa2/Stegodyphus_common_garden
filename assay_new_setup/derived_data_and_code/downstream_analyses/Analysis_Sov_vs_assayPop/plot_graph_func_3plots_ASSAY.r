#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript


# function for plotting 3 barplots with results


plot_graph_now <- function(){
#SNPS
if( (COEF_ONLY==T & SIGN_ONLY==T) | (COEF_ONLY==T & COEF_AND_SIG==T) | (SIGN_ONLY==T & COEF_AND_SIG==T) ){message("ERROR! Choose just one dataset.")}
if( (COEF_ONLY==F & SIGN_ONLY==F & COEF_AND_SIG==F) ){message("ERROR! Choose at least one dataset.")}
if( (SIGN_ONLY==T) ){message("ERROR! Choose at least one dataset.")}
if( (COEF_AND_SIG==T) ){message("ERROR! Choose at least one dataset.")}
if( COEF_ONLY==T ) {	dataset_TOUSE <- dataset_above_perctile_sub		;			dataset_TOUSE_EXP <- dataset_exp_above_perctile_sub	}
#if( SIGN_ONLY==T ) {	dataset_TOUSE <- dataset_below_SIG_sub		;				dataset_TOUSE_EXP <- dataset_exp_below_SIG_sub}
#if( COEF_AND_SIG==T ) {	dataset_TOUSE <- dataset_above_perctile_andSIG_sub		;	dataset_TOUSE_EXP <- dataset_exp_above_perctile_andSIG_sub	}
if( expected_theoretical_number==T ) {	dataset_TOUSE_EXP <- dataset_exp_above_perctile_theory }
snpsmpMAX <- round(max(dataset_TOUSE[,1:3]),digits=-1)
microMAX <- (max(dataset_TOUSE[,4:ncol(dataset_TOUSE)])+1)
bari <- barplot(dataset_TOUSE[,1],		col=colorsare[1],
		ylab=paste0("N genes > ",as.numeric(percentile_is1)/100,"th percentile of null dist"),		main=mainis[1],
		cex.main=1,		ylim=c(0,snpsmpMAX), names.arg=rownames(dataset_TOUSE), yaxt="n",xpd=NA)
axis(side=2)
print(bari)
print(bari[c(2,4)])
#axis(side=1, labels=c("PC1","PC3","PC5"), at=c(1,3,5), lwd=0)
axis(side=1, labels=c("PC2","PC4"),at=bari[c(2,4)], lwd=0) ###################HERE IS THE ISSUE
#axis(side=1, labels=c("PC2","PC4"), at=c(2,4),lwd.ticks=0)
#mtext(annotateis[1],side=3,line=2,adj=adjustis,cex=1.2,outer=F)	
legend("bottom",	inset=c(0,-leg_inset),		 title="Corrected for:",
		legend=c( paste0(	"Pop structure")),
		fill=c(colorsare[1]),		bty="n", #y.intersp=1.2,
		border=c("black"),		xpd=NA,		horiz=T) #adjust inset to move legend further out of the plot
legend("bottom",inset=c(0,-(leg_inset+.1)),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,1])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, t(t(dataset_TOUSE_EXP[,1]))[rowis,each],	x1=first+.5, t(t(dataset_TOUSE_EXP[,1]))[rowis,each],	col=col_lines, lwd=1.5)
}	}	#abline(h=number_exp, col="orange", lwd=2)

#SMPS
bari <- barplot(t(dataset_TOUSE[,2:3]),		col=colorsare[2:3],
		main=mainis[2], ylab=paste0("N genes > ",as.numeric(percentile_is1)/100,"th percentile of null dist"),	
		cex.main=1,		ylim=c(0,snpsmpMAX), names.arg=rownames(dataset_TOUSE),beside=T, yaxt="n",xpd=NA)
axis(side=2,  labels=T) #at=snpsmp_axlab,)		
#mtext(annotateis[2],side=3,line=2,adj=adjustis,cex=1.2,outer=F)
legend("bottom",	inset=c(0,-leg_inset-leg_inset_micro), title="Corrected for:",
		legend=c( paste0(	"Pop structure"),	paste0(	"Cis genetic variation")), 
		text.col="white",		border="white",		#y.intersp=1.2, 	
		fill=c(colorsare[2],colorsare[3]), 
		bty="n",		xpd=NA,		ncol=1) #horiz=T)
legend("bottom",	inset=c(0,-leg_inset-leg_inset_micro),	title="Corrected for:",
		legend=c(paste0("Pop structure"),	paste0("Cis genetic variation")), 
		border=c("black","black"),		#y.intersp=1.2,
		fill=c(colorsare[2],colorsare[3]), 
		bty="n",		xpd=NA,		ncol=1,	#horiz=T,
		density=c(NA,Density_is[3]),		angle=c(angle_is[2],angle_is[3]))
legend("bottom",inset=c(0,-(leg_inset+.1)-leg_inset_micro),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,2:3])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, t(dataset_TOUSE_EXP[,2:3])[rowis,each],	x1=first+.5, y1=t(dataset_TOUSE_EXP[,2:3])[rowis,each],	col=col_lines, lwd=1.5)
		#print(t(dataset_TOUSE_EXP[,2:3])[rowis,each])
}	}		#abline(h=number_exp, col="orange", lwd=2)
#Micro cex
'
if( (COEF_ONLY=T & SIGN_ONLY=T) | (COEF_ONLY=T & COEF_AND_SIG=T) | (SIGN_ONLY=T & COEF_AND_SIG=T) ){message("ERROR! Choose just one dataset.")}
if( COEF_ONLY=T ) {dataset_TOUSE <- dataset_above_perctile_sub}
if( SIGN_ONLY=T ) {dataset_TOUSE <- dataset_below_SIG_sub}
if( COEF_AND_SIG=T ) {dataset_TOUSE <- dataset_above_perctile_andSIG_sub}
'

dataset_TOUSE[c(4,6)]
names_micro <- gsub("singleASV_","ASVs",colnames(dataset_TOUSE[c(4,6)]))
names_micro <- gsub("absolute","",names_micro)
names_micro <- gsub("relative","",names_micro)
names_micro <- gsub("genuslevel_","Genera",names_micro)
bari <- barplot(t(dataset_TOUSE[,c(4,6)]),		col=colorsare[c(4,5)],		ylim=c(0,microMAX),		beside=T,	main=mainis[3],		cex.main=1, ylab=paste0("N taxa (ASVs or genera) > null dist"),xpd=NA)			
barplot(t(dataset_TOUSE[,c(4,6)]),		col=c(colorsare[4],"white"),
	ylim=c(0,microMAX), names.arg=rownames(dataset_TOUSE),beside=T,
	density=c(		Density_is[4],					Density_is[5]),		
	angle=c(	angle_is[4],				angle_is[5]),
	add=T	,xpd=NA)	
#mtext(annotateis[3],side=3,line=2,adj=adjustis,cex=1.2,outer=F)
legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),	title= "Corrected for pop structure:", ncol=1,
	legend=c(		paste0(names_micro[1]),			paste0(names_micro[2])),
	text.col="white",		border="white",				fill=c(	colorsare[4],		colorsare[5]), 
	bty='n',				xpd=NA)#,			y.intersp=1.2)#,	text.width=3.8) #inset=c(0,inset_is-1), only when keyword
legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),		title= "Corrected for pop structure:", ncol=1,
	legend=c(		paste0(names_micro[1]),			paste0(names_micro[2])	),
	fill=c(	colorsare[4],			"white"		), 
	border=c("black","black"	),		bty="n",		
	density=c(		Density_is[4],			Density_is[5]	),		
	angle=c(		angle_is[4],			angle_is[5]	),	
	xpd=NA)#,	text.width=3.8) #y.intersp=1.2,		
legend("bottom",inset=c(0,-(leg_inset+.1)-leg_inset_micro),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,4:7])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, t(dataset_TOUSE_EXP[,c(4,6)])[rowis,each],	x1=first+.5, t(dataset_TOUSE_EXP[,c(4,6)])[rowis,each],	col=col_lines, lwd=1.5)
}	}	}	

	
	
plot_graph_now_microonly <- function(){
	if( (COEF_ONLY==T & SIGN_ONLY==T) | (COEF_ONLY==T & COEF_AND_SIG==T) | (SIGN_ONLY==T & COEF_AND_SIG==T) ){message("ERROR! Choose just one dataset.")}
	if( (COEF_ONLY==F & SIGN_ONLY==F & COEF_AND_SIG==F) ){message("ERROR! Choose at least one dataset.")}
	if( COEF_ONLY==T ) {	dataset_TOUSE <- dataset_above_perctile_sub		;			dataset_TOUSE_EXP <- dataset_exp_above_perctile_sub	}
	#if( SIGN_ONLY==T ) {	dataset_TOUSE <- dataset_below_SIG_sub		;				dataset_TOUSE_EXP <- dataset_exp_below_SIG_sub}
	#if( COEF_AND_SIG==T ) {	dataset_TOUSE <- dataset_above_perctile_andSIG_sub		;	dataset_TOUSE_EXP <- dataset_exp_above_perctile_andSIG_sub	}
	if( expected_theoretical_number==T ) {	dataset_TOUSE_EXP <- dataset_exp_above_perctile_theory }
	snpsmpMAX <- round(max(dataset_TOUSE[,1:3]),digits=-1)
	microMAX <- (max(dataset_TOUSE[,4:ncol(dataset_TOUSE)])+1)
	#Micro
	dataset_TOUSE[4:7]
	names_micro <- gsub("singleASV_","ASVs ",colnames(dataset_TOUSE[4:7]))
	names_micro <- gsub("absolute","Absolute Abund.",names_micro)
	names_micro <- gsub("relative","Relative Abund.",names_micro)
	names_micro <- gsub("genuslevel_","Genera ",names_micro)
	bari <- barplot(t(dataset_TOUSE[,4:7]),		col=colorsare[c(4,5,6,7)],		ylim=c(0,microMAX),		beside=T,	main=mainis[3],		cex.main=1, ylab=paste0("N taxa (ASVs or genera) > null dist"),xpd=NA)			
	barplot(t(dataset_TOUSE[,4:7]),		col=c(colorsare[4],colorsare[6],"white","white"),
		ylim=c(0,microMAX), names.arg=rownames(dataset_TOUSE),beside=T,
		density=c(		Density_is[4],					Density_is[5],					Density_is[6],					Density_is[7]	),		
		angle=c(	angle_is[4],				angle_is[5],				angle_is[6], 				angle_is[7]),
		add=T,	xpd=NA	)	
	legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),	title= "Corrected for pop structure:", ncol=2,
		legend=c(		paste0(names_micro[1]),			paste0(names_micro[2]),			paste0(names_micro[3]),			paste0(names_micro[4])	),
		text.col="white",		border="white",				fill=c(	colorsare[4],		colorsare[5],		colorsare[6],		colorsare[7]		), 
		bty='n',				xpd=NA)#,			y.intersp=1.2)#,	text.width=3.8) #inset=c(0,inset_is-1), only when keyword
	legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),		title= "Corrected for pop structure:", ncol=2,
		legend=c(		paste0(names_micro[1]),			paste0(names_micro[2]), 		paste0(names_micro[3]),			paste0(names_micro[4])		),
		fill=c(	colorsare[4],			colorsare[6],			"white", 			"white"			), 
		border=c("black","black","black","black"	),		bty="n",		
		density=c(		Density_is[4],			Density_is[5],				Density_is[6],				Density_is[7]	),		
		angle=c(		angle_is[4],			angle_is[5],				angle_is[6], 				angle_is[7]		),	
		xpd=NA)#,	text.width=3.8) #y.intersp=1.2,		
	legend("bottom",inset=c(0,-(leg_inset+.08)-leg_inset_micro),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
	#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,4:7])))
	for( each in 1:ncol(bari) ){
		for( rowis in 1:nrow(bari)	){
			first <- bari[rowis,each]
			segments(	first-.5, t(dataset_TOUSE_EXP[,4:7])[rowis,each],	x1=first+.5, t(dataset_TOUSE_EXP[,4:7])[rowis,each],	col=col_lines, lwd=1.5)
}	} 	}



# modify to fit!

plot_graph_now_mrm <- function(){
#SNPS
if( (COEF_ONLY==T & SIGN_ONLY==T) | (COEF_ONLY==T & COEF_AND_SIG==T) | (SIGN_ONLY==T & COEF_AND_SIG==T) ){message("ERROR! Choose just one dataset.")}
if( (COEF_ONLY==F & SIGN_ONLY==F & COEF_AND_SIG==F) ){message("ERROR! Choose at least one dataset.")}
if( COEF_ONLY==T ) {	dataset_TOUSE <- dataset_above_perctile_sub		;			dataset_TOUSE_EXP <- dataset_exp_above_perctile_sub	}
#if( SIGN_ONLY==T ) {	dataset_TOUSE <- dataset_below_SIG_sub		;				dataset_TOUSE_EXP <- dataset_exp_below_SIG_sub}
#if( COEF_AND_SIG==T ) {	dataset_TOUSE <- dataset_above_perctile_andSIG_sub		;	dataset_TOUSE_EXP <- dataset_exp_above_perctile_andSIG_sub	}
if( expected_theoretical_number==T ) {	dataset_TOUSE_EXP <- dataset_exp_above_perctile_theory }
snpsmpMAX <- round(max(dataset_TOUSE[,1:2]),digits=-1)
microMAX <- (max(dataset_TOUSE[,3:ncol(dataset_TOUSE)])+1)

bari <- barplot(dataset_TOUSE[,1],		col=colorsare[1], cex.names = .9,
		ylab=paste0("N genes > ",as.numeric(percentile_is1)/100,"th percentile of null dist"),		main=mainis[1],
		cex.main=1,		ylim=c(0,snpsmpMAX), names.arg=rownames(dataset_TOUSE), yaxt="n",xpd=NA)
axis(side=2)
#mtext(annotateis[1],side=3,line=2,adj=adjustis,cex=1.2,outer=F)	
legend("bottom",	inset=c(0,-leg_inset),	
		legend=c( paste0(	"Nucleotide variation")),
		fill=c(colorsare[1]),		bty="n", #y.intersp=1.2,
		border=c("black"),		xpd=NA,		horiz=T) #adjust inset to move legend further out of the plot
legend("bottom",inset=c(0,-(leg_inset+.1)),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,1])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, dataset_TOUSE_EXP[rowis,2],	x1=first+.5, dataset_TOUSE_EXP[rowis,2],	col=col_lines, lwd=1.5)
}	}	#abline(h=number_exp, col="orange", lwd=2)
#SMPS
bari <- barplot(dataset_TOUSE[,2],		col=colorsare[2],  cex.names = .9,
		main=mainis[2], ylab=paste0("N genes > ",as.numeric(percentile_is1)/100,"th percentile of null dist"),	
		cex.main=1,		ylim=c(0,snpsmpMAX), names.arg=rownames(dataset_TOUSE),beside=T, yaxt="n",xpd=NA)
axis(side=2,  labels=T) #at=snpsmp_axlab,)		
#mtext(annotateis[2],side=3,line=2,adj=adjustis,cex=1.2,outer=F)
'legend("bottom",	inset=c(0,-leg_inset-leg_inset_micro), 
		legend=c( paste0(	"DNA methylation variation")), 
		text.col="white",		border="white",		#y.intersp=1.2, 	
		fill=c(colorsare[2]), 
		bty="n",		xpd=NA,		ncol=1) #horiz=T)'
legend("bottom",	inset=c(0,-leg_inset),
		legend=c(paste0("DNA methylation variation")), 
		border=c("black","black"),		#y.intersp=1.2,
		fill=c(colorsare[2]), 
		bty="n",		xpd=NA,		ncol=1,	#horiz=T,
		)
legend("bottom",inset=c(0,-(leg_inset+.1)),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,2:3])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, dataset_TOUSE_EXP[rowis,2],	x1=first+.5, y1=dataset_TOUSE_EXP[rowis,2],	col=col_lines, lwd=1.5)
		#print(t(dataset_TOUSE_EXP[,2:3])[rowis,each])
}	}		#abline(h=number_exp, col="orange", lwd=2)
#Micro cex

dataset_TOUSE[3:6]
names_micro <- gsub("singleasvs_","ASVs ",colnames(dataset_TOUSE[3:6]))
names_micro <- gsub("absolute","Absolute",names_micro)
names_micro <- gsub("relative","Relative",names_micro)
names_micro <- gsub("genuslevel_","Genera ",names_micro)
bari <- barplot(t(dataset_TOUSE[,3:6]),		col=colorsare[3:6],		ylim=c(0,microMAX),		beside=T,	main=mainis[3],		cex.main=1, ylab=paste0("N taxa (ASVs or genera) > null dist"),xpd=NA, axisnames=F )			
barplot(t(dataset_TOUSE[,3:6]),		col=c(colorsare[3],colorsare[4],"white","white"),
	ylim=c(0,microMAX), names.arg=rownames(dataset_TOUSE),beside=T,  cex.names = .9,
	density=c(		Density_is[3],					Density_is[4],			Density_is[5],					Density_is[6]),
	angle=c(	angle_is[3],					angle_is[4],			angle_is[5],					angle_is[6]),
	add=T	,xpd=NA)	
#mtext(annotateis[3],side=3,line=2,adj=adjustis,cex=1.2,outer=F) colorsare[4],colorsare[6],"white","white"
legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),	ncol=2,
	legend=c(		paste0(names_micro[1]),			paste0(names_micro[2]),	paste0(names_micro[3]),			paste0(names_micro[4])),
	text.col="white",		border="white",				fill=c(	colorsare[3],		colorsare[4], 	colorsare[5],		colorsare[6]), 
	bty='n',				xpd=NA)#,			y.intersp=1.2)#,	text.width=3.8) #inset=c(0,inset_is-1), only when keyword
legend("bottom",inset=c(0,-leg_inset-leg_inset_micro),		ncol=2,
	legend=c(	paste0(names_micro[1]),			paste0(names_micro[2]),	paste0(names_micro[3]),			paste0(names_micro[4])	),
	fill=c(	colorsare[3],			colorsare[4],			"white"	,	"white"		), 
	border=c("black","black"	),		bty="n",		
	density=c(			Density_is[3],					Density_is[4],			Density_is[5],					Density_is[6]	),		
	angle=c(	angle_is[3],					angle_is[4],			angle_is[5],					angle_is[6]		),	
	xpd=NA)#,	text.width=3.8) #y.intersp=1.2,		
legend("bottom",inset=c(0,-(leg_inset+.1)-leg_inset_micro),	legend="N expected", 	lwd=1.5,	bty="n",	xpd=NA,	col=col_lines)
#number_exp <- unique(as.vector(as.matrix(dataset_TOUSE_EXP[,4:7])))
for( each in 1:ncol(bari) ){
	for( rowis in 1:nrow(bari)	){
		first <- bari[rowis,each]
		segments(	first-.5, t(dataset_TOUSE_EXP[,3:6])[rowis,each],	x1=first+.5, t(dataset_TOUSE_EXP[,3:6])[rowis,each],	col=col_lines, lwd=1.5)
}	}	}	



	# Test color blindness

library(colorspace)
colorblind_test <- function(colname_is=colname_is1) {
	if (fileout=="pdf"){
			pdf(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_TRI.pdf"),
				height=heightis,width=widthis)
		}
		if (fileout=="png"){
			png(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_TRI.png"),
				height=heightis,width=widthis,units="in",res=400)
		}
	colorsare <<- tritan(colorsare); col_lines <- tritan(col_lines)
		par(mfrow=c(1,3),oma=c(2.5,.5,0,0),mar=c(5.1,4.1,4.1,1.1))
	plot_graph_now(); dev.off()

	if (fileout=="pdf"){
			pdf(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_DEU.pdf"),
				height=heightis,width=widthis)
		}
		if (fileout=="png"){
			png(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_DEU.png"),
				height=heightis,width=widthis,units="in",res=400)
		}
	colorsare <<- deutan(colorsare); col_lines <- deutan(col_lines)
		par(mfrow=c(1,3),oma=c(2.5,.5,0,0),mar=c(5.1,4.1,4.1,1.1))
	plot_graph_now(); dev.off()
 
	if (fileout=="pdf"){
			pdf(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_PRO.pdf"),
				height=heightis,width=widthis)
		}
		if (fileout=="png"){
			png(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_PRO.png"),
				height=heightis,width=widthis,units="in",res=400)
		}
	colorsare <<- protan(colorsare); col_lines <- protan(col_lines)
		par(mfrow=c(1,3),oma=c(2.5,.5,0,0),mar=c(5.1,4.1,4.1,1.1))
	plot_graph_now(); dev.off()

	if (fileout=="pdf"){
			pdf(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_DESA.pdf"),
				height=heightis,width=widthis)
		}
		if (fileout=="png"){
			png(file=paste0("./figures_for_paper/sum_up_plots/",
				"COR_Coef_plots_sunimp2","_",analysis_type,dataset_is,"_",colname_is,"_DESA.png"),
				height=heightis,width=widthis,units="in",res=400)
		}
	colorsare <<- desaturate(colorsare); col_lines <- desaturate(col_lines)
		par(mfrow=c(1,3),oma=c(2.5,.5,0,0),mar=c(5.1,4.1,4.1,1.1))
	plot_graph_now(); dev.off()
}

	