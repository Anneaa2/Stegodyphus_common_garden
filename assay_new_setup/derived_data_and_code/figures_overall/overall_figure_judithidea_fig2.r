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

#fig_types=	"individuel_panels"# "overview" # "individuel_panels"
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

#par( mfrow=c(4,3))
#col_vec =  c("red", "orange", "yellow","green",
#				"turquoise","blue","purple","magenta4","maroon","tomato3")
				
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



#################################################### IDEA 2 #################################################

    # violin plots with boxplots embedded
library(ggplot2)
library(cowplot)

####################################
#------------- CTmax---------------#
####################################
ToD_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/ToD/ToD_nest-temp_mean.csv", header=T,stringsAsFactors=T) 
jitter_move=seq(-.75,.75,length.out=4)
jitter_vec <- ToD_final$Population  ; levels(jitter_vec) <- jitter_move
cols_vec <- ToD_final$Population  ; levels(cols_vec) <- popcols
violin_df_ct <- data.frame(ToD_final$Population, paste(ToD_final$Population,ToD_final$Treatment,sep="_"), ToD_final$ToD_mean, ToD_final$Treatment, ToD_final$Treatment+as.numeric(as.character(jitter_vec)), cols_vec)  
colnames(violin_df_ct) <- c("Population","Poptemp","CTmax","Acclimation_temperature","Acclimation_temperature_mod", "color_vector")
linear_ctmax <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_select_mod_NOTMEANC_intc_slope.csv")	# import linear estimates

CCR_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/CCR/CCR_nest-temp_mean.csv", header=T,stringsAsFactors=T) 
jitter_vec <- CCR_final$Population  ; levels(jitter_vec) <- jitter_move
cols_vec <- CCR_final$Population  ; levels(cols_vec) <- popcols
violin_df_ccr <- data.frame(CCR_final$Population, paste(CCR_final$Population,CCR_final$Treatment,sep="_"), CCR_final$CCR_mean, CCR_final$Treatment, CCR_final$Treatment+as.numeric(as.character(jitter_vec)), cols_vec)  
colnames(violin_df_ccr) <- c("Population","Poptemp","CCRTemp","Acclimation_temperature","Acclimation_temperature_mod", "color_vector")
linear_ccr <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/results/ccr_select_mod_NOTMEANC_intc_slope.csv")	# import linear estimates

	# Betta (K) Otavi Stampriet different from zero --> break other lines
'Betta - intercept and trendline estimate
[1] 48.1001310  0.1040479
Karasburg - intercept and trendline estimate
[1] 49.62960290  0.02210014
Otavi - intercept and trendline estimate
[1] 48.57027707  0.08615041
Stampriet - intercept and trendline estimate
[1] 49.31309213  0.04165179'
line_type_vec <- c(1,2,1,1)
ctmax_plot <- ggplot( aes(x=Acclimation_temperature_mod, y=CTmax, group=Poptemp, fill=Population, colour=Population), 
			data=violin_df_ct ) + 
    theme_bw() 
for( pop_is in levels(ToD_final$Population) ){
	pop_nr <- which( levels(ToD_final$Population) == pop_is )
	message(pop_is, " - intercept and trendline estimate")
	# print(coef(summary(lm(CTmax~Acclimation_temperature, data=violin_df_ct[ which( violin_df_ct$Population == pop_is ), ])))[2,1:2]) # printing standard error
	# cf <- coef(lm(CTmax~Acclimation_temperature, data=violin_df_ct[ which( violin_df_ct$Population == pop_is ), ]))
	cf <- linear_ctmax[ grep(pop_is, linear_ctmax[,1]), 2]
	print(cf)
	ctmax_plot <- ctmax_plot + geom_abline(intercept=cf[1], slope=cf[2], linetype = line_type_vec[pop_nr], color=levels(cols_vec)[pop_nr], lwd=.6)	}
ctmax_plot <- ctmax_plot +
	geom_line(lwd=0, show.legend = T, alpha=0) + #only to add legend with lines (thus plotting nothing with no color)
	geom_violin(width=.85, position = position_dodge(), show.legend=F) +	# add violins
	scale_colour_manual( values=alpha(c(popcols,rep("black",6)),.7)) +	# change color of lines on violins
    geom_boxplot(width=.1, colour="grey25", alpha=.5, position = position_dodge(), outlier.shape = NA, lwd=.1, show.legend=F) + # change fill alpha using alpha
	scale_fill_manual( values=alpha(c(popcols,rep("black",6)),.18)) + # change violin fill color
	scale_x_continuous(breaks=seq(15,29,2), labels=seq(15,29,2)) +  
	theme(axis.text.x=element_text(colour="transparent", family="Helvetica")) + #c("black","transparent","black","transparent","black","black","transparent","black"))) +
	theme(legend.position="none", legend.title=element_blank(), legend.spacing.x =unit(.35,"lines")) +
	theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
	ylab(expression(paste("CTmax (", degree,"C)"))) +
	theme(plot.margin=unit(c(0,.3,0,.4), "lines"), axis.title.y=element_text(size=10.5, family="Helvetica")) +
	guides( color=guide_legend(override.aes = list(alpha=1, lwd=.7))) +	# changes legend alpha and lwd
	coord_cartesian(xlim = c(14.5,29.5))

	# (B) Karasburg (O) Stampriet different from zero --> break other lines
'Betta - intercept and trendline estimate:
[1] 12.30775464  0.02775537
Karasburg - intercept and trendline estimate:
[1] 7.6325820 0.2005924
Otavi - intercept and trendline estimate:
[1] 12.78086642  0.03271787
Stampriet - intercept and trendline estimate:
[1] 9.3483056 0.1383383'
line_type_vec <- c(2,1,2,1)
ccr_plot <-ggplot( aes(x=Acclimation_temperature_mod, y=CCRTemp, group=Poptemp, fill=Population, colour=Population), 
			data=violin_df_ccr ) + 
    theme_bw()
for( pop_is in levels(CCR_final$Population) ){
	pop_nr <- which( levels(CCR_final$Population) == pop_is )
	message(pop_is, " - intercept and trendline estimate:")
	#print(coef(summary(lm(CCRTemp~Acclimation_temperature, data=violin_df_ccr[ which( violin_df_ccr$Population == pop_is ), ])))[2,1:2]) # printing standard error
	#cf <- coef(lm(CCRTemp~Acclimation_temperature, data=violin_df_ccr[ which( violin_df_ccr$Population == pop_is ), ]))
	cf <- linear_ccr[ grep(pop_is, linear_ccr[,1]), 2]
	print(cf)
	ccr_plot <- ccr_plot + geom_abline(intercept=cf[1], slope=cf[2], linetype = line_type_vec[pop_nr], color=levels(cols_vec)[pop_nr], lwd=.6)	}
ccr_plot <- ccr_plot +
    geom_violin(width=.85, position = position_dodge()) +
	scale_colour_manual( values=alpha(c(popcols,rep("black",6)),.7)) +
	geom_boxplot(width=.1, colour="grey25", alpha=.5, position = position_dodge(), outlier.shape = NA, lwd=.1) +
	scale_fill_manual( values=alpha(c(popcols,rep("black",6)),.18)) +
	scale_x_continuous(breaks=seq(15,29,2), labels=seq(15,29,2)) +  
	theme(axis.text.x=element_text(colour=c("black","transparent","black","transparent","black","black","transparent","black"), family="Helvetica")) +
	theme(legend.position="none") +
	theme(axis.title.x=element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
	ylab(expression(paste("CCRTemp (", degree,"C)"))) +
	theme(plot.margin=unit(c(0,.3,0,.4), "lines"), axis.title=element_text(size=10.5, family="Helvetica")) +
	coord_cartesian(xlim = c(14.5,29.5))

if (	output_type == "pdf"	)	{ 	
    pdf(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_2_phenoTemp1_idea2.pdf"),width=7,height=5) 		
}else if ( 	output_type == "png"	){
    png(paste0("assay_new_setup/derived_data_and_code/figures_overall/",coltest_type,"figure_2_phenoTemp1_idea2.png"),width=7,height=5,units="in",res=400) 	
}
common_legend <- get_legend(ctmax_plot + theme(legend.position="bottom"))
combined_plot <- plot_grid(common_legend, ctmax_plot, ccr_plot, ncol=1, rel_heights = c(.2,1,1), labels=c(" ","a","b"), label_y=1.06)
ggdraw( add_sub(combined_plot, label=expression(paste("Acclimation temperature (", degree,"C)")),x=.55,y=.5, size=10.5))
dev.off()






