## This script will make figure 1, which consist of three elements:
	# 1: a map of the populations on the left.
	# 2: three temperature plots: min, mean, max over the year. in the midsection
	# 3: a neighbour joining tree to the right. imported from picture file.
	

# conda activate assay_magick1; R 
## updated magick in another environment, since I could not install it in assay_copy.
# librariet etc:
#setwd("/faststorage/project/spider2/assay_study_anne/figures_overall/")
library(tidyr)			# data magic
library(dendextend)		# plotting phylogeny
library(spData)		#dat for map
data(world)			#map
library(tmap) 		#map
library(sf)			#stuff for map
library(sp)			#raster object for map
library(magick)		# importing pictures from pdf
library(grid)		# final figure layout
library(gridExtra)	#final figure layout
library(scales)		# function for alpha color
library(phylogram)		# plotting phylogeny
library(car)		# stats stuff
#library(ape)			# phylogeny stuff

# functions:
	colMax <- function(tdata) sapply(tdata, max, na.rm = TRUE)
	colMIN <- function(tdata) sapply(tdata, min, na.rm = TRUE)

# adjustable variables:
output_type = "pdf"
	
#plot map

height_is=3.5
width_is=7


# create spacial object dataframe
#popcols=c("steelblue1","blue","tan4","green3","red") # B, K, N, S, E/O
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")
lati <- read.csv(list.files("/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/nest_interpolation_data/",pattern = paste0("latitude_pops.csv"),full.names = T,recursive = F),header=T,sep=";",row.names=1)
longi <- read.csv(list.files("/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/nest_interpolation_data/",pattern = paste0("longitude_pops.csv"),full.names = T,recursive = F),header=T,sep=";",row.names=1)
names_are <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Gobabis","Otavi")


coordi_dat <- cbind(rownames(lati),lati,longi,names_are)		; colnames(coordi_dat) <- c("names1","latitude","longitude","names_are")
coordi_dat <- coordi_dat[ -which( rownames(coordi_dat) == "G"),]

coords_are <- coordi_dat[,c("longitude","latitude")]
data_is <- droplevels(data.frame(coordi_dat,popcols))
data_is$popcols <- factor(data_is$popcols, levels=popcols)
crs <- CRS("+proj=longlat +ellps=WGS84") #("+init=epsg:28992") #WGS84
crs <- CRS("+proj=longlat +datum=WGS84") #("+init=epsg:28992") #WGS84
	#"+proj=longlat +datum=WGS84"
spdf <- SpatialPointsDataFrame(coords      = coords_are,
                               data        = data_is, 
                               proj4string = crs)
							   
Africa<-world[which(world$continent=="Africa"),]
SouthernAfrica<-Africa[c(which(Africa$name_long=="Botswana"),
							which(Africa$name_long=="South Africa"),
							which(Africa$name_long=="Namibia"),
							which(Africa$name_long=="Lesotho"),
							which(Africa$name_long=="Swaziland")),]
							#which(Africa$name_long=="Zimbabwe"),
							#which(Africa$name_long=="Mozambique"))]


PopMap <- tm_shape(SouthernAfrica)+
	tm_borders(alpha=.5)+
	tm_text("name_long",just=c(.65,-.3),alpha=.8,size=.65)+ # iso_a3
	tm_graticules(alpha=.1)+#,labels.size=1.1)+
	tm_layout(frame="black",frame.lwd=2.4)+#,outer.margins=c(0.1, 0.1, 0, 0))+
tm_shape(spdf)+
	tm_dots(col="popcols",size=.4, palette=popcols, legend.show=F) + 
	tm_text("names_are",col="popcols",just=c(1.25,.9),size=.71, palette=popcols, legend.col.show=F)
#PopMap

tmap_save(PopMap,paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part1_figure1.pdf"),height=height_is,width=width_is/2,units="in")#,outer.margins = c(.021,.01,.039,.02))  # Save map in appropriate format
tmap_save(PopMap,paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part1_figure1.png"),height=height_is,width=width_is/2,units="in")#,outer.margins = c(.021,.01,.039,.02))  # Save map in appropriate format






## import data for weather:
# import data:
	diris<-"/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/nest_interpolation_data/means_thinplate/"
	#enddir<-"/home/anneaa/spider2/faststorage/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/graphs/"
		#//uni.au.dk/Users/au308918/Documents/Spiderlab/Sources_of_variation/weather_data_SOV/NewLocClim/means_thinplate/"
		
	files<- list.files(diris,pattern = paste0("_daily"),full.names = T,recursive = T)
	files<-files[!grepl("allpoints",files)]
	files<-files[!grepl(".pdf",files)]
	files<-files[!grepl("dist.csv",files)]
	files<-files[grepl("T",files)]
	files<-files[!grepl("PET",files)]
	' basename(files)
	[1] "Max. Temp._daily_means.csv" "Mean Temp._daily_means.csv"
	[3] "Min. Temp._daily_means.csv"
	'
	titlesare=c("Maximum","Mean","Minimum")
	
### plot 3 weather
	#popcols=c("steelblue1","blue","green3","red") # B, K, S, E/O
if (	output_type == "pdf"	)	{ 	
	#pdf(paste0("/faststorage/project/spider2/assay_study_anne/figures_overall/part3_figure1.pdf"),width=width_is/4,height=height_is) 		
	pdf(paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part3_figure1_mod.pdf"),width=width_is/4,height=height_is) 		
}else if ( 	output_type == "png"	){
	#png(paste0("/faststorage/project/spider2/assay_study_anne/figures_overall/part3_figure1.png"),width=width_is/4,height=height_is,units="in",res=400) 	
	png(paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part3_figure1_mod.png"),width=width_is/4,height=height_is,units="in",res=400) 	
}

par(mfrow=c(4,1),cex=.6, cex.lab=.9, cex.axis=.8, mar=c(.5,3.5,1,1), oma=c(3+1/7,0,0+1.3+.065,0), mgp=c(3,.8,0))
var_df=list()
for (data_file in files){
	typeis=sapply(strsplit(basename(data_file),"_"),"[[",1)
	#Maybe adjust names
	message(titlesare[which(files==data_file)])
	message(typeis)
	tdata=read.csv(data_file, header=T,row.names=1)
	colnames(tdata)[8]<-gsub("_E","_O",colnames(tdata)[8])
	tdata <- tdata[,grep("date|day_nr|_B|_S|_O|_K",colnames(tdata))]
	popcols_4 <- popcols[grep("B|S|O|K",names(popcols))]
	
	maxofcols<- colMax(tdata[,3:ncol(tdata)])
	minofcols<- colMIN(tdata[,3:ncol(tdata)])
	highestnr<-max(maxofcols)
	lowestnr<-min(minofcols)
	  
	## work on removing x axis on two plots here
	plot(tdata[,2],tdata[,3], col=popcols_4[1],#ylab = type, xlab = xlab_ved[which(files==data_file)],
		lwd=3, ylim = c(lowestnr-2,highestnr+2),type = "l",ylab='',xlab='', xaxt='n',, yaxt='n', cex.axis=.9, las=2) #
	axis(side=2, labels=T, tcl=-.3, las=2) 
	axis(side=1,labels=F, tcl=-.3) 	
	title(titlesare[which(files==data_file)], line = -1,cex.main=.9,font=1)
	
	for (datacol in 4:ncol(tdata)) {
		lines(tdata[,datacol],col=popcols_4[datacol-2],lwd=3)
	}
	if( titlesare[which(files==data_file)] == "Maximum" || titlesare[which(files==data_file)] == "Minimum" ){
		var_df[[titlesare[which(files==data_file)]]] <- tdata
	}
}
tdata <- var_df[[1]]
tdata[,3:ncol(tdata)] <- var_df[["Maximum"]][,3:ncol(var_df[["Maximum"]])] - var_df[["Minimum"]][,3:ncol(var_df[["Minimum"]])]
maxofcols<- colMax(tdata[,3:ncol(tdata)])	; minofcols<- colMIN(tdata[,3:ncol(tdata)])
highestnr<-max(maxofcols)	; lowestnr<-min(minofcols)
plot(tdata[,2],tdata[,3], col=popcols_4[1],#ylab = type, xlab = xlab_ved[which(files==data_file)],
	lwd=3, ylim = c(lowestnr-2,highestnr+2),type = "l",ylab='',xlab='', xaxt='n', yaxt='n', cex.axis=.9) #
axis(side=2, labels=T, tcl=-.3, las=2)	; 	axis(side=1,labels=T, tcl=-.3)
title("Variation", line = -1,cex.main=.9,font=1)
for (datacol in 4:ncol(tdata)) {
		lines(tdata[,datacol],col=popcols_4[datacol-2],lwd=3, ylim = c(lowestnr-2,highestnr+2))
	}
title(xlab="Day of the year",xpd=NA, line=2)

par(new=TRUE,mfrow=c(1,1))
plot(0,type='n',axes=FALSE,ann=FALSE)
title(ylab=expression(paste("Temperature (", degree,"C)")),xpd=NA, line=2.5, cex.lab=.59)
dev.off()

# get stats:
data_list <- lapply(as.list(files), read.csv, header=T, row.names=1)
data_list <- append(data_list,list(data.frame(data_list[[1]][,c(1,2)], data_list[[1]][,-c(1,2)]-data_list[[3]][,-c(1,2)])))
names(data_list) <- c("Max temperature", "Mean temperature", "Min temperature", "Var temperature")
data_list_sub <- lapply(data_list, function (x) x[, grep("date|day_nr|_B|_S|_E|_K",colnames(x))])
data_list_sub <- lapply(data_list_sub, function (x){ colnames(x) <- c("date", "day_nr", "B", "K", "S", "O"); return(x)})
#make into long format:
data_list_sub <- lapply(data_list_sub, function(x) gather(x, population, temperature, "B":"O", factor_key=TRUE))
#lapply(data_list_sub, function(x) aggregate(temperature ~ population,x,median))
data_list_sub <- lapply(data_list_sub, function(x){ medians <- aggregate(temperature ~ population,x,median); x$medians <- x$population; levels(x$medians) <- medians$temperature; return(x)})
# residuals
data_list_sub <- lapply(data_list_sub, function(x){ x$res <- abs(x$temperature-as.numeric(as.character(x$medians))); return(x)})

# if I were to test (which we decided against - because what biological meaning does significance have?)
### # test homogeniety of variance (if deviate from normality)
	#shapiro.test()
	# deviate from normality

# If I should test if Variances are similar?
 # Bartlets:
 	#Performs Bartlett's test of the null that the variances in each of the groups (samples) are the same.
	# assumes normality
# levenes test is more robust to deviations from normality
# leveneTest()		# Levene is an anova on residuals, we can run the anova on residual ourself and do tukeyhsd
	# TukeyHSD(aov())
# flingers test is nonparametric, thus does not assume normality
	# Performs a Fligner-Killeen (median) test of the null that the variances in each of the groups (samples) are the same.
# fligner.test()


#### Test similarity of means
# t test assumes normality and equal variance (pairwise)
# welch assumes non-equal variance 
# anova does the same but can handle more groups (still normality)
# wilcoxon is nonparametrical
# kruskall-wallis is for groups, and nonparametrical
	#kruskal.test()
	# pairwise.wilcox.test()



###################################
###########	 phylogeny  ###########
###################################

#popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols

if (	output_type == "pdf"	)	{ 	
	pdf(paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part2_figure1.pdf"),width=width_is/4,height=height_is) 		
}else if ( 	output_type == "png"	){
	png(paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part2_figure1.png"),width=width_is/4,height=height_is,units="in",res=400) 	
}
par(xpd=NA, mar=c(2,1,0,3.5), cex.axis=.58)
x <- read.dendrogram(text = "(((Betta:0.00089274,Otavi:0.00090302):0.00010593,(Karasburg:0.00092387,Stampriet:0.00097664):0.00021567):0.0001,Ndumo:0.00137902); ")
xlims_are <- c(0.00137902,0)
toOBKS <- 0.0001		;	toOB <- 0.00010593	;		toKS <- 0.00021567

pop_order <- vector()
for( each in labels(x) ){
	pop_order <- append(pop_order,grep(each,names(popcols)))
}
popcols_phylo <- popcols[pop_order]
rev(x) %>% 
	set("labels_col", rev(popcols_phylo)) %>% 
	set("labels_cex", .75) %>% 
	plot(horiz=T, edgePar=list(lwd=1),xpd=NA, axes=F)
	# axis, a third length FIX 
xlims_new <- xlims_are/3
xlims_new2 <- c(xlims_new[1],			xlims_new[2]+((xlims_new[1]-xlims_new[2])/2),			xlims_new[2])
factor_xlim <- xlims_new[1]
axis(1,at=xlims_new2+factor_xlim, labels=c(NA,format(round(xlims_new2[1],digits=4), scientific=F),NA),lwd.ticks=0, mgp=c(3,-.1,0), line=-.2)
#axis(1,at=xlims_new2+factor_xlim, labels=c(NA,0.0005,NA),lwd.ticks=0, mgp=c(3,-.1,0), line=-.2)
#nodelabels("200",2)
text(xlims_are[1]-toOBKS-toOB-0.00017,4.5,labels="100", cex=.48)
text(xlims_are[1]-toOBKS-toKS-0.00017,2.5,labels="100", cex=.48)

dev.off()


# Combine plots: Should fit half a page:
panel.left <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part1_figure1.pdf")				# read in map and pca plots as picture files)
#panel.mid <- image_read("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/assay_phylogeny.png")				# read in map and pca plots as picture files)
#panel.mid <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/phylogeny_long.pdf")				# read in map and pca plots as picture files)
panel.mid <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part2_figure1.pdf")				# read in map and pca plots as picture files)
#panel.right <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part3_figure1.pdf")		#magick
panel.right <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/part3_figure1_mod.pdf")		#magick
panel.blank <- image_read_pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/empty_fig.pdf")		#magick

img1<-rasterGrob(panel.left, interpolate=TRUE)    # Make them into raster objects
img2<-rasterGrob(panel.mid, interpolate=TRUE)		# grid
img3<-rasterGrob(panel.right, interpolate=TRUE)		# grid
blank_img <-rasterGrob(panel.blank, interpolate=TRUE)		# grid


text_hjust=.65		; 		Text_size=9
pdf("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/Figure1_mod.pdf",height=3.5,width=7)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(img1,img2,img3,ncol=3,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("a",vjust=-13,hjust=35)										#grid
grid.text("b",vjust=-13,hjust=-1)
grid.text("c",vjust=-13,hjust=-23)
dev.off()

png("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/figures_overall/Figure1_mod.png",height=3.5,width=7, units="in",res=400)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(img1,img2,img3,ncol=3,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("a",vjust=-13,hjust=35)										#grid
grid.text("b",vjust=-13,hjust=-1)
grid.text("c",vjust=-13,hjust=-23)
dev.off()



