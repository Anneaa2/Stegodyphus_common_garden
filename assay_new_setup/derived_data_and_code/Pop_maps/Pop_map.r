#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript

#R

library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(ggplot2) # tidyverse data visualization package
library(shiny)   # for web applications
data(World)
#library(gridGraphics)
#library(gridExtra)
#library(cowplot)
#library(grid)
#library(magick) #for combining in panels

outdir<-"/faststorage/project/spider2/assay_study_anne/Pop_maps/"


# create spacial object dataframe
#get data for spacial object:
popcols=c("steelblue1","blue","green3","red") # B, K, N, S, G, E/O

lati <- read.csv(list.files("/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/nest_interpolation_data/",pattern = paste0("latitude_pops.csv"),full.names = T,recursive = F),header=T,sep=";",row.names=1)
rownam <- rownames(lati);rownam <- rownam[-c(3,5)];colis<- colnames(lati); lati <- as.data.frame(lati[-c(3,5),]);rownames(lati)<-rownam; colnames(lati)<- colis

longi <- read.csv(list.files("/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/climate_var_values/nest_interpolation_data/",pattern = paste0("longitude_pops.csv"),full.names = T,recursive = F),header=T,sep=";",row.names=1)
colis1<- colnames(longi); longi <- as.data.frame(longi[-c(3,5),]); rownames(longi)<-rownam; colnames(longi)<- colis1

names_are <- c("Betta", "Karasburg", "Stampriet", "Otawi")
coordi_dat <- cbind(rownames(lati),lati,longi,names_are)

# Create spacial object with datum for plotting on map
coords_are <- coordi_dat[,c("longitude","latitude")]
data_is <- coordi_dat;data_is$popcols=popcols
crs <- CRS("+proj=longlat +ellps=WGS84")
crs <- CRS("+proj=longlat +datum=WGS84")
spdf <- SpatialPointsDataFrame(coords      = coords_are,
                               data        = data_is, 
                               proj4string = crs)
							   
							   
Africa<-World[which(World$continent=="Africa"),]
Namibia_map<-Africa[which(Africa$name=="Namibia"),]


PopMap <- tm_shape(Namibia_map)+
	tm_borders(alpha=.5)+
	tm_text("name",just=c(.5,-.3),alpha=.8)+
	tm_graticules(alpha=.1)+
	tm_layout(frame="black",frame.lwd=1.8)+
	#tm_xlab("Longitude",space=1.1)+
	#tm_ylab("Latitude",space=-.2)+
tm_shape(spdf)+
	tm_dots(size=.8,col="popcols") + #col=popcols,size=.5)+
	tm_text("names_are",col="popcols",just=c(1,-1.4))#,auto.placement=T)
PopMap


tmap_save(PopMap,paste0(outdir,"pop_map_assay.png"),height=8,width=6)
tmap_save(PopMap,paste0(outdir,"pop_map_assay.pdf"),height=8,width=6)

