#!/usr/bin/env Rscript
##################!/usr/bin/env conda run -n my_env Rscript
# cd /home/anneaa/spider2/faststorage/assay_study_anne/metabolites ; 

# conda activate metabo_env ; R
# source("/faststorage/project/spider2/assay_study_anne/metabolites/data/analyzing_metabolome.r")
# source("assay_new_setup/derived_data_and_code/metabolites/data/analyzing_metabolome.r")

write_to_picturetype="pdf"    #
# write_to_picturetype="png"    #

library(MetaboAnalystR)
library(tools)
#install.packages("viridis")
library(viridis)
library(scales)
library(magick)		# importing pictures from pdf
library(grid)		# final figure layout
library(gridExtra)	#final figure layout
library(vegan)
library(pls) #functions for metaboalalyst
library(callr)
# setwd("/faststorage/project/spider2/assay_study_anne")	
source("assay_new_setup/derived_data_and_code/metabolites/data/plot_functions_metabo.r")
##
## OBS maybe setup plotting functions in another script and run it up here. Then employ functions below.
## plot both PCA and PLSDAs with location and treatment

metabofilesare <- c("/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CCR.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CCR.txt",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")

popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("B", "K", "N", "S", "O") #c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")

#AQ_CTmax.txt"
#ORG_CTmax.txt

setwd("/faststorage/project/spider2/assay_study_anne")	
datname_vec <- vector(); popsig_vec <-  vector(); tempsig_vec <-  vector()
pop_temp_sig_vec  <-  vector() ; intsig_vec <-  vector() ; pop_temp_int_sig_vec <-  vector()
#metabofile=metabofilesare[5]
for ( metabofile in metabofilesare ) {
  tryCatch({
    #-------------- Make metadata file
    overall_metabo_type <- unlist(strsplit( 	rev(unlist(strsplit(metabofile,"/")))[2] 	, "_"))[1]
	if( overall_metabo_type == "NMR" ){
		datname <- paste0(overall_metabo_type,"_",file_path_sans_ext(basename(metabofile))) #AQ_CCR
	}else if(overall_metabo_type == "LCMS"){
		datname <- paste0(overall_metabo_type,"_",rev(unlist(strsplit(file_path_sans_ext(basename(metabofile)),"_")))[1]) #AQ_CCR
	}
    dirnameis <- dirname(metabofile) #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten"
    new_dir_name <- paste0(dirnameis,"/",datname,"_results")  #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CCR_results"

    dir.create(new_dir_name)
    setwd(new_dir_name)
    
	if( overall_metabo_type == "NMR" ){
		datasub <- read.table(metabofile, header=T,sep=",")
		datasub1 <- datasub[,-(5:ncol(datasub))]	;    datasub <- datasub[,-c(3,4)]
		print(datname)	;    dim(datasub)
		write.table(datasub, file=paste0(datname,"_mainDat.txt"),row.names = F)
		write.table(datasub1, file=paste0(datname,"_meta.txt"),row.names = F)
		rm(datasub) 	; rm(datasub1)
		main_dat_file <- paste0(datname,"_mainDat.txt")
		meta_dat_file <-  paste0(datname,"_meta.txt")
	}else if( overall_metabo_type == "LCMS"){
		main_dat_file <- metabofile
		meta_dat_file <- dir(path=dirnameis, pattern=paste0("metadata_lcms_onlyNamed_",rev(unlist(strsplit(datname,"_")))[1]),full.name=T)
	}
	
    mSet<-InitDataObjects("pktable", "ts", FALSE) # peak intestisty table, timeseries/two factor, not paired
    mSet<-SetDesignType(mSet, "multi")
    mSet<-Read.TextDataTs(mSet, main_dat_file, "rowts")
    mSet<-ReadMetaData(mSet, meta_dat_file)	;    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet)	;	    mSet<-SanityCheckMeta(mSet, 1)
    mSet<-SetDataTypeOfMeta(mSet)	;	    mSet<-SanityCheckData(mSet)
    message("Sanity check passed, now normalizing.")
    mSet <- FilterVariable(mSet, "none", "F", 25)
    mSet <- PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "SumNorm", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)
      # normalised by sum, and pareto scaled
    mSet<-PlotNormSummary(mSet, paste0(datname,"_norm_0_"), "pdf", 72, width=NA)
    mSet<-PlotSampleNormSummary(mSet, paste0(datname,"_snorm_0_"), "pdf", 72, width=NA)
    write.csv(mSet$dataSet$norm, file=paste0("normalized_data_",datname,".csv"))
    message("Normalization has been done and plotted. ")

    mSet<-ANOVA2.Anal(mSet, 0.05, "fdr", "multi", 1, 1)
    mSet<-PlotANOVA2(mSet, paste0(datname,"_aov_abint_"), "pdf", 72, width=NA)
    message("Performed Two-way Anova ")

    # get color and pch vectors
	
	if (overall_metabo_type == "NMR"){
		tempcols <- plasma(NROW(levels(mSet$dataSet$meta.info$Acclimation))) ; names(tempcols) <- c("15","19","23","25","29")
		popcols1 <- popcols[which(names(popcols) %in% levels(mSet$dataSet$meta.info$Population ))]
	}else if(overall_metabo_type == "LCMS"){
		tempcols <- plasma(NROW(levels(mSet$dataSet$meta.info$Acclimation))) ; names(tempcols) <- c("19","29")
		popcols1 <- popcols[which(names(popcols) %in%  substr(levels(mSet$dataSet$meta.info$Population ),1,1))]
	}
    mSet <- PCA.Anal(mSet)
    mSet <- PlotPCAPairSummary(mSet, paste0("pca_pair_0_overall_"), "pdf", 72, width=NA, 5)
    mSet <- PlotPCAScree(mSet, paste0("pca_scree_0_overall_"), "pdf", 72, width=NA, 5)
    mSet <- PlotPCA2DScore(mSet, paste0("pca_score2d_0_overall_"), "pdf", 72, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPCABiplot(mSet, paste0("pca_biplot_0_overall_"), "pdf", 72, width=NA, 1,2)
    #mSet<-Perform.ASCA(mSet, a=1, b=1, x=2, res=2)
    #mSet<-CalculateImpVarCutoff(mSet, spe.thresh = 0.05, lev.thresh = 0.9)

    mSet <- PLSR.Anal(mSet, reg=TRUE)
    #mSet <- PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5) # modify if want right colors
    mSet <- PlotPLS2DScore_mod(mSet, imgName = paste0("pls_pop_score2d_0_overall_"), factor_to_use="pop", col_vector=popcols1, "pdf", 72, width=NA, inx1=1,inx2=2,reg=0.95,show=1,grey.scale=0)
    mSet <- PlotPLS2DScore_mod(mSet, imgName = paste0("pls_temp_score2d_0_overall_"), factor_to_use="temp", col_vector=tempcols, "pdf", 72, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPLS2DScore_mod(mSet, imgName = paste0("pls_pop_score2d_0_overall_"), factor_to_use="pop", col_vector=popcols1, "png", 400, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPLS2DScore_mod(mSet, imgName = paste0("pls_temp_score2d_0_overall_"), factor_to_use="temp", col_vector=tempcols, "png", 400, width=NA, 1,2,0.95,1,0)
    #mSet <- PlotPLS2DScore(mSet, imgName = "pls_pop_score2d_0_", "png", 400, width=NA, 1,2,0.95,1,0)
    #mSet <- PlotPLS2DScore(mSet, imgName = "pls_temp_score2d_0_", "png", 400, width=NA, 1,2,0.95,1,0)
    message("PLS-DA analysis has been run.")

    # make new dataframe with data (if making asca):
    'pop_sig <- read.csv("Sig_features_Model_a.csv")
    temp_sig <- read.csv("Sig_features_Model_b.csv")
    int_sig <- read.csv("Sig_features_Model_ab.csv")'

    # make new dataframe with data (if making ANOVA):
    ANOVA_dat <- read.csv("anova_between_sbj.csv")
    pop_sig <- NROW(which(ANOVA_dat$Population.adj.p. <.05))
    temp_sig <- NROW(which(ANOVA_dat$Acclimation.adj.p. <.05))
    pop_temp_sig <- NROW(intersect(which(ANOVA_dat$Acclimation.adj.p. <.05), which(ANOVA_dat$Population.adj.p. <.05) ))
    int_sig <- NROW(which(ANOVA_dat$Interaction.adj.p. <.05))
    pop_temp_int_sig <- NROW(intersect( intersect(which(ANOVA_dat$Acclimation.adj.p. <.05), which(ANOVA_dat$Population.adj.p. <.05) ), which(ANOVA_dat$Interaction.adj.p. <.05)))

    datname_vec <- append(datname_vec,datname)
    popsig_vec <- append(popsig_vec,pop_sig)
    tempsig_vec <- append(tempsig_vec,temp_sig)
    pop_temp_sig_vec <- append(pop_temp_sig_vec,pop_temp_sig)
    intsig_vec <- append(intsig_vec,int_sig)
    pop_temp_int_sig_vec <- append(pop_temp_int_sig_vec,pop_temp_int_sig)

    #PCA of pop-effects
    sig_features_pop <- ANOVA_dat[which(ANOVA_dat$Population.adj.p. <.05),1]
    if (NROW(sig_features_pop) > 0){
      Pca_dat_pop_features <- (mSet$dataSet$norm[,which(colnames(mSet$dataSet$norm) %in% sig_features_pop)])
      write.csv(Pca_dat_pop_features,file=paste0("Signi_features_pop_",datname,".csv"),row.names = T)
      
   
      pca_pop_features <- prcomp(Pca_dat_pop_features)
      colvec1 = factor_is = substr(mSet$dataSet$cls,1,1)
      col_vector1 = popcols1[as.character(colvec1)] # creates named vector of colors - use as input.
      
      write.csv(data.frame(rownames(Pca_dat_pop_features),factor_is,mSet$dataSet$facB, tempcols[as.character(mSet$dataSet$facB)], col_vector1),file=paste0(new_dir_name,"/","Signi_features_effectANDcolor_pop_",datname,".csv"),row.names = F)
      
      #pdf(paste0(datname,"_pca_featuers_pop_plot.pdf"),width=5,height=5)
      pdf(paste0(datname,"_pca_featuers_pop_plot_mod.pdf"),width=2.5,height=5)
        pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
      dev.off()
      #png(paste0(datname,"_pca_featuers_pop_plot.png"),width=5,height=5,units="in",res=400)
      png(paste0(datname,"_pca_featuers_pop_plot_mod.png"),width=2.5,height=5,units="in",res=400)
        pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
      dev.off()

      # Also plot PLS-DA for only 
      message("Also making PLS-DA for only Population related features")
       # note, standardize the cls, to minimize the impact of categorical to numerical impact
      reg=0.95
      if(reg){         cls <- scale(as.numeric(mSet$dataSet$cls))[,1];
      }else{        cls <- model.matrix(~mSetObj$dataSet$cls-1);      }
      comp.num <- dim(Pca_dat_pop_features)[1]-1
      datmat <- as.matrix(Pca_dat_pop_features)
        #as.matrix(mSetObj$dataSet$norm);
      mSet$analSet$plsr_POP <- pls::plsr(cls~datmat, method='oscorespls')#, ncomp=comp.num);
      mSet$analSet$plsr_POP$reg <- reg;
      mSet$analSet$plsr_POP$loading.type <- "all";
      
      PlotPLS2DScore_mod_for_pop_temp(mSet, manual_mset_analysis="plsr_POP",imgName = paste0("pls_popDat_pop_score2d_0_"), width=3, factor_to_use="pop", col_vector=popcols1, "pdf", 72, 1,2,0.95,1,0)
      PlotPLS2DScore_mod_for_pop_temp(mSet, manual_mset_analysis="plsr_POP",imgName = paste0("pls_popDat_pop_score2d_0_"), width=3, factor_to_use="pop", col_vector=popcols1, "png", 400, 1,2,0.95,1,0)
		}

    sig_features_temp <- ANOVA_dat[which(ANOVA_dat$Acclimation.adj.p. <.05),1]
    if (NROW(sig_features_temp) > 0){
      Pca_dat_temp_features <- (mSet$dataSet$norm[,which(colnames(mSet$dataSet$norm) %in% sig_features_temp)])
      write.csv(Pca_dat_temp_features,file=paste0("Signi_features_temp_",datname,".csv"),row.names = T)
      #pca_pop_features <- mSet$analSet$pca
      pca_pop_features <- prcomp(Pca_dat_temp_features)

      colvec1 = factor_is = mSet$dataSet$facB
      #levels(colvec1)[order(levels(colvec1))] <- popcols1[order(names(popcols1))]
      col_vector1 <- tempcols[as.character(colvec1)] # creates named vector of colors - use as input.
      
      write.csv(data.frame(rownames(Pca_dat_temp_features),factor_is,mSet$dataSet$cls,col_vector1, popcols1[as.character(mSet$dataSet$cls)]),file=paste0("Signi_features_effectANDcolor_temp_",datname,".csv"),row.names = F)
      
      #pdf(paste0(datname,"_pca_featuers_temp_plot.pdf"),width=5,height=5)
      pdf(paste0(datname,"_pca_featuers_temp_plot_mod.pdf"),width=2.5,height=5)
        pca_plot_metabo(PCA_element=pca_pop_features, factor_to_use=factor_is, color_vector=col_vector1)
      dev.off()
      #png(paste0(datname,"_pca_featuers_temp_plot.png"),width=5,height=5,units="in",res=400)
      png(paste0(datname,"_pca_featuers_temp_plot_mod.png"),width=2.5,height=5,units="in",res=400)
        pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
      dev.off()

       # Also plot PLS-DA for only 
      message("Also making PLS-DA for only Temperature related features")
       # note, standardize the cls, to minimize the impact of categorical to numerical impact
      reg=0.95
      if(reg){         cls <- scale(as.numeric(mSet$dataSet$cls))[,1];
      }else{        cls <- model.matrix(~mSetObj$dataSet$cls-1);      }
      comp.num <- dim(Pca_dat_temp_features)[1]-1
      datmat <- as.matrix(Pca_dat_temp_features)
        #as.matrix(mSetObj$dataSet$norm);
      mSet$analSet$plsr_TEMP <- pls::plsr(cls~datmat, method='oscorespls')#, ncomp=comp.num);
      mSet$analSet$plsr_TEMP$reg <- reg;
      mSet$analSet$plsr_TEMP$loading.type <- "all";
      
      PlotPLS2DScore_mod_for_pop_temp(mSet, manual_mset_analysis="plsr_TEMP",imgName = paste0("pls_tempDat_temp_score2d_0_"), width=3, factor_to_use="temp", col_vector=tempcols, "pdf", 72, 1,2,0.95,1,0)
      PlotPLS2DScore_mod_for_pop_temp(mSet, manual_mset_analysis="plsr_TEMP",imgName = paste0("pls_tempDat_temp_score2d_0_"), width=3, factor_to_use="temp", col_vector=tempcols, "png", 400, 1,2,0.95,1,0)
    }
    print(new_dir_name)
    message("Number of peaks/features:",NROW(mSet$dataSet$prenorm.feat.nms),"\n")
    rm(mSet)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  setwd("/faststorage/project/spider2/assay_study_anne")
}
count_sig_dat <- data.frame(   datname_vec,        popsig_vec,     tempsig_vec,     pop_temp_sig_vec,   intsig_vec,   pop_temp_int_sig_vec)
#write.csv(count_sig_dat,file=paste0(dirnameis,"/Significance_table_asca_metabolites_NMR.csv"),row.names = F)
write.csv(count_sig_dat,file=paste0("assay_new_setup/derived_data_and_code/metabolites/results/","Significance_table_ANOVA_metabolites_ALL.csv"),row.names = F)



###########################
### Plot all metabolites with effect of pop or temp
############################

plotit_simple_func <- function(metabolite_name, dat_vector, temp_vector, pop_col_vector, pop_vector, scatter_vector, overall_metabo_type ){
  plot(temp_vector+scatter_vector, dat_vector, xaxt='n', xpd=NA, 
    main=metabolite_name , col=pop_col_vector, xlab=NA, ylab=NA, las=2)
  if( overall_metabo_type == "LCMS" ){
    axis(1, at=c(21, 23, 25, 27), labels=NA, tcl=-.25, col.ticks="grey70")
    axis(1, at=c(19, 29), labels=NA)
  }else{
    axis(1, at=c(17, 21, 27), labels=NA, tcl=-.25, col.ticks="grey70")
    axis(1, at=c(15,19,23,25,29), labels=NA)
  }
  for ( pop_is in c("B","K","O","S") ){
    dat_b <- dat_vector[ which( pop_vector == pop_is ) ]
    temp_b <- temp_vector[ which( pop_vector == pop_is ) ]
    abline( lm( dat_b ~ temp_b ), col=unique(pop_col_vector[ which( pop_vector == pop_is ) ]) )
  }
}

# write_to_picturetype="pdf"    #
# write_to_picturetype="png"    #

  # files must have been through above loop, otherwise files used here does not exist.
metabofilesare <- c("/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CCR.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CCR.txt",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("B", "K", "N", "S", "O") #c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")
popcols <- popcols[grep("N",names(popcols), invert=T)]
if( any(names(popcols) == "B") ) {names(popcols)[which(names(popcols) == "B")] <- "Betta"}
if( any(names(popcols) == "K") ) {names(popcols)[which(names(popcols) == "K")] <- "Karasburg"}
if( any(names(popcols) == "S") ) {names(popcols)[which(names(popcols) == "S")] <- "Stampriet"}
if( any(names(popcols) == "O") ) {names(popcols)[which(names(popcols) == "O")] <- "Otavi"}

setwd("/faststorage/project/spider2/assay_study_anne")	
datname_vec <- vector(); popsig_vec <-  vector(); tempsig_vec <-  vector()
pop_temp_sig_vec  <-  vector() ; intsig_vec <-  vector() ; pop_temp_int_sig_vec <-  vector()
#metabofile=metabofilesare[5]
for ( metabofile in metabofilesare ) {
    #-------------- get names of type and swich to folder
  overall_metabo_type <- unlist(strsplit( 	rev(unlist(strsplit(metabofile,"/")))[2] 	, "_"))[1]
	
  if( overall_metabo_type == "NMR" ){
		datname <- paste0(overall_metabo_type,"_",file_path_sans_ext(basename(metabofile))) #AQ_CCR
	}else if(overall_metabo_type == "LCMS"){
		datname <- paste0(overall_metabo_type,"_",rev(unlist(strsplit(file_path_sans_ext(basename(metabofile)),"_")))[1]) #AQ_CCR
	}
  dirnameis <- dirname(metabofile) #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten"
  new_dir_name <- paste0(dirnameis,"/",datname,"_results")  #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CCR_results"
  setwd(new_dir_name)

  # loop over temp significant and pop significant
  # get files with signi_features
  files_to_iterate <- list.files(pattern="Signi_features_[p|t]" )
  for ( temp_pop in files_to_iterate ) {
    intensity_dat <- read.csv(file=temp_pop,row.names = 1)
    intensity_dat <- intensity_dat[,order(colnames(intensity_dat))]
    if( grepl("pop", temp_pop) ){ pop_or_temp <- "pop" 
    }else if ( grepl("temp", temp_pop) ){ pop_or_temp <- "temp" }
    datname_2 <- paste0(datname,"_",pop_or_temp,"Sign")
    intensity_meta <- read.csv( file=list.files(pattern=paste0("Signi_features_effectANDcolor_",pop_or_temp)), row.names=1 ) 
      # plotting loop
    for ( metabolite in colnames(intensity_dat) ){
      nr_is <- which(colnames(intensity_dat)==metabolite)
      if ( nr_is %in% seq(1,ncol(intensity_dat),20) ){
        if ( write_to_picturetype == "pdf" ){
          pdf( file=paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/results/plots_trends_significant_",datname_2,"_",nr_is,".pdf"),width=8.5, height=5.5 )
        }else if ( write_to_picturetype == "png" ){
          png( file=paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/results/plots_trends_significant_",datname_2,"_",nr_is,".png"),width=8.5, height=5.5, units="in", res=400 ) }
        # png()
        par(mfrow=c(4,5), mar=c(1.1, 2.1,1.8,1.1), oma=c(3,2.1,1.8,0))
        if ( (nr_is + 20) > ncol(intensity_dat) ){ 
          number_rows = ceiling( (ncol(intensity_dat) - nr_is + 1) / 5 )
          par(mfrow=c(number_rows,5)) 
        }
      }
      data_vector <- intensity_dat[ ,which(colnames(intensity_dat)==metabolite) ]
      if ( pop_or_temp == "pop" ){ popcolumn_is <- 1 ; tempcolumn_is <- 2 
      }else if ( pop_or_temp == "temp" ){ popcolumn_is <- 2 ; tempcolumn_is <- 1 }
      temper_vector <- intensity_meta[,tempcolumn_is]
      pop_colorista <- intensity_meta[,4]
      pop_vectorista <- intensity_meta[,popcolumn_is]
      if ( is.na(unique(pop_colorista)) ){
          pop_colorista <- as.factor(pop_vectorista)    ;     levels(pop_colorista) <- popcols[ sort(names(popcols)) ]    ; pop_colorista <- as.character(pop_colorista)  }
      if ( nchar(pop_vectorista)[1] > 1 ){      pop_vectorista <- substr(pop_vectorista,1,1)      }
      scatter_vectorista <- as.factor(pop_vectorista)   ; levels(scatter_vectorista) <- seq(-.5,.5,length.out=4)
      scatter_vectorista <- as.numeric( as.character( scatter_vectorista) )
      
      plotit_simple_func( metabolite, dat_vector=data_vector, temp_vector=temper_vector, pop_col_vector=pop_colorista, pop_vector=pop_vectorista, scatter_vector=scatter_vectorista, overall_metabo_type=overall_metabo_type )
      if( nr_is %in% seq(1,ncol(intensity_dat),5)){
        title(ylab="Metabolite Intensity", xpd=NA)
      }
      if( ( nr_is %% 20 ) %in% c(0,seq(16,20,1))    ||     nr_is %in% c((ncol(intensity_dat)-4):ncol(intensity_dat)) ){ # if last row
        title(xlab=expression(paste("Acclimation Temperature (", degree,"C)")), xpd=NA)
        if( overall_metabo_type == "LCMS" ){
          axis(1, at=c(19,29), labels=c(19,29), tick=F, gap.axis=.4)
        }else{
          axis(1, at=c(15,19,25,29), labels=c(15,19,25,29), tick=F, gap.axis=.4)
          axis(1, at=c(23), labels=c(23), tick=F, gap.axis=.4)
        }
      }
      if ( nr_is %in% seq(0,ncol(intensity_dat),20) || nr_is == ncol(intensity_dat) ){ # if 20 plots are drawn, or it is last column
        par(new=T, mfrow=c(1,1))
        plot(0,type='n',axes=FALSE,ann=FALSE)
        #title(main=datname_2, xpd=NA, line=5.5)
        legend("top", legend=names(popcols), lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9, inset=c(0,-.19))
        dev.off()
      }
    }
  } 
}

setwd("/faststorage/project/spider2/assay_study_anne")	






###########################
### Plot all! LCMS metabolites
############################



plotit_simple_func <- function(metabolite_name, dat_vector, temp_vector, pop_col_vector, pop_vector, scatter_vector, overall_metabo_type ){
  plot(temp_vector+scatter_vector, dat_vector, xaxt='n', xpd=NA, 
    main=metabolite_name , col=pop_col_vector, xlab=NA, ylab=NA, las=2)
  if( overall_metabo_type == "LCMS" ){
    axis(1, at=c(21, 23, 25, 27), labels=NA, tcl=-.25, col.ticks="grey70")
    axis(1, at=c(19, 29), labels=NA)
  }else{
    axis(1, at=c(17, 21, 27), labels=NA, tcl=-.25, col.ticks="grey70")
    axis(1, at=c(15,19,23,25,29), labels=NA)
  }
  for ( pop_is in c("B","K","O","S") ){
    dat_b <- dat_vector[ which( pop_vector == pop_is ) ]
    temp_b <- temp_vector[ which( pop_vector == pop_is ) ]
    abline( lm( dat_b ~ temp_b ), col=unique(pop_col_vector[ which( pop_vector == pop_is ) ]) )
  }
}

# write_to_picturetype="pdf"    #
# write_to_picturetype="png"    #

  # files must have been through above loop, otherwise files used here does not exist.
metabofilesare <- c("/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/ORG_CCR.txt",
                    "/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/AQ_CCR.txt",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
					"/faststorage/project/spider2/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("B", "K", "N", "S", "O") #c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")
popcols <- popcols[grep("N",names(popcols), invert=T)]
if( any(names(popcols) == "B") ) {names(popcols)[which(names(popcols) == "B")] <- "Betta"}
if( any(names(popcols) == "K") ) {names(popcols)[which(names(popcols) == "K")] <- "Karasburg"}
if( any(names(popcols) == "S") ) {names(popcols)[which(names(popcols) == "S")] <- "Stampriet"}
if( any(names(popcols) == "O") ) {names(popcols)[which(names(popcols) == "O")] <- "Otavi"}

setwd("/faststorage/project/spider2/assay_study_anne")	
datname_vec <- vector(); popsig_vec <-  vector(); tempsig_vec <-  vector()
pop_temp_sig_vec  <-  vector() ; intsig_vec <-  vector() ; pop_temp_int_sig_vec <-  vector()
#metabofile=metabofilesare[5]
for ( metabofile in metabofilesare ) {
    #-------------- get names of type and swich to folder
  overall_metabo_type <- unlist(strsplit( 	rev(unlist(strsplit(metabofile,"/")))[2] 	, "_"))[1]
	
  if( overall_metabo_type == "NMR" ){
		datname <- paste0(overall_metabo_type,"_",file_path_sans_ext(basename(metabofile))) #AQ_CCR
	}else if(overall_metabo_type == "LCMS"){
		datname <- paste0(overall_metabo_type,"_",rev(unlist(strsplit(file_path_sans_ext(basename(metabofile)),"_")))[1]) #AQ_CCR
	}
  dirnameis <- dirname(metabofile) #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten"
  new_dir_name <- paste0(dirnameis,"/",datname,"_results")  #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CCR_results"
  setwd(new_dir_name)

  # change to get "normalized_data_LCMS_CCR"
  files_to_iterate <- list.files(pattern="normalized_data" )
  #for ( temp_pop in files_to_iterate ) {
    intensity_dat <- read.csv(file=files_to_iterate,row.names = 1)
    intensity_dat <- intensity_dat[,order(colnames(intensity_dat))]
    #if( grepl("pop", temp_pop) ){ pop_or_temp <- "pop" 
    #}else if ( grepl("temp", temp_pop) ){ pop_or_temp <- "temp" }
    #datname_2 <- paste0(datname,"_",pop_or_temp,"Sign")
    metafile <- list.files(pattern=paste0("Signi_features_effectANDcolor_"))[1]
    intensity_meta <- read.csv( file=metafile, row.names=1 ) 
      # plotting loop
    for ( metabolite in colnames(intensity_dat) ){
      nr_is <- which(colnames(intensity_dat)==metabolite)
      if ( nr_is %in% seq(1,ncol(intensity_dat),20) ){
        if ( write_to_picturetype == "pdf" ){
          pdf( file=paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/results/plots_trends_all_",datname,"_",nr_is,".pdf"),width=8.5, height=5.5 )
        }else if ( write_to_picturetype == "png" ){
          png( file=paste0("/home/anneaa/spider2/faststorage/assay_study_anne/assay_new_setup/derived_data_and_code/metabolites/results/plots_trends_all_",datname,"_",nr_is,".png"),width=8.5, height=5.5, units="in", res=400 ) }
        # png()
        par(mfrow=c(4,5), mar=c(1.1, 2.1,1.8,1.1), oma=c(3,2.1,1.8,0))
        if ( (nr_is + 20) > ncol(intensity_dat) ){ 
          number_rows = ceiling( (ncol(intensity_dat) - nr_is + 1) / 5 )
          par(mfrow=c(number_rows,5)) 
        }
      }
      data_vector <- intensity_dat[ ,which(colnames(intensity_dat)==metabolite) ]
      if ( grepl("pop", metafile) ){ popcolumn_is <- 1 ; tempcolumn_is <- 2 
      }else if ( grepl("temp", metafile) ) { popcolumn_is <- 2 ; tempcolumn_is <- 1 }
      temper_vector <- intensity_meta[,tempcolumn_is]
      pop_colorista <- intensity_meta[,4]
      pop_vectorista <- intensity_meta[,popcolumn_is]
      if ( is.na(unique(pop_colorista)) ){
          pop_colorista <- as.factor(pop_vectorista)    ;     levels(pop_colorista) <- popcols[ sort(names(popcols)) ]    ; pop_colorista <- as.character(pop_colorista)  }
      if ( nchar(pop_vectorista)[1] > 1 ){      pop_vectorista <- substr(pop_vectorista,1,1)      }

      scatter_vectorista <- as.factor(pop_vectorista)   ; levels(scatter_vectorista) <- seq(-.5,.5,length.out=4)
      scatter_vectorista <- as.numeric( as.character( scatter_vectorista) )
      
      plotit_simple_func( metabolite, dat_vector=data_vector, temp_vector=temper_vector, pop_col_vector=pop_colorista, pop_vector=pop_vectorista, scatter_vector=scatter_vectorista, overall_metabo_type=overall_metabo_type )
      if( nr_is %in% seq(1,ncol(intensity_dat),5)){
        title(ylab="Metabolite Intensity", xpd=NA)
      }
      if( ( nr_is %% 20 ) %in% c(0,seq(16,20,1))    ||     nr_is %in% c((ncol(intensity_dat)-4):ncol(intensity_dat)) ){ # if last row
        title(xlab=expression(paste("Acclimation Temperature (", degree,"C)")), xpd=NA)
        if( overall_metabo_type == "LCMS" ){
          axis(1, at=c(19,29), labels=c(19,29), tick=F, gap.axis=.4)
        }else{
          axis(1, at=c(15,19,25,29), labels=c(15,19,25,29), tick=F, gap.axis=.4)
          axis(1, at=c(23), labels=c(23), tick=F, gap.axis=.4)
        }
      }
      if ( nr_is %in% seq(0,ncol(intensity_dat),20) || nr_is == ncol(intensity_dat) ){ # if 20 plots are drawn, or it is last column
        par(new=T, mfrow=c(1,1))
        plot(0,type='n',axes=FALSE,ann=FALSE)
        #title(main=datname_2, xpd=NA, line=5.5)
        legend("top", legend=names(popcols), lty=1, lwd=1.1, col=popcols, xpd=NA, bty="n",horiz=T, cex=.9, inset=c(0,-.19))
        dev.off()
      }
    }
  #} 
}

setwd("/faststorage/project/spider2/assay_study_anne")	


#################################
# Make panel of PLSDA plots:
###############################


# Combine plots: Should fit a page:

# CCR panel
pan1<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
pan2<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_AQ_CTmax_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

#pan3<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_ORG_CTmax_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
pan4<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_ORG_CTmax_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

pan5<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/LCMS_CTmax_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
pan6<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/LCMS_CTmax_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

  # CCR panel
panb1<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
panb2<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_AQ_CCR_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

#panb3<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_ORG_CCR_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
panb4<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/NMR_Kirsten/NMR_ORG_CCR_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

panb5<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/LCMS_CCR_results/pls_popDat_pop_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)
panb6<- image_read_pdf("assay_new_setup/derived_data_and_code/metabolites/data/LCMS_tobias/LCMS_CCR_results/pls_tempDat_temp_score2d_0_dpi72.pdf")				# read in map and pca plots as picture files)

panel.blank <- image_read_pdf("assay_new_setup/derived_data_and_code/figures_overall/empty_fig.pdf")		#magick



img1<-rasterGrob(pan1, interpolate=TRUE)    # Make them into raster objects
img2<-rasterGrob(pan2, interpolate=TRUE)		# grid
#img3<-rasterGrob(pan3, interpolate=TRUE)		# grid
img4<-rasterGrob(pan4, interpolate=TRUE)    # Make them into raster objects
img5<-rasterGrob(pan5, interpolate=TRUE)		# grid
img6<-rasterGrob(pan6, interpolate=TRUE)		# grid

imgb1<-rasterGrob(panb1, interpolate=TRUE)    # Make them into raster objects
imgb2<-rasterGrob(panb2, interpolate=TRUE)		# grid
#imgb3<-rasterGrob(panb3, interpolate=TRUE)		# grid
imgb4<-rasterGrob(panb4, interpolate=TRUE)    # Make them into raster objects
imgb5<-rasterGrob(panb5, interpolate=TRUE)		# grid
imgb6<-rasterGrob(panb6, interpolate=TRUE)		# grid


blank_img <-rasterGrob(panel.blank, interpolate=TRUE)		# grid


# Panel 1 CTmax
pdf("assay_new_setup/derived_data_and_code/metabolites/results/metab_ctmax_plsda_panel.pdf",height=7.5,width=5.5)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(img1,img2,blank_img,img4,img5,img6,ncol=2, vp=viewport(width=0.9, height=0.9)) #,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("NMR - AQ",vjust=-20.5, hjust=-2.3, rot=90)		# when rotated - horizontal and vertical is opposite?								#grid
grid.text("NMR - ORG",vjust=-20.5, hjust=0.39, rot=90)
grid.text("LCMS",vjust=-20.5, hjust=5.35, rot=90)
grid.text("Population effect",vjust=-29, hjust=1.5)
grid.text("Temperature effect",vjust=-29, hjust=-.4)
grid.text("a", vjust=-27, hjust=26)
grid.text("b", vjust=-27, hjust=-1)
grid.text("c", vjust=-8, hjust=-1)
grid.text("d", vjust=11, hjust=26)
grid.text("e", vjust=11, hjust=-1)
dev.off()

png("assay_new_setup/derived_data_and_code/metabolites/results/metab_ctmax_plsda_panel.png",height=7.5, width=5.5, units="in",res=400)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(img1,img2,blank_img,img4,img5,img6,ncol=2, vp=viewport(width=0.9, height=0.9)) #,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("NMR - AQ",vjust=-21.5, hjust=-2.5, rot=90)		# when rotated - horizontal and vertical is opposite?								#grid
grid.text("NMR - ORG",vjust=-21.5, hjust=0.39, rot=90)
grid.text("LCMS",vjust=-21.5, hjust=5.35, rot=90)
grid.text("Population effect",vjust=-30.5, hjust=1.5)
grid.text("Temperature effect",vjust=-30.5, hjust=-.4)
grid.text("a", vjust=-28, hjust=26)
grid.text("b", vjust=-28, hjust=-1)
grid.text("c", vjust=-8, hjust=-1)
grid.text("d", vjust=11, hjust=26)
grid.text("e", vjust=11, hjust=-1)
dev.off()

# panel 2 CCR

pdf("assay_new_setup/derived_data_and_code/metabolites/results/metab_ccr_plsda_panel.pdf",height=7.5, width=5.5)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(imgb1,imgb2,blank_img,imgb4,imgb5,imgb6,ncol=2, vp=viewport(width=0.9, height=0.9)) #,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("NMR - AQ",vjust=-20.5, hjust=-2.3, rot=90)		# when rotated - horizontal and vertical is opposite?								#grid
grid.text("NMR - ORG",vjust=-20.5, hjust=0.39, rot=90)
grid.text("LCMS",vjust=-20.5, hjust=5.35, rot=90)
grid.text("Population effect",vjust=-29, hjust=1.5)
grid.text("Temperature effect",vjust=-29, hjust=-.4)
grid.text("a", vjust=-27, hjust=26)
grid.text("b", vjust=-27, hjust=-1)
grid.text("c", vjust=-8, hjust=-1)
grid.text("d", vjust=11, hjust=26)
grid.text("e", vjust=11, hjust=-1)
dev.off()

png("assay_new_setup/derived_data_and_code/metabolites/results/metab_ccr_plsda_panel.png",height=7.5, width=5.5, units="in",res=400)   # Opden pdf for writing and plot objects on grid according to width
grid.arrange(imgb1,imgb2,blank_img,imgb4,imgb5,imgb6,ncol=2, vp=viewport(width=0.9, height=0.9)) #,widths=unit(c(7/2,7/4,7/4), c("in")))		#gridExtra
grid.text("NMR - AQ",vjust=-21.5, hjust=-2.5, rot=90)		# when rotated - horizontal and vertical is opposite?								#grid
grid.text("NMR - ORG",vjust=-21.5, hjust=0.39, rot=90)
grid.text("LCMS",vjust=-21.5, hjust=5.35, rot=90)
grid.text("Population effect",vjust=-30.5, hjust=1.5)
grid.text("Temperature effect",vjust=-30.5, hjust=-.4)
grid.text("a", vjust=-28, hjust=26)
grid.text("b", vjust=-28, hjust=-1)
grid.text("c", vjust=-8, hjust=-1)
grid.text("d", vjust=11, hjust=26)
grid.text("e", vjust=11, hjust=-1)
dev.off()






