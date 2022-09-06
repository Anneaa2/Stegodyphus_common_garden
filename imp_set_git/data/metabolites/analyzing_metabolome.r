#!/usr/bin/env Rscript
##################!/usr/bin/env conda run -n my_env Rscript
# cd /home/anneaa/spider2/faststorage/assay_study_anne/metabolites ; conda activate metabo_env ; R
# source("/faststorage/project/spider2/assay_study_anne/metabolites/data/analyzing_metabolome.r")

library(MetaboAnalystR)
library(tools)
#install.packages("viridis")
library(viridis)
library(scales)
library(vegan)
library(pls) #functions for metaboalalyst
library(callr)
setwd("/faststorage/project/spider2/assay_study_anne/metabolites/data")
source("/faststorage/project/spider2/assay_study_anne/metabolites/data/plot_functions_metabo.r")
##
## OBS maybe setup plotting functions in another script and run it up here. Then employ functions below.
## plot both PCA and PLSDAs with location and treatment

metabofilesare <- c("/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/ORG_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CTmax.txt",
                    "/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/ORG_CCR.txt",
                    "/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CCR.txt",
					"/faststorage/project/spider2/assay_study_anne/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
					"/faststorage/project/spider2/assay_study_anne/metabolites/data/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")

popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("B", "K", "N", "S", "O") #c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")

#AQ_CTmax.txt"
#ORG_CTmax.txt
setwd("/faststorage/project/spider2/assay_study_anne/metabolites/data")
	
datname_vec <- vector(); popsig_vec <-  vector(); tempsig_vec <-  vector()
pop_temp_sig_vec  <-  vector() ; intsig_vec <-  vector() ; pop_temp_int_sig_vec <-  vector()
#metabofile=metabofilesare[5]
for ( metabofile in metabofilesare ) {
  tryCatch({
    #-------------- Make metadata file
    overall_metabo_type <- unlist(strsplit( 	rev(unlist(strsplit(metabofile,"/")))[2] 	, "_"))[1]
	if(overall_metabo_type == "NMR"){
		datname <- paste0(overall_metabo_type,"_",file_path_sans_ext(basename(metabofile))) #AQ_CCR
	}else if(overall_metabo_type == "LCMS"){
		datname <- paste0(overall_metabo_type,"_",rev(unlist(strsplit(file_path_sans_ext(basename(metabofile)),"_")))[1]) #AQ_CCR
	}
    dirnameis <- dirname(metabofile) #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten"
    new_dir_name <- paste0(dirnameis,"/",datname,"_results")  #"/faststorage/project/spider2/assay_study_anne/metabolites/data/NMR_Kirsten/AQ_CCR_results"

    dir.create(new_dir_name)
    setwd(new_dir_name)
    
	if(overall_metabo_type == "NMR"){
		datasub <- read.table(metabofile, header=T,sep=",")
		datasub1 <- datasub[,-(5:ncol(datasub))]	;    datasub <- datasub[,-c(3,4)]
		print(datname)	;    dim(datasub)
		write.table(datasub, file=paste0(dirnameis,"/",datname,"_mainDat.txt"),row.names = F)
		write.table(datasub1, file=paste0(dirnameis,"/",datname,"_meta.txt"),row.names = F)
		rm(datasub) 	; rm(datasub1)
		main_dat_file <- paste0(dirnameis,"/",datname,"_mainDat.txt")
		meta_dat_file <-  paste0(dirnameis,"/",datname,"_meta.txt")
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
    mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_overall_", "pdf", 72, width=NA, 5)
    mSet <- PlotPCAScree(mSet, "pca_scree_0_overall_", "pdf", 72, width=NA, 5)
    mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_overall_", "pdf", 72, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPCABiplot(mSet, "pca_biplot_0_overall_", "pdf", 72, width=NA, 1,2)
    #mSet<-Perform.ASCA(mSet, a=1, b=1, x=2, res=2)
    #mSet<-CalculateImpVarCutoff(mSet, spe.thresh = 0.05, lev.thresh = 0.9)

    mSet <- PLSR.Anal(mSet, reg=TRUE)
    #mSet <- PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5) # modify if want right colors
    mSet <- PlotPLS2DScore_mod(mSet, imgName = "pls_pop_score2d_0_overall_", factor_to_use="pop", col_vector=popcols1, "pdf", 72, width=NA, 1,2,0.95,1,0)
		#OBS for LCMS this does not plot points
    mSet <- PlotPLS2DScore_mod(mSet, imgName = "pls_temp_score2d_0_overall_", factor_to_use="temp", col_vector=tempcols, "pdf", 72, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPLS2DScore_mod(mSet, imgName = "pls_pop_score2d_0_overall_", factor_to_use="pop", col_vector=popcols1, "png", 400, width=NA, 1,2,0.95,1,0)
    mSet <- PlotPLS2DScore_mod(mSet, imgName = "pls_temp_score2d_0_overall_", factor_to_use="temp", col_vector=tempcols, "png", 400, width=NA, 1,2,0.95,1,0)
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
		
		#pca_pop_features <- mSet$analSet$pca
		pca_pop_features <- prcomp(Pca_dat_pop_features)
				
		colvec1 = factor_is = substr(mSet$dataSet$cls,1,1)
		col_vector1 = popcols1[as.character(colvec1)] # creates named vector of colors - use as input.
		
		write.csv(data.frame(rownames(Pca_dat_pop_features),factor_is,mSet$dataSet$facB, tempcols[as.character(mSet$dataSet$facB)], col_vector1),file=paste0("Signi_features_effectANDcolor_pop_",datname,".csv"),row.names = F)
		
		pdf(paste0(datname,"_pca_featuers_pop_plot.pdf"),width=5,height=5)
		pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
		dev.off()
		png(paste0(datname,"_pca_featuers_pop_plot.png"),width=5,height=5,units="in",res=400)
		pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
		dev.off()
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
		
		pdf(paste0(datname,"_pca_featuers_temp_plot.pdf"),width=5,height=5)
		pca_plot_metabo(pca_pop_features, factor_is, color_vector=col_vector1)
		dev.off()
		png(paste0(datname,"_pca_featuers_temp_plot.png"),width=5,height=5,units="in",res=400)
		pca_plot_metabo(pca_pop_features, factor_is, col_vector1)
		dev.off()
    }

    print(new_dir_name)
    message("Number of peaks/features:",NROW(mSet$dataSet$prenorm.feat.nms),"\n")
    rm(mSet)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
count_sig_dat <- data.frame(   datname_vec,        popsig_vec,     tempsig_vec,     pop_temp_sig_vec,   intsig_vec,   pop_temp_int_sig_vec)

#write.csv(count_sig_dat,file=paste0(dirnameis,"/Significance_table_asca_metabolites_NMR.csv"),row.names = F)
write.csv(count_sig_dat,file=paste0("/faststorage/project/spider2/assay_study_anne/metabolites/data","/Significance_table_ANOVA_metabolites_ALL.csv"),row.names = F)










####################################
#
###################################
#metabofile <- c("\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_all1.csv")
'metabofile <- c("\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_onlyNamed.csv",
                "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
                "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")
metaDAT_file <- c("\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/metadata_lcms_onlyNamed.csv",
                  "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/metadata_lcms_onlyNamed_CTmax.csv",
                  "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/metadata_lcms_onlyNamed_CCR.csv")

metabofile <- c("\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_onlyNamed_CTmax.csv",
                "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/data_norm_bysum_onlyNamed_CCR.csv")
metaDAT_file <- c("\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/metadata_lcms_onlyNamed_CTmax.csv",
                  "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/metadata_lcms_onlyNamed_CCR.csv")

LCMS_dir <- "\\\\uni.au.dk/Users/au308918/Documents/Spiderlab/assay_study/data/metabolite/LCMS_tobias/"
setwd(LCMS_dir)
names(popcols) <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi") #c("B", "K", "N", "S", "O") #

#-------------- Make metadata file
datname_vec <- vector()
sig_names_list <- list()
popsig_vec <-  vector()
tempsig_vec <-  vector()
pop_temp_sig_vec  <-  vector()
intsig_vec <-  vector()
pop_temp_int_sig_vec <-  vector()
for(filenr in c(1:3) ){
  tryCatch({
    if( grepl("CTmax",metabofile[filenr]) ){datname <- "LCMS_CTmax";  dir.create("CTmax_performance",showWarnings=F); setwd(paste0(LCMS_dir,"CTmax_performance/"))}
    if( grepl("CCR",metabofile[filenr]) ) {datname <- "LCMS_CCR";  dir.create("CCR_performance",showWarnings=F); setwd(paste0(LCMS_dir,"CCR_performance/"))}
    if( grepl("CTmax",metabofile[filenr])==F && grepl("CCR",metabofile[filenr])==F ){ datname <- "LCMS_allPerf"  ;  dir.create("All_performance",showWarnings=F); setwd(paste0(LCMS_dir,"All_performance/"))}

    mSet<-InitDataObjects("pktable", "ts", FALSE) # peak intestisty table, timeseries/two factor, not paired
    mSet<-SetDesignType(mSet, "multi")
    mSet<-Read.TextDataTs(mSet, metabofile[filenr], "rowts")
    mSet<-ReadMetaData(mSet, metaDAT_file[filenr])
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet)
    mSet<-SanityCheckMeta(mSet, 1)
    mSet<-SetDataTypeOfMeta(mSet)
    mSet<-SanityCheckData(mSet)
    message("Sanity check passed, now normalizing.")
    mSet <- FilterVariable(mSet, "none", "F", 25)
    mSet <- PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "SumNorm", "NULL", "ParetoNorm", ratio=FALSE, ratioNum=20)
    # normalised by sum, and pareto scaled
    mSet<-PlotNormSummary(mSet, paste0(datname,"_norm_0_"), "pdf", 72, width=NA)
    mSet<-PlotSampleNormSummary(mSet, paste0(datname,"_snorm_0_"), "pdf", 72, width=NA)
    message("Normalization has been done and plotted. ")
    #meta.vec2 <- [population, treatment]
    mSet<-ANOVA2.Anal(mSet, 0.05, "fdr", "multi", 1, 1)
    mSet<-PlotANOVA2(mSet, paste0(datname,"_aov_abint_"), "pdf", 72, width=NA)
    message("Performed Two-way Anova ")

    # make new dataframe with data (if making ANOVA):
    ANOVA_dat <- read.csv("anova_between_sbj.csv")

    pop_sig <- (which(ANOVA_dat$population.adj.p. <.05))
    temp_sig <- (which(ANOVA_dat$treatment.adj.p. <.05))
    pop_temp_sig <- (intersect(which(ANOVA_dat$treatment.adj.p. <.05), which(ANOVA_dat$population.adj.p. <.05) ))
    int_sig <- (which(ANOVA_dat$Interaction.adj.p. <.05))
    pop_temp_int_sig <- (intersect( intersect(which(ANOVA_dat$treatment.adj.p. <.05), which(ANOVA_dat$population.adj.p. <.05) ), which(ANOVA_dat$Interaction.adj.p. <.05)))

    sig_names_list <- append(sig_names_list, list( ANOVA_dat[pop_sig,1],  ANOVA_dat[temp_sig,1],  ANOVA_dat[pop_temp_sig,1],  ANOVA_dat[int_sig,1]))
    names(sig_names_list)[(NROW(sig_names_list)-3):NROW(sig_names_list)] <- c(paste0(datname,"_Pop_Sign"),  paste0(datname,"_Temp_Sign"),  paste0(datname,"_PopTemp_Sign"),  paste0(datname,"_Interaction_Sign"))
    datname_vec <- append(datname_vec,datname)
    popsig_vec <- append(popsig_vec,NROW(pop_sig))
    tempsig_vec <- append(tempsig_vec,NROW(temp_sig))
    pop_temp_sig_vec <- append(pop_temp_sig_vec,NROW(pop_temp_sig))
    intsig_vec <- append(intsig_vec,NROW(int_sig))
    pop_temp_int_sig_vec <- append(pop_temp_int_sig_vec,NROW(pop_temp_int_sig))

    #PCA of pop-effects
    sig_features_pop <- ANOVA_dat[which(ANOVA_dat$population.adj.p. <.05),1]
    if (NROW(sig_features_pop) > 0){
      Pca_dat_pop_features <- (mSet$dataSet$norm[,which(colnames(mSet$dataSet$norm) %in% sig_features_pop)])

      pop_order <- vector()
      for( each in mSet$dataSet$meta.info$population ){	pop_order <- append(pop_order,grep(each,names(popcols)))	}
      popcols2 <- popcols[pop_order]

      pca_pop_features <- prcomp(Pca_dat_pop_features)
      summa <- summary(pca_pop_features)
      xlimsare <- c(-max(abs(c(pca_pop_features$x[,1],pca_pop_features$x[,2]))),max(abs(c(pca_pop_features$x[,1],pca_pop_features$x[,2]))))

      pdf(paste0(datname,"_pca_featuers_pop_plot.pdf"),width=5,height=5)
      par(mfrow=c(1,2), mar=c(4.7,4.2,4,2.5), cex=.7)
      #p1
      plot(pca_pop_features$x[,1], pca_pop_features$x[,2], col=popcols2, pch=1,
           xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"),
           cex=0.8, xlim=xlimsare,ylim=xlimsare, main=datname)
      ordiellipse(pca_pop_features,popcols2,conf=0.8, col=levels(as.factor(popcols2)))
      #p2
      plot(pca_pop_features$x[,2], pca_pop_features$x[,3], col=popcols2, pch=1,
           xlab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"),
           cex=0.8, xlim=xlimsare,ylim=xlimsare, main=datname)
      ordiellipse(pca_pop_features,popcols2,conf=0.8, col=levels(as.factor(popcols2)))
      #scree
      screeplot(pca_pop_features, type="l", main=paste0("Screeplot ",datname))
      #bi
      biplot(pca_pop_features, main=paste0("Biplot ",datname))
      dev.off()

      png(paste0(datname,"_pca_featuers_pop_plot.png"),width=5,height=5,units="in",res=400)
      par(mfrow=c(1,2), mar=c(4.7,4.2,4,2.5), cex=.7)
      #p1
      plot(pca_pop_features$x[,1], pca_pop_features$x[,2], col=popcols2, pch=1,
           xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"),
           cex=0.8, xlim=xlimsare,ylim=xlimsare, main=datname)
      ordiellipse(pca_pop_features,popcols2,conf=0.8, col=levels(as.factor(popcols2)))
      #p2
      plot(pca_pop_features$x[,2], pca_pop_features$x[,3], col=popcols2, pch=1,
           xlab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"),
           cex=0.8, xlim=xlimsare,ylim=xlimsare, main=datname)
      ordiellipse(pca_pop_features,popcols2,conf=0.8, col=levels(as.factor(popcols2)))
      #scree
      screeplot(pca_pop_features, type="l", main=paste0("Screeplot ",datname))
      #bi
      biplot(pca_pop_features, main=paste0("Biplot ",datname))
      dev.off()
    }

    rm(mSet)
    setwd(LCMS_dir)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
  setwd(LCMS_dir)

  count_sig_dat <- data.frame(   datname_vec,        popsig_vec,     tempsig_vec,     pop_temp_sig_vec,   intsig_vec,   pop_temp_int_sig_vec)
write.csv(count_sig_dat,file="LCMS_Significance_table_ANOVA_metabolites_LCMS.csv",row.names = F)

  # write feature names (Feature LCSM IDs)
sig_names_list <- sig_names_list[-c(1:4)]
write.xlsx(sig_names_list, file = "LCMS_significant_feature_name.xlsx")


  # take the names of the significant ones and grep their real name from file.
feature_info_dat <- read.csv("data_norm_bysum_onlyNamed_incl_feature_info.csv", header=T, stringsAsFactors = T)


merge_func <- function(x){
  merged_df <- merge(x,feature_info_dat,by.x=1,by.y=4, all.x=T, all.y=F)
  return(merged_df)
}
extract_two_cols <- function(x){
  extracted_cols <- x[,1:2]
  return(extracted_cols)
}

sig_names_list_large <- lapply(sig_names_list,merge_func)
sig_names_list_ID <- lapply(sig_names_list_large,extract_two_cols)

write.xlsx(sig_names_list_ID, file = "LCMS_significant_feature_ID.xlsx")'
