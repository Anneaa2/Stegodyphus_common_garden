#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript

#srun --mem=16g --pty /bin/bash
#conda activate sources_of_var; R
#source("figures_for_paper/make_data_lists.r")

#setwd("/faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/")


t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
rgb.val <- col2rgb(color)	## Get RGB values for named color
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name) ## Make new color using input color as base and alpha set by transparency
invisible(t.col) ## Save the color
}

# other script data:


# Datafiles:



## Get data into lists


## MAke distributions of correlation coefficients
		#make lists for data
listofdists <- list() #Make list to add correlation vectors to.
listofp <- list() #Make list to add correlation vectors to.
list_overlap <- list()
list_env <- list()

data_to_list <- function(nameofdata_is, 
						#files_are_snps1=files_are_snps, 
						files_are_popstruck1=files_are_popstruck#, 
						#files_are_cisvar1=files_are_cisvar, 
						#files_are_micro1=files_are_micro, 
						#files_are_micro_genus1=files_are_micro_genus, 
						#files_are_mdrm_snp1=files_are_mdrm_snp,
						#files_are_mdrm1=files_are_mdrm, 
						#files_are_mdrm_micro1=files_are_mdrm_micro, 
						#files_are_mdrm_micro_genus1=files_are_mdrm_micro_genus
						)
						{
	if ( analysis_type== "mantel"){
		files_all<-c(files_are_popstruck1)
		# files_all<-files_all[c(1,2,3,4,5)]
		
		# files_all<-c(files_are_snps1, files_are_popstruck1, files_are_cisvar1, files_are_micro1,files_are_micro_genus1)
		# files_all<-files_all[c(1,2,3,4,5)] #files_all_cov10<-files_all[c(1,2,4,8)] 
	}# else if (analysis_type== "mrm"){
	#	files_all<-c(files_are_mdrm_snp1, files_are_mdrm1, files_are_mdrm_micro1,	files_are_mdrm_micro_genus1 )
	#	files_all<-files_all[c(1,2,3,4)]  #files_all_cov10<-files_all[c(1,3,6)] 
	#}
	print(files_all)
	counter=0
	for (file in files_all){
		nameofdatais <- nameofdata_is #"RealCor_COR"  "RandCor_COR"
		N_tails = 3
		#if ( grepl("NOCORRECT",file)  ) { N_tails = N_tails+1 }
		#if ( grepl("PERMrand",file) ) {  N_tails = N_tails+1  }
		
		#if ( grepl("micro",file) &  grepl("mrm",file)  ) {  N_tails = N_tails+1  }
		#if ( grepl("micro",file) &  grepl("mantel",file) ) {  N_tails = N_tails+0  }
			# read in file
		data <- read.csv(file, header=T,stringsAsFactors=T, row.names=1)
			#correct for analysis type of input
		if ( grepl("micro",file) ) {	#sapply(strsplit(file,"[_.]"),"[[",7) == "micro" ){ # if micro
		#	for ( typeofis in c("relative","absolute") ) {
		#		counter=counter+1
		#		for ( axisis in unique(data$env_vector)){
					
		#			data2 <- data[which(data$env_vector == as.character(axisis)),]	
					#print(counter)
		#			print(typeofis)
		#			data1 <- data2[grep(typeofis, data2$micro_vector),]
					
		#			if ( grepl("mrm",file) ) { #sapply(strsplit(file,"[_.]"),"[[",6)=="mrm"  || sapply(strsplit(file,"[_.]"),"[[",5)=="mrm" ){
		#				data_name1 <- sapply(strsplit(file,"_|\\."),tail,N_tails)[1,]
		#				data_name <- paste0(data_name1,"_",typeofis)
		#				datanameis <- paste0(analysis_type,"_",nameofdatais,"_",data_name,"_",sapply(strsplit(as.character(axisis),"_",),"[[",3))
		#				listofdists[[datanameis]] <<- data1$rsquared
		#				listofp[[datanameis]] <<- data1$pval
		#				list_overlap[[datanameis]] <<- data1[,1]
		#				list_env[[datanameis]] <<- data1[,2]
		#			} else{
		#				data_name1 <- sapply(strsplit(file,"_|\\."),tail,N_tails)[1,]
		#				data_name <- paste0(data_name1,"_",typeofis)
		#				datanameis <- paste0(analysis_type,"_",nameofdatais,"_",data_name,"_",sapply(strsplit(as.character(axisis),"_",),"[[",3))
		#				listofdists[[datanameis]] <<- data1$r_vector # assign to environment variable by <<-
		#				listofp[[datanameis]] <<- data1$sig_vector
		#				list_overlap[[datanameis]] <<- data1[,1]
		#				list_env[[datanameis]] <<- data1[,2]
		#		}	}
		#		print(datanameis)
		#	}
		}else{ #if smp or snp
			counter=counter+1
			print(counter)
			data_name <- paste0(	sapply(strsplit(file,"_|\\."),tail,N_tails)[1:2,] , collapse="_")
			for ( axisis in unique(data$env_vector)){
				data2 <- data[which(data$env_vector == as.character(axisis)),]
				#if ( sapply(strsplit(file,"[_.]"),"[[",6)=="mrm"  || sapply(strsplit(file,"[_.]"),"[[",5)=="mrm" ){
				#	datanameis <- paste0(analysis_type,"_",nameofdatais,"_",data_name,"_",sapply(strsplit(as.character(axisis),"_",),"[[",3))
				#	listofdists[[datanameis]] <<- data2$rsquared
				#	listofp[[datanameis]] <<- data2$pval
				#	list_overlap[[datanameis]] <<- data2[,1]
				#	list_env[[datanameis]] <<- data2[,2]
				#} else{
					datanameis <- paste0(analysis_type,"_",nameofdatais,"_",data_name,"_",sapply(strsplit(as.character(axisis),"_",),"[[",3))
					listofdists[[datanameis]] <<- data2$r_vector
					listofp[[datanameis]] <<- data2$sig_vector
					list_overlap[[datanameis]] <<- data2[,1]
					list_env[[datanameis]] <<- data2[,2]
				#}	
			}
			print(datanameis)
}	} 	}





 # Real data 
if ( analysis_type== "mantel"){
	#files_are_popstruck <- list.files(pattern="pca_sunimp2_partial_mantel_res_.*_smps_refpop.csv",recursive=T,path="./genewise_SMPS_res10/", full.names=T) #4 files
	files_are_popstruck <- list.files(	path="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/",
										pattern=glob2rx("PCA_partial_mantel_res_genes_smps_refpop_d*.csv"),
										recursive=F, full.names=T) #4 files

	data_to_list("RealCor_COR")
} #else if (analysis_type== "mrm"){
#	files_are_mdrm <- list.files(pattern="pca_sunimp2_mrm_res_genes_smp.csv",recursive=T,path="./genewise_SMPS_res10/", full.names=T) #4 files
#	data_to_list("RealCor_COR")
#}

## Works

 # Permuted data
if ( analysis_type== "mantel"){
	#files_are_popstruck <- list.files(pattern="pca_sunimp2_partial_mantel_res_.*_smps_refpop_PERMrand.csv",recursive=T,path="./genewise_SMPS_res10/", full.names=T) #4 files
	files_are_popstruck <- list.files(	path="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results",
										pattern=glob2rx("PCA_partial_mantel_res_genes_smps_refpop_PERMrand_d*.csv"),
										recursive=F, full.names=T) #4 files
	
	data_to_list("RandCor_COR")
} #else if (analysis_type== "mrm"){
#	files_are_mdrm <- list.files(pattern="pca_sunimp2_mrm_res_genes_smp_PERMrand.csv",recursive=T,path="./genewise_SMPS_res10/", full.names=T) #4 files
#	data_to_list("RandCor_COR")
#}




