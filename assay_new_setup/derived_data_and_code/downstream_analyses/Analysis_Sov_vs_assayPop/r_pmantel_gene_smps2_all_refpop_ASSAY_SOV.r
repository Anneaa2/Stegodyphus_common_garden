#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript


###################							###################
#######			This is for SMPS - methylation			#######
#######		Partial Mantel Tests on distance matrices 	#######
#################							###################

# set this to run for all SMP files


#install.packages('vegan', repos='http://cran.us.r-project.org')
library(vegan) #could also be ecodist
library(usedist)
library(Matrix)

	# Parse arguments
argsis<-commandArgs(trailingOnly=TRUE)
directory_is=argsis[1]	
	# "/faststorage/project/spider2/AAL_20190722_SOV_correlations/"
	# directory_is=
	# NOT USED?
env_dir=argsis[2]	
	#  "/faststorage/project/spider2/AAL_20190722_SOV_correlations/environmental_dist_matrixes/reduce_env_jan20/pca_red/"
	# env_dir="../AAL_20190722_SOV_correlations/environmental_dist_matrixes/reduce_env_jan20/pca_red/"
	# env_dir="../AAL_20190722_SOV_correlations/environmental_dist_matrixes/reduce_env_june20_imp/"
runname=argsis[3]	
	# "pca1"
	# runname="PCA"

SMPS_GENE_FILE=argsis[4]	# this one changes
	# SOV - ORIG: SMPS_GENE_FILE="/faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/gene_pairwise_fst_smps_d5D0.2perc.txt"
	# SMPS_GENE_FILE="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/WML_files/gene_pairwise_fst_smps_d10D30_SOV_wild.txt"

snp_all_file=argsis[5]		
	# "/faststorage/project/spider2/AAL_20190722_SOV_correlations/january_run/snp_smp_data/fst_snp/snp_genes_all_simple.fst"
	# snp_all_file="/faststorage/project/spider2/AAL_20190722_SOV_correlations/january_run/snp_smp_data/fst_snp/snp_genes_all_simple.fst"
outfolder=argsis[6]		
	#  "/faststorage/project/spider2/AAL_20190722_SOV_correlations//correlation_results//genewise_SMPS_res/"
	#  outfolder="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/env_correlation_results/"

message("Directory is: ", directory_is)
message("Folder with environmental distances is: ", env_dir)
message("The name of this run is: ", runname)
message("The file of SMPs in genes is: ", SMPS_GENE_FILE)
message("The file with pop struct to control for, is: ", snp_all_file)
message("Outfolder is: ", outfolder)

############################################################################################################
####################		 Distance matrixes for env variables "....dist.csv"			####################
############################################################################################################

#env_dir=paste0(directory_is,"/environmental_dist_matrixes/reduce_env_2/dist_matx_18/")
	## List env_var files:
#env_files<- list.files(env_dir,pattern = "_dist.csv",full.names = T,recursive = F) # 17 files
env_files<- list.files(env_dir,pattern = "_DIST_",full.names = T,recursive = F)
	#read env file # 16 files 


################################## read pop structure file##########################################################3
## File with all snps - not just the one passing threshold
##################################

			#(not just the ones passing threshold)

#snp_all_file="/faststorage/project/spider2/dumicola_genome/AAL_20190606_SOV_snps_filtering/popoolation_run/snp_fst_snps/fst_perbase_snpmean_simply.fst"
		# OBS ETOSHA BECOMES OTAWI
		# snp_nr	1:2		1:3		1:4		1:5		1:6		2:3		2:4		2:5		2:6		3:4		3:5		3:6		4:5		4:6		5:6
		#	"xx"		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst		fst
		# sync files order :	"A",		"B",		"C",		"D",			"O",			"N" 			(from lins sync files - see /faststorage/project/spider2/dumicola_genome/Lin_20180103_mapping_v2/run_aln.sh)
		# 						A=gobabis, 	B=Etosha,	C=Betta, 	D=Karasburg, 	E=Stampriet, 	N=Ndumo			(From sandbjerg slides)
		#									(OTAWI)
		
		# snp_nr	1:2			1:3			1:4			1:5			1:6			2:3			2:4			2:5			2:6			3:4			3:5			3:6			4:5			4:6			5:6
		#	1		fst			fst			fst			fst			fst			fst			fst			fst			fst			fst			fst			fst			fst			fst			fst
		#			G:E[2]		G:B[3]		G:K[4]		G:S[5]		G:N[6]		E:B[7]		E:K[8]		E:S[9]		E:N[10]		B:K[11]		B:S[12]		B:N[13]		K:S[14]		K:N[15]		S:N[16]
		

tempfile_snps_all<-read.csv(snp_all_file,header = F,sep="\t",stringsAsFactors = F,strip.white = F)
		# mean each column and construct the distance file as a matrix yourself
nrow(tempfile_snps_all)
#		#36999
mean_fst_per_pair_all<-colMeans(tempfile_snps_all[,-1])
		#V1           V2           V3           V4           V5           V6	         V7           V8           V9          V10          V11          V12			V13          V14          V15          V16
		#1.311481e+06 5.842154e-02 7.418002e-02 6.087309e-02 4.598772e-02 1.188041e-01	6.016511e-02 9.695195e-02 9.146080e-02 1.305493e-01 1.184908e-01 1.139207e-01	1.459124e-01 4.401869e-02 1.319756e-01 1.254495e-01
		
		
		#create distance matrix
pop_struct_dist1=matrix(0,6,6)
pop_struct_dist1[lower.tri(pop_struct_dist1,diag=F)]<-mean_fst_per_pair_all
pop_struct_dist1[upper.tri(pop_struct_dist1)]<-t(pop_struct_dist1)[upper.tri(pop_struct_dist1)]

		#enter names:
names_vec=c("G","O","B","K","S","N")
colnames(pop_struct_dist1)<-names_vec
rownames(pop_struct_dist1)<-names_vec
#
		# Get them in the right order: --> Env files order: B K N S G E
	#pop_struct_dist<-as.matrix(dist_subset(pop_struct_dist1,c("B","K","N","S","G","O")))
	# Order them to get similar order.
pop_struct_dist<-as.matrix(dist_subset(pop_struct_dist1,sort(c("B","K","N","S","G","O"))))
		# It has already been exported in snp-script
# remove rows and columns with G and N
pop_struct_dist <- pop_struct_dist[ grep("G|N",rownames(pop_struct_dist),invert=T), grep("G|N",colnames(pop_struct_dist),invert=T) ]

write.csv(file=paste0(outfolder,"/distance_matrices/",runname,"_popStructure_fst_matrix_all_genes_snps.csv"),pop_struct_dist)
	


################################### read each snp from GENE file: ##################################################
## File with all sMps per gene
###################################


#Methylation naming:	 "gene"		"PopA"		"PopB"		"PopC"		"PopD"		"PopE"		"PopN"
# Methylation number:					1			2			3			4			5			6
# Env_naming:						Gobabis		Etosha		Betta		Karasburg	Stampriet	Ndumo

	### Translation of numbers: 1=Gobabis, 2=Etosha, 3=Betta, 4=Karasburg, 5=Stampriet, 6=Ndumo
		#		B		K			N			S			G			E
		#B		0	snp[11]		snp[13]		snp[12]		snp[3]		snp[7]
		#K	snp[11]		0		snp[15]		snp[14]		snp[4]		snp[8]
		#N	snp[13]		snp[15]		0		snp[16]		snp[6]		snp[10]
		#S	snp[12]		snp[14]		snp[16]		0		snp[5]		snp[9]
		#G	snp[3]		snp[4]		snp[6]		snp[5]		0		snp[2]
		#E	snp[7]		snp[8]		snp[10]		snp[9]		snp[2]		0
		
	#	1		2		3		4		5		6		7		8		9		10		11		12		13		14		15		16
	# 	g1		1:2     1:3     1:4     1:5     1:6     2:3     2:4     2:5     2:6     3:4     3:5     3:6     4:5     4:6     5:6
	
	
	# SMPS_GENE_FILE="/faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/gene_pairwise_fst_smps_d10D30.txt"
tempfile_smps<-read.csv(SMPS_GENE_FILE,header = T,sep="\t",stringsAsFactors = F,strip.white = F)

#tempfile_smps<-cbind(rownames(tempfile_smps),tempfile_smps)

	# Subsetting to not run genes with no sMps:
tempfile_smps_sub<-tempfile_smps[which(rowSums(tempfile_smps[,2:ncol(tempfile_smps)]) != 0),]
message("Number of genes to run: before removing genes with no snps: ",nrow(tempfile_smps)," and after removal: ", nrow(tempfile_smps_sub))
'
Number of genes to run: before removing genes with no snps: 32562 and after removal: 32236
'
	


#calculate mean FST for pops.
mean_fst_SMPS<-colMeans(tempfile_smps[,-1])

SMP_dist1 <- matrix(0,4,4)
SMP_dist1[lower.tri(SMP_dist1,diag=F)] <- mean_fst_SMPS
SMP_dist1[upper.tri(SMP_dist1)] <- t(SMP_dist1)[upper.tri(SMP_dist1)]

names_vec <- sort(c("O","B","K","S"))
colnames(SMP_dist1) <- names_vec
rownames(SMP_dist1) <- names_vec
	# Get them in the right order: --> Env files order: B K N S G E
# SMP_dist<-as.matrix(dist_subset(SMP_dist1,c("B","K","N","S","G","O")))
SMP_dist <- SMP_dist1
#write

datatype_name <- gsub(".txt","",			gsub("^.*smps_", "", basename(SMPS_GENE_FILE))	)
write.csv(file=paste0(outfolder,"/distance_matrices/",runname,"_mean_fst_matrix_all_genes_smps_",datatype_name,".csv"),SMP_dist)





##########################################################
# Read in SNPS file#######################################
##########################################################

#tempfile_sNps<-read.csv(snp_gene_file,header = F,sep="\t",stringsAsFactors = F,strip.white = F)
#rownames(tempfile_sNps)<-tempfile_sNps[,1]; tempfile_sNps<-tempfile_sNps[,-1]


	###############################					##################################
	###########################	 Running Mantel tests 	##############################
	############################### 				##################################

dir.create( paste0(outfolder,"/distance_matrices/", gsub("^d[0-9]{2}D[0-9]{2}_", "", datatype_name)) )	
	
r_vector=vector() ; sig_vector=vector() ; perm_vector=vector() ; smp_vector=vector() ; env_vector=vector()
r_vector1=vector() ; sig_vector1=vector() ; perm_vector1=vector() ; smp_vector1=vector() ; env_vector1=vector()
r_vector2=vector() ; sig_vector2=vector() ; perm_vector2=vector() ; smp_vector2=vector() ; env_vector2=vector()
iter=0

for (file in env_files) {
	iter=iter+1
	tempfile_env<-read.csv(file,row.names=1,header = T,stringsAsFactors = F,strip.white = T)
	predictor_var=as.matrix(tempfile_env)

	# subset predictor/env variable to 4 pops (and order them):
	predictor_var <- as.matrix(dist_subset(predictor_var, sort(colnames(predictor_var))))
	predictor_var <- predictor_var[ grep("G|N",rownames(predictor_var),invert=T), grep("G|N",colnames(predictor_var),invert=T) ]

		########################### each snp correlated
	for (line in 1:nrow(tempfile_smps_sub)) {
		gene_nameis<-tempfile_smps_sub[line,1]
		smp<-as.numeric(as.matrix(tempfile_smps_sub[line,-1]))
				
			# make correct structure of snp_line
		dependent_var1=matrix(0,4,4)
		dependent_var1[lower.tri(dependent_var1,diag=F)]<-smp#[2:NROW(smp)]
		dependent_var1[upper.tri(dependent_var1)]<-t(dependent_var1)[upper.tri(dependent_var1)] 
			#enter names:
		names_vec=sort(c("O","B","K","S"))
		colnames(dependent_var1)<-names_vec ; rownames(dependent_var1)<-names_vec
			# Get them in the right order: Env files order: B K N S G E
		#dependent_var<-as.matrix(dist_subset(dependent_var1,c("B","K","N","S","G","O")))
			# already ordered everythin alphabetically, including env_vars (predictor)
		dependent_var <- dependent_var1
		if (iter==1){
			write.csv(file=paste0(outfolder,"/distance_matrices/", gsub("^d[0-9]{2}D[0-9]{2}_", "", datatype_name),"/", runname,"_gene_dist_smp_", datatype_name, "_", gene_nameis,".csv"),dependent_var)
		}
		
		man_result<-mantel.partial(dependent_var,predictor_var,pop_struct_dist,method="pearson",permutations = 999)
		
		#if (is.na(man_result$statistic) || is.na(man_result$signif)){
		#	message("r or p-value is NA.")
		#	print(man_result)
		#	print(dependent_var)
		#	print(predictor_var)
		#	print(pop_struct_gene)
		#}
			#smp number to traceback
		smp_vector<-append(smp_vector,paste0("gene_smp_",tempfile_smps_sub[line,1]))
		r_vector<-append(r_vector,man_result$statistic) #r statistic
		sig_vector<-append(sig_vector,man_result$signif) 
		perm_vector<-append(perm_vector,man_result$permutations) 
			#environmental variable correlated with
		env_vector<-append(env_vector,sapply(strsplit(basename(file),"_dist"),"[",1))
		
		#alternative test 1:random permutation
		predictor_var1 <- predictor_var; predictor_var1[lower.tri(predictor_var1)] <- sample(as.dist(predictor_var)) ; predictor_var1 <- forceSymmetric(predictor_var1,uplo="L")
		man_result1<-mantel.partial(dependent_var,predictor_var1,pop_struct_dist,method="pearson",permutations = 999)
		
		smp_vector1<-append(smp_vector1,paste0("gene_smp_",tempfile_smps_sub[line,1]))
		r_vector1<-append(r_vector1,man_result1$statistic) #r statistic
		sig_vector1<-append(sig_vector1,man_result1$signif) #r statistic
		perm_vector1<-append(perm_vector1,man_result1$permutations) #r statistic
			#environmental variable correlated with
		env_vector1<-append(env_vector1,sapply(strsplit(basename(file),"_dist"),"[",1))
		
		#alternative test 2:random permutation once
		#if( line == 1 ){
		#	predictor_var2 <- predictor_var; predictor_var2[lower.tri(predictor_var2)] <- sample(as.dist(predictor_var)) ; predictor_var2 <- forceSymmetric(predictor_var2,uplo="L")
		#}
		#man_result2<-mantel.partial(dependent_var,predictor_var2,pop_struct_dist,method="pearson",permutations = 999)
			
		#smp_vector2<-append(smp_vector2,paste0("gene_smp_",tempfile_smps_sub[line,1]))
		#r_vector2<-append(r_vector2,man_result2$statistic) #r statistic
		#sig_vector2<-append(sig_vector2,man_result2$signif) #r statistic
		#perm_vector2<-append(perm_vector2,man_result2$permutations) #r statistic
			#environmental variable correlated with
		#env_vector2<-append(env_vector2,sapply(strsplit(basename(file),"_dist"),"[",1))
		
		message("Finished ",line/nrow(tempfile_smps_sub)*100," % of current environmental variable: ", iter," of ",NROW(env_files),". Var is: ",basename(file))
	}
}

man_result_df<-data.frame(smp_vector,env_vector,r_vector,sig_vector,perm_vector)
man_result_df1<-data.frame(smp_vector1,env_vector1,r_vector1,sig_vector1,perm_vector1)
#man_result_df2<-data.frame(smp_vector2,env_vector2,r_vector2,sig_vector2,perm_vector2)
message("results stored to dataframe")

	#write.csv(file=paste0(directory_is,"/genewise_smps_env/results/partial_mantel_res_genes_smps1_all.csv"),man_result_df)
write.csv(file=paste0(outfolder,"/", runname,"_partial_mantel_res_genes_smps_refpop_", datatype_name,".csv"),man_result_df)
write.csv(file=paste0(outfolder,"/", runname,"_partial_mantel_res_genes_smps_refpop_PERMrand_", datatype_name,".csv"),man_result_df1)
#write.csv(file=paste0(outfolder,"/",runname,"_partial_mantel_res_genes_smps_refpop_PERMonce.csv"),man_result_df2)
message(paste0("outfile is here: ", outfolder,"/", runname,"_partial_mantel_res_genes_smps_refpop_", datatype_name,".csv"))
				#/faststorage/project/spider2/dumicola_genome/AAL_20190722_SOV_correlations/genewise_smps_env/			
	# data: snp_nr	environmental_parameter		r_statistics	Significance	nr_permutations
#paste0(outfolder,"/", runname,"_partial_mantel_res_genes_smps_refpop_", datatype_name,".csv")

	###########################													##############################
	###########################	 Filter on r threshold and/or significance? 	##############################
	########################### 												##############################
'
thressie=read.csv(paste0(directory_is,"/thresholds/",runname,"_mantel_r_threshold_99quant.csv"), header=T,stringsAsFactors=F,strip.white=F,row.names=1)
threshold=thressie[1,1]

	# Keep the ones above threshold:
man_result_df_passthreshold<-man_result_df[which(man_result_df$r_vector > threshold),]

write.csv(file=paste0(outfolder,"/",runname,"_partial_mantel_res_genes_smps_passthresh_refpop.csv"),man_result_df_passthreshold)

	# Keep the ones that pass significance
man_result_df_passthreshold_sign <-man_result_df_passthreshold[which(man_result_df_passthreshold$sig_vector < 0.05),]

write.csv(file=paste0(outfolder,"/",runname,"_partial_mantel_res_genes_smps_passthresh_sign_refpop.csv"),man_result_df_passthreshold_sign)

'

# To get criteria for R:
# Run part.mantel (fx. 20000) with pairwise fst values simulated (numbers between 1 and 0) and out env. factors (real values), controlling for real pop structure (genome wide fst).
# use fx 99% quantile of the 100000 r values as threshold. (20000 x 5 env. factors in fishers case)
# if a snp is above this threshold, it is considered to be associated with the env. var.

#Assumes linear relation between env. factors and allele frequencies







