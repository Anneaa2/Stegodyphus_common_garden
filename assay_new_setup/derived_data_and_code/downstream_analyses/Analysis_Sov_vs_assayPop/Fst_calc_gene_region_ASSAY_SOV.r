#!/home/anneaa/miniconda3/envs/sources_of_var/bin/Rscript


#=========== Get genes where methylation patterns are found in both SOV and Assay

# Modified to analyse both Assay and SOV subset data. Original script file: 
	# /home/anneaa/spider2/faststorage/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/GO_get_gene_region_sunrem.sh


# THis script calculates pairwise Fst for Assay populations per temperature (5 pop series), and SOV populations (B,O,K,S)
# The pairwise FST will be setup as a distance matrix, and correlated to the climate PCA from SOV study.
	# The Genes being found in SOV, but NOT in Assay --> Must be genes that can change methylation level
	# The genes found in both SOV and Assay --> Must be stable genes, that does not change with changed environment.




#-----------Coverage of methylated and unmethylated CpG's.
# BASH
# To use for calculating methylation level difference per gene among populations.
# Using Jesper's measurement (methylation coverage ratio).
# Gene, coverage of methy, coverage of unmethy, count of cpg sites.


#for file in `ls gene_region/Pop*d*D30*.gene.cpg`
#do
#	sort -s -k 13,13 ${file} | awk '
#	 BEGIN{FS="\t";OFS="\t";m=0;u=0;count=0;getline;gene=$13;m=m+$5;u=u+$6;count++}
#	 $13!=gene{print gene,m,u,count;m=0;u=0;gene=$13;count=0;m=m+$5;u=u+$6;count++;next}
#	 {m=m+$5;u=u+$6;count++}
#	 END{print gene,m,u,count;m=0;u=0;count=0}' > ${file/.cpg/.cpg.sum}
#done

	# m = methylated , u = unmethylated

# Count number of Cs in each gene also (count variable)



#----------------- Calculate methylation coverage ratio:
# (R).

#conda activate sources_of_var; R

# function for calculating Fst
	## Fst pairwise
	frq2fst<-function(freq){
		#message("freq is:", freq)
		he<-2*freq-2*(freq^2)
		#message("he is:",he)
		hs<-rowMeans(he)
		#message("hs is:",hs)
		mean.freq<-rowMeans(freq)
		ht<-2*mean.freq-2*(mean.freq^2)
		#message("ht is:",ht)
		fst<-1-hs/ht
		fst[is.na(fst)]<-0
		#message("fst is:",fst)
		fst
	}


	# THESE FILES already exist for both Assay and SOV
		# assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/
			# patternis="*15_d10D32*gene.cpg.sum"
			# patternis="*19_d10D32*gene.cpg.sum"
			# patternis="*23_d10D32*gene.cpg.sum"
			# patternis="*25_d10D32*gene.cpg.sum"
			# patternis="*29_d10D32*gene.cpg.sum"
		# ../Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/
			# patternis="d10D30*gene.cpg.sum"
		# sync files order :	"A",		"B",		"C",		"D",			"O",			"N" 			(from lins sync files - see /faststorage/project/spider2/dumicola_genome/Lin_20180103_mapping_v2/run_aln.sh)
		# 		A=gobabis, 	B=Etosha,	C=Betta, 	D=Karasburg, 	E=Stampriet, 	N=Ndumo			(From sandbjerg slides)
pattern_list=c("*15_d10D32*gene.cpg.sum", "*19_d10D32*gene.cpg.sum", "*23_d10D32*gene.cpg.sum",
			"*25_d10D32*gene.cpg.sum", "*29_d10D32*gene.cpg.sum", "*d10D30*gene.cpg.sum")
for( patternis in pattern_list ){
	# patternis="*d10D3*gene.cpg.sum"
	# patternis="*d10D30*gene.cpg.sum"
	# patternis="*29_d10D32*gene.cpg.sum"

	files <- list.files("assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG",pattern=glob2rx(patternis),full.names=T)
	files <- append(files, list.files("../Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region",pattern=glob2rx(patternis),full.names=T))
	if( NROW(grep("SOV",files)) > 0 ){
		files <- files[ grep("PopA|PopN",files,invert=T) ]	
	}

	data<-read.table(files[1],stringsAsFactors=F)
	c_inreg<-data.frame(data[,1],data[,4],stringsAsFactors=F)
	region_methyl_reads<-data.frame(data[,1],data[,2],stringsAsFactors=F)
	region_total_reads<-data.frame(data[,1],data[,2]+data[,3],stringsAsFactors=F)
	data<-data.frame(data[,1],data[,2]/(data[,2]+data[,3]),stringsAsFactors=F)
	for(i in 2:length(files))
	{
		bdata<-read.table(files[i],stringsAsFactors=F)
		bc_inreg<-data.frame(bdata[,1],bdata[,4],stringsAsFactors=F)
		B_region_methyl_reads<-data.frame(bdata[,1],bdata[,2],stringsAsFactors=F)
		B_region_total_reads<-data.frame(bdata[,1],bdata[,2]+bdata[,3],stringsAsFactors=F)
		bdata<-data.frame(bdata[,1],bdata[,2]/(bdata[,2]+bdata[,3]),stringsAsFactors=F)
		data<-merge(data,bdata,by.x=1,by.y=1,sort=F)
		c_inreg<-merge(c_inreg,bc_inreg,by.x=1,by.y=1,sort=F)
		region_methyl_reads<-merge(region_methyl_reads,B_region_methyl_reads,by.x=1,by.y=1,sort=F)
		region_total_reads<-merge(region_total_reads,B_region_total_reads,by.x=1,by.y=1,sort=F)
	}
	if( NROW(grep("SOV",files)) > 0 ){
		names(data)<-c("Gene","PopO","PopB","PopK","PopS")	# files are called B,C,D,E. Renaming them in data.
	}else{
		names(data)<-c("Gene","PopB","PopK","PopO","PopS")	# file order (O was formerly named N)
	}
	names(c_inreg) <- names(data)
	names(region_methyl_reads) <- names(data)
	names(region_total_reads) <- names(data)
	data <- data[,c(1,order(names(data[2:5]))+1)]
	# 					B=Etosha/Otavi,	C=Betta, 	D=Karasburg, 	E=Stampriet,	(From sandbjerg slides)
	c_inreg <- c_inreg[,c(1,order(names(c_inreg[2:5]))+1)]
	region_methyl_reads <- region_methyl_reads[,c(1,order(names(region_methyl_reads[2:5]))+1)]
	region_total_reads <- region_total_reads[,c(1,order(names(region_total_reads[2:5]))+1)]

	c_inreg$Mean <- rowMeans(c_inreg[,-1])

	#data10<-data
	#region_names_vector10 <- data[,1]
	#papnamesare10 <- names(data)[-1]
	if( NROW(grep("SOV",files[1])) > 0 ){
		coverage_thesh <- sapply(strsplit(files[1],"/|\\.|_|B|N|S|K"),"[[",16)
		name_file <- paste0(coverage_thesh,"_SOV_wild")
	}else{
		temp_series <- sapply(strsplit(files[1],"/|\\.|_|B|N|S|K"),"[[",14)
		coverage_thesh <- sapply(strsplit(files[1],"/|\\.|_|B|N|S|K"),"[[",15)
		name_file <- paste0(coverage_thesh,"_ASSAY_",temp_series)
	}
	outfolder = "assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/WML_files"
	write.table(data, file=paste0(outfolder,"/all_pops_mcr_gene_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	write.table(c_inreg, file=paste0(outfolder,"/all_pops_nrCpg_gene_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	write.table(region_methyl_reads, file=paste0(outfolder,"/all_pops_methreads_gene_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	write.table(region_total_reads, file=paste0(outfolder,"/all_pops_totalreads_gene_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

	nrow(data)
	# [1] 31948
	NROW(which(rowSums(data[,-1])==0))
	#[1] 118


	# Calculate fst?
	# use dist?
	a_matrix<-data[,-1]
	pairwise_fst_methylation<- data.frame(matrix(NA,nrow=nrow(a_matrix)))[-1]
	for (popnum in c(1:3)){
		print(popnum)
		for (poppair in popnum+1:3){
			if ( poppair > 4){break}
			print(poppair)
			fst_calc <- frq2fst(a_matrix[,c(popnum,poppair)])
			pairwise_fst_methylation<-cbind(pairwise_fst_methylation,fst_calc)
			#print(ncol(pairwise_fst_methylation))
			popname<-sapply(strsplit(names(a_matrix),"Pop"),"[[",2)[popnum]; popname2<-sapply(strsplit(names(a_matrix),"Pop"),"[[",2)[poppair]
			colnames(pairwise_fst_methylation)[ncol(pairwise_fst_methylation)]<- paste0(popnum,poppair,"_fst_",popname,popname2)
		}
	}

	#name_file=sapply(strsplit(patternis,".cov"),"[[",1)
	name_file
	pairwise_fst_methylation <- cbind(data$Gene,pairwise_fst_methylation)
	colnames(pairwise_fst_methylation)[1]<-"Gene"
	write.table(pairwise_fst_methylation,paste0(outfolder,"/gene_pairwise_fst_smps_",name_file,".txt"),sep="\t",row.names=T,col.names=T,quote=F)
}

















# postponed for matrix scripts


# Make distance-like matrices with FST values in:

#library(usedist)

#for (rowis in 1:nrow(pairwise_fst_methylation)){
#	lineis <- pairwise_fst_methylation[rowis,2:ncol(pairwise_fst_methylation)]
#	if ( sum(lineis) == 0 ) {
#		next
#	}
#	if ( sum(lineis) != 0 ) {
#		matis=matrix(0,6,6)
#		matis[lower.tri(matis,diag=F)]<-lineis
#		matis=matrix(matis,6,6)
#		matis[upper.tri(matis)]<-t(matis)[upper.tri(matis)] 
#		colnames(matis) <-  sapply(strsplit(names(data)[-1],"Pop"),"[[",2) ; rownames(matis) <- sapply(strsplit(names(data)[-1],"Pop"),"[[",2)
		
		# Match Env files order: B K N S G E
#		matis1<-as.matrix(dist_subset(matis,c("B","K","N","S","G","O")))
		#print(matis1)
#		write.csv(file=paste0("gene_region/gene_fst_distmatrix/","gene_fst_dist_smp_",pairwise_fst_methylation[rowis,1],".csv"),matis1)
#	}
#}

#------------------------------

# Number of gene regions with sum > 0
#gene_region/gene_fst_distmatrix/ |wc -l
# 4194






