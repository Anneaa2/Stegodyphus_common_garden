#!/home/anneaa/miniconda3/envs/assay/bin/Rscript


###################							###################
#######					This Scripts calculates			#######
#######		 	Weighted methylation level per region	#######
#################							###################

	# Parse arguments
argsis<-commandArgs(trailingOnly=TRUE)
directory_is=argsis[1]
	# directory_is="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/estimate_nest_methylation_variability/cov_files/"	
file_is=argsis[2]
	# file_is="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/estimate_nest_methylation_variability/cov_files//T52_1_d10_gene.cov.sum"
extension_is=argsis[3]  
	# extension_is="cov.sum"

## MAIN

patternis=sub("(.*_){2}","\\1",basename(file_is))	 #"d10_gene.cov.sum"  - finds the pattern of any character a number of times, followed by _. this should be found twice. The \\1 is a backreference, recalling the pattern found. (d10_ in this case)
	#paste(sapply(strsplit(basename(file_is),split="_"),"[[",2),sapply(strsplit(basename(file_is),split="_"),"[[",3),sep="_")
print(patternis)

files<-dir(directory_is,pattern=patternis,full.names=T)
	#[1] "smps/tss_region//B15_d10D32_tss.cpg.sum"
	#[2] "smps/tss_region//B19_d10D32_tss.cpg.sum"
	#[2] "smps/tss_region//B19_d10D32_gene.cpg.sum"
#	assay_new_setup/derived_data_and_code/methylation/smps/gene_region/K19_d10D32_gene.chg.sum

data<-read.table(files[1],stringsAsFactors=F)
	# number of c's in region
c_inreg<-data.frame(data[,1],data[,4],stringsAsFactors=F)
	# number of reads methylated
region_methyl_reads<-data.frame(data[,1],data[,2],stringsAsFactors=F)
	# Total reads
region_total_reads<-data.frame(data[,1],data[,2]+data[,3],stringsAsFactors=F)
	# calculate weighted methylation level
data<-data.frame(data[,1],data[,2]/(data[,2]+data[,3]),stringsAsFactors=F)
	# get namings
name_popgroup=sub("(.*_[0-9]).*","\\1",basename(file_is))
	# sapply(strsplit(files[1],split="/|//|_"),"[[",12)
regionis="Gene"
	#sapply(strsplit(files[1],split="//|/|_|.c"),"[[",14)

for(i in 2:length(files))
{
	name_popgroup<-append(name_popgroup,	sub("(.*_[0-9]).*","\\1",basename(files[i]))	)
		#sapply(strsplit(files[i],split="/|//|_"),"[[",12)
	bdata<-read.table(files[i],stringsAsFactors=F)
	bc_inreg<-data.frame(bdata[,1],bdata[,4],stringsAsFactors=F)
	B_region_methyl_reads<-data.frame(bdata[,1],bdata[,2],stringsAsFactors=F)
	B_region_total_reads<-data.frame(bdata[,1],bdata[,2]+bdata[,3],stringsAsFactors=F)
	bdata<-data.frame(bdata[,1],bdata[,2]/(bdata[,2]+bdata[,3]),stringsAsFactors=F)
	
		# Merge the data
	data<-merge(data,bdata,by.x=1,by.y=1,sort=F, all=T)
	c_inreg<-merge(c_inreg,bc_inreg,by.x=1,by.y=1,sort=F, all=T)
	region_methyl_reads<-merge(region_methyl_reads,B_region_methyl_reads,by.x=1,by.y=1,sort=F, all=T)
	region_total_reads<-merge(region_total_reads,B_region_total_reads,by.x=1,by.y=1,sort=F, all=T)
}
head(data)
names(data)<-c(regionis,name_popgroup)
#30414 genes
names(c_inreg)<-c(regionis,name_popgroup)
names(region_methyl_reads)<-c(regionis,name_popgroup)
names(region_total_reads)<-c(regionis,name_popgroup)

c_inreg$Mean <- rowMeans(c_inreg[,-1])

# Save results files:

name_file=gsub(".cov.sum","",patternis)
	# "d10_gene"
	# paste0(	sapply(strsplit(patternis,paste0(".",extension_is)),"[[",1) ,	"_"	, sapply(strsplit(extension_is,".sum"),"[[",1) 	)
						# d5D32_gene										# chh or cpg or chg
	# patternis="d5D32_gene.cpg.sum"
	# name_file= d5D32_gene_chh
write.table(data, file=paste0(directory_is,"/wmlvl_allnests_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(c_inreg, file=paste0(directory_is,"/nrCpg_allnests_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(region_methyl_reads, file=paste0(directory_is,"/methreads_allnests_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(region_total_reads, file=paste0(directory_is,"/totalreads_allnests_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
print(name_file)
message("Total number of regions: ", nrow(data),"         Number of regions with no methylation: ", NROW(which(rowSums(data[,-1])==0)))
	# Total number of regions: 30414         Number of regions with no methylation: 49

warnings()