#!/home/anneaa/miniconda3/envs/assay/bin/Rscript


###################							###################
#######					This Scripts calculates			#######
#######		 	Weighted methylation level per region	#######
#################							###################

	# Parse arguments
argsis<-commandArgs(trailingOnly=TRUE)
directory_is=argsis[1]		#smps/gene_region/	
							#directory_is="assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/CPG/"
file_is=argsis[2]		#B15_d5D32_gene.cpg.sum
						#file_is="B15.deduplicated_d10D32_gene.cpg.sum"
extension_is=argsis[3]  #chh.sum
						#extension_is="deduplicated_d10D32_gene.cpg.sum"

# grep("g10062",rownames(count_meth))
## MAIN

patternis=paste(sapply(strsplit(file_is,split="_"),"[[",2),sapply(strsplit(file_is,split="_"),"[[",3),sep="_")
#patternis="d5D32_gene.cpg.sum"
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
name_popgroup=sapply(strsplit(files[1],split="/|//|_|\\."),"[[",16)
regionis=sapply(strsplit(files[1],split="//|/|_|.c"),"[[",13)


for(i in 2:length(files))
{
	name_popgroup<-append(name_popgroup,	sapply(strsplit(files[i],split="/|//|_|\\."),"[[",16)	)

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
	head(data)
}
names(data)<-c(regionis,name_popgroup)
#35305 genes

# 						A=gobabis, 	B=Etosha,	C=Betta, 	D=Karasburg, 	E=Stampriet, 	N=Ndumo			(From sandbjerg slides)
#									(OTAWI)
#names_vec=c("G","O","B","K","S","N")
names(data)<-c(regionis,name_popgroup)
names(c_inreg)<-c(regionis,name_popgroup)
names(region_methyl_reads)<-c(regionis,name_popgroup)
names(region_total_reads)<-c(regionis,name_popgroup)

c_inreg$Mean <- rowMeans(c_inreg[,-1])

# Save results files:
#name_file=paste0(	sapply(strsplit(patternis,paste0(".",extension_is)),"[[",1) ,	"_"	, sapply(strsplit(extension_is,".sum"),"[[",1) 	)
name_file=paste0( sapply(strsplit(extension_is,".sum"),"[[",1) 	)
						# d5D32_gene										# chh or cpg or chg
	# patternis="d5D32_gene.cpg.sum"
	# name_file= d5D32_gene_chh
write.table(data, file=paste0(directory_is,"/wmlvl_allpops_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(c_inreg, file=paste0(directory_is,"/nrCpg_allpops_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(region_methyl_reads, file=paste0(directory_is,"/methreads_allpops_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(region_total_reads, file=paste0(directory_is,"/totalreads_allpops_",name_file,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
print(name_file)
message("Total number of regions: ", nrow(data),"         Number of regions with no methylation: ", NROW(which(rowSums(data[,-1])==0)))

warnings()
