#!/home/anneaa/miniconda3/envs/assay_copy/bin/Rscript

#-------------------- Parse arguments
argsis<-commandArgs(trailingOnly=TRUE)
directory_is=argsis[1]		#smps/gene_region/		
		# "assay_new_setup/derived_data_and_code/methylation/smps/gene_region/"
file_is=argsis[2]		#B15_d5D32_gene_inputDSS.cpgf
		#S15_d10D32_gene_inputDSS.chg.cpgf
assay_directory=argsis[3]

dir.create(paste0(directory_is,"/DSS_output/"))
outdirectory=paste0(directory_is,"/DSS_output/")
				# directory_is="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/dss_files_SOV/"
				# file_is="PopA.d10D30.cov.gene_inputDSS.cpgf"
				# assay_directory="assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG/"

assay_temp_to_get=25

#-------------------- Setup
library(DSS)
require(bsseq)



#-------------------- Parse file name
	# Get the name that is shared for all sample files:
#patternis=paste(sapply(strsplit(file_is,split="_"),"[[",2),sapply(strsplit(file_is,split="_"),"[[",3),	paste0(	sapply(strsplit(file_is,split="_|\\."),"[[",4) , "." , sapply(strsplit(file_is,split="_|\\."),"[[",5)	),sep="_")
patternis=paste(sapply(strsplit(file_is,split="_|\\."),"[[",4),sapply(strsplit(file_is,split="_"),"[[",2), sep="_")
	# patternis="gene_inputDSS.cpgf"
	# Adjust to your needs
files_are <- list.files(path=directory_is,pattern=glob2rx(paste0("*",patternis)),full.names=T) # "$" makes sure the string ends
	# subset to have correct pops: 

files_are_assay <- list.files(path=assay_directory, pattern=glob2rx(paste0("*",assay_temp_to_get,"*",patternis)),full.names=T)
files_are <- append(files_are,files_are_assay)
	# split for output naming:
outname=	paste0(	sapply(strsplit(patternis,"_input"),"[[",1),	"_" ,	toupper( sapply(strsplit(patternis,"\\.|f"),"[[",2))	)

	# d5D32_gene
	
	
#-------------------- Import files into list of dataframes
samplelist=vector()
varlist=list()
dict = list(A = 'G', B = 'O', C = 'B', D = 'K', E = 'S', N = 'N')
	# Modify SOV sample names:
	# sync files order :	"A",		"B",		"C",		"D",			"O",			"N" 			(from lins sync files - see /faststorage/project/spider2/dumicola_genome/Lin_20180103_mapping_v2/run_aln.sh)
		# 					A=gobabis, 	B=Etosha,	C=Betta, 	D=Karasburg, 	E=Stampriet, 	N=Ndumo			(From sandbjerg slides)
for (each in files_are){
	message(each)
	sub_name <- sapply(strsplit(each,"//"), "[[", 2) # primary split on file name
	sample_name=sapply(strsplit(sub_name,"\\.|_"),"[[",1) #Gets the sample name (uniq) --> Adjust to your needs
	if( grepl("Pop",sample_name) ){
		sample_name <- gsub("Pop","",sample_name)
		sample_name1 <- unlist( sapply(1:6, function(i) replace(sample_name, sample_name == names(dict[i]), dict[i])) )# changes into corresponsing dictonary name
		if ( NROW( unique( sample_name1 ) ) == 1 ){ sample_name2 <- sample_name1[1] }else{ sample_name2 <- sample_name1 [ grep(sample_name,sample_name1, invert=T) ] }
		if ( grepl( "G|N", sample_name2 ) ){ next }
		sample_name <- paste0(sample_name2,"_t15")
	}else if( grepl("15|19|23|25|29",sample_name) ){
		sample_name <- gsub("15|19|23|25|29","_t17",sample_name)
		if ( grepl("N",sample_name) ){ sample_name <- gsub("N", "O", sample_name) }
	}
	nam <- paste("data_", sample_name, sep = "")	# assigns a name to a variable for use in next step
	assign(nam,read.table(each,header=T)) # this assigns the read-in table to a variable with the name that is inside the nam variable - equivalent to data_POPtreat
	varlist[[nam]]<- get(nam)	# get the data now assigned to nam and add to list varlist
	samplelist <- append(samplelist,sample_name) # add also the sample name to a vector
	rm(nam)
}

#-------------------- Make design matrix

poplist <- sapply(strsplit(samplelist,"_"),"[[",1) # split samplenames to pop
templist <- sapply(strsplit(samplelist,"_"),"[[",2) # split samplenames to pop

design<-as.data.frame(cbind(poplist,templist)); colnames(design)<- c("population","samplingtime") # make design matrix


#-------------------- Make BSobject needed for DSS
BSobj <- makeBSseqData( varlist, samplelist)
BSobj
'
An object of type "BSseq" with
  36006 methylation loci
  8 samples
has not been smoothed
All assays are in-memory
'


#-------------------- Fit the model in DSS

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~population+samplingtime)


#-------------------- Test the effect of each term in the model and write out data

terms_are = c("samplingtime","population")
	#terms_are = c("Temperature","population","Interaction")
		#### First term:
term=terms_are[1]  # adjust the term to fit your data
test_1 = DMLtest.multiFactor(DMLfit, term=term) 
head(test_1)
		##       chr     pos     stat        pvals         fdrs
		## 1273 chr1 2930315 5.280301 1.289720e-07 0.0006448599
		## 4706 chr1 3321251 5.037839 4.708164e-07 0.0011770409

	# write all results to file
write.table(test_1, file=paste0(outdirectory,"/DSS_result_addit_Assay",assay_temp_to_get,"_",outname,"_",term,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	#test_1<-read.table(file=paste0(directory_is,"/DSS_result_addit_",outname,"_",term,".txt"),header=T,sep="\t")
	# pick out results that are significant and below FDR 0.1
test_1_sub <- test_1[ which(test_1$pvals < 0.05 & test_1$fdrs < 0.1),]
write.table(test_1_sub, file=paste0(outdirectory,"/DSS_result_addit_Assay",assay_temp_to_get,"_",outname,"_",term,"_significant.txt"),row.names=F,col.names=T,quote=F,sep="\t")


		#### Second term:
term=terms_are[2]  # adjust the term to fit your data
test_2 = DMLtest.multiFactor(DMLfit, term=term)
head(test_2)
	# write all results to file
write.table(test_2, file=paste0(outdirectory,"/DSS_result_addit_Assay",assay_temp_to_get,"_",outname,"_",term,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
	# pick out results that are significant and below FDR 0.1
test_2_sub <- test_2[ which(test_2$pvals < 0.05 & test_2$fdrs < 0.1),]
write.table(test_2_sub, file=paste0(outdirectory,"/DSS_result_addit_Assay",assay_temp_to_get,"_",outname,"_",term,"_significant.txt"),row.names=F,col.names=T,quote=F,sep="\t")

		## Combine terms table
require("plyr")

### export counts table

effect_of_tab <- data.frame(NROW(test_1_sub),NROW(test_2_sub))
	#effect_of_tab <- c(NROW(test_1_sub),NROW(test_2_sub),NROW(test_3_sub))
colnames(effect_of_tab) <- terms_are
write.table(effect_of_tab, file=paste0(outdirectory,"/DSS_result_counteffects_Assay",assay_temp_to_get,"_",outname,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

# make a merged results table:

merged_table <- join_all(list(test_1,test_2), by ="chr" , type = 'left')
merged_table <- merged_table[,-which(colnames(merged_table) == "pos")]

colnames(merged_table)<-c("region_name",
		paste0("stat_",substr(terms_are[1],1,1)),paste0("pvalue_",substr(terms_are[1],1,1)), paste0("fdr_",substr(terms_are[1],1,1)),
		paste0("stat_",substr(terms_are[2],1,1)),paste0("pvalue_",substr(terms_are[2],1,1)), paste0("fdr_",substr(terms_are[2],1,1)))


write.table(merged_table, file=paste0(outdirectory,"/DSS_result_all_Assay",assay_temp_to_get,"_",outname,".txt"),row.names=F,col.names=T,quote=F,sep="\t")


#Maybe call DMRs:
	#**Step 4*. DMRs for multifactor design can be called using {callDMR} function:

#callDMR(test_pop, p.threshold=0.05)
	##      chr   start     end length nCG   areaStat
	## 33  chr1 2793724 2793907    184   5  12.619968
	## 413 chr1 3309867 3310133    267   7 -12.093850


