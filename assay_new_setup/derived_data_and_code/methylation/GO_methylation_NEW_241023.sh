



## ------------------ Get methylation coverage files
	# Lin have already run the mapping and methylation calling.
	# Documentation of how here: /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/run.sh
	# Coverage files can be found here: ls /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/ | grep "[SBKN][0-9][0-9].cov$"

### NEW deduplicated files run below:

#srun --mem=16g --pty /bin/bash

conda activate assay_new

#subset files for plotting

#for file in `ls /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/ | grep "[SBKN][0-9][0-9].cov$"`
#do
#awk 'BEGIN{srand()}{if(rand()<=.001)print $0}' /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/${file} > assay_new_setup/derived_data_and_code/methylation/smps/${file/.cov/_sub.cov}
#done

# ----------------- plot coverage in R and calculate top percent threshold cutoffs
# Plot (R).

'files<-dir(,"cov$")
png("assay_new_setup/derived_data_and_code/methylation/smps/Cdepth.tiff",height=700,width=1200)
par(mfrow=c(4,5),mar=c(2,2,3,1))
total_depths=vector()
for(file in files)
{
	a<-read.table(file)
	a<-a[,5]+a[,6]
	total_depths<-append(total_depths,a)
	i<-a<100
	a<-a[i]
	plot(table(a),main=sub(".sub.cov$","",file))
	grid()
}
dev.off()

percentscore=length(total_depths)/100*1
numberis=head(tail(sort(total_depths),percentscore));numberis
#[1] 32 32 32 32 32 32

percentscore=length(total_depths)/100*0.5
numberis=head(tail(sort(total_depths),percentscore));numberis
#[1] 39 39 39 39 39 39

percentscore=length(total_depths)/100*0.2
numberis=head(tail(sort(total_depths),percentscore));numberis
#[1] 61 61 61 61 61 61

'
#-------------------------------------- Back in bash, Apply coverage thresholds:
	

#===========    # UNZIP
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/new_files_241023/ | grep "[SBKN][0-9][0-9].deduplicated.bismark.cov.gz$"
    #/faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/ | grep "[SBKN][0-9][0-9].cov$"`
do
	echo $file #sbatch assay_new_setup/derived_data_and_code/methylation//coverage_filter_awk.sh /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/$file assay_new_setup/derived_data_and_code/methylation//smps/${file/.cov/_d5D32.cov} 5
    cp assay_new_setup/derived_data_and_code/methylation/smps/new_files_241023/$file assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file
    sbatch --account=spider2 --wrap="gunzip assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file"
done

#=========== Add strand and context information to the coverage file.
	# relatively time and memory expensive
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/*deduplicated.bismark.cov`
do
    sbatch --account=spider2 --mem=100G --wrap="Rscript assay_new_setup/derived_data_and_code/methylation/addSC.r $file ${file/bismark./}"
done

#===========# COVERAGE FILTER - cut at low 10 and high 32 (corresponding to top 1%)
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/ | grep "[SBKN][0-9][0-9].deduplicated.cov$"`
    #/faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/ | grep "[SBKN][0-9][0-9].cov$"`
do
	#sbatch assay_new_setup/derived_data_and_code/methylation//coverage_filter_awk.sh /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/$file assay_new_setup/derived_data_and_code/methylation//smps/${file/.cov/_d5D32.cov} 5
    sbatch assay_new_setup/derived_data_and_code/methylation/coverage_filter_awk.sh assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/${file/.cov/_d10D32.cov} 10
        # script infile outfile lowcut
done

#===========    # REZIP
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/ | grep "[SBKN][0-9][0-9].deduplicated.bismark.cov$"`
do
	#sbatch --acechoicount=spider2 --wrap="gzip -9 assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file"
	rm assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file
	# or just remove them!
done


## ---------------- Intersect to regions:

    # for picking out TSS and upstreams region, see bottom of script.
    
## ---------------- Intersect to regions:

# Genes and tss regions were extracted as in : /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/GO_overview.sh

#look here: less /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/


##### Gene region:

#cd /faststorage/project/spider2/assay_study_anne/methylation/

mkdir assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/ | grep "d10D32.cov$"`
do
	sbatch --account=spider2 --mem=50G --wrap="bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/$file -b /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/gene_pos_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/${file/.cov/_gene.cov}"
done



###################################################
## ------------------------- Extract cpgs CHH and CHG
###################################################

#conda activate assay

for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/*region`
do 
    mkdir $directory/CPG
	mkdir $directory/CHH
	mkdir $directory/CHG
	for file in `ls $directory | grep "_d10D32_" | grep ".cov$"`
	do
		sbatch --account=spider2 --wrap="grep "CpG" $directory/${file} > $directory/CPG/${file/.cov/.cpg}"
		sbatch --account=spider2 --wrap="grep "CHH" $directory/${file} > $directory/CHH/${file/.cov/.chh}"
		sbatch --account=spider2 --wrap="grep "CHG" $directory/${file} > $directory/CHG/${file/.cov/.chg}"
	done
done


## -------------------------- Calculate weighted methylation level for all the picked out regions
			# also check number of C's in each region retained after filtering
			
			
# To use for calculating methylation level difference per gene among populations.
# Using Jesper's measurement (methylation coverage ratio).
# Gene, coverage of methy, coverage of unmethy.

for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/C[A-Z][A-Z]`
do 
	for file in `ls $directory | grep ".d10D32_" | grep ".cpg$\|.chh$\|.chg$"`
	do
		echo $directory/$file
		extension=${file##*.}
		echo $extension
		echo ${file/.$extension/.$extension.sum}		
		sbatch assay_new_setup/derived_data_and_code/methylation/calculates_sums.sh $directory $file	
	done
done
	# Outfile loooks like this: Gene	methylated_reads		unmethylated_reads		count_of_Cs_in_region
	# m = methylated , u = unmethylated

############################################################
## Employ the scripts that calculate methylation level

#sacct -S 2020-03-11T15:20 --format="JobID,JobName,CPUTime,MaxRSS,ReqMem,NodeList,NNodes,State,TotalCPU,AveCPU,NCPUS"

for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/C[A-Z][A-Z]`
do 
	for file in `ls $directory | grep ".d10D32_" | grep ".cpg.sum$\|.chh.sum$\|.chg.sum$" |head -n1`
	do
		echo $directory/$file
		extension=${file#*.}
		echo $extension
		sbatch assay_new_setup/derived_data_and_code/methylation/cal_weighted_mlvl.sh $directory $file	
	done
done

mv cal_weighted_mlvl.sh*out assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/outfiles/

## -------------------------------- Make files ready for input into DSS for modelling

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

#cd /faststorage/project/spider2/assay_study_anne/methylation/

chmod a+x assay_new_setup/derived_data_and_code/methylation/make_input_for_dss.sh

for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/C[A-Z][A-Z]`
do 
	for file in `ls $directory | grep ".d10D32_" | grep ".cpg$\|.chh$\|.chg$"`
	do
		echo $directory/$file
		sbatch assay_new_setup/derived_data_and_code/methylation/make_input_for_dss.sh $directory $file	
	done
done
	# Outfile (_inputDSS.cpg) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads
	# Outfile (_inputDSS.cpg.counts) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads	count_of_Cs_in_region
# Check here what happens. if we get the right folders (CHH CPG CHG)

mv make_input_for_dss.sh*out assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/outfiles/

############################################
# Run binominal model in R
############################################

#srun --mem=10g --pty /bin/bash
#conda activate assay; R



# ------------- R code to run DSS for methylation data per region


chmod a+x assay_new_setup/derived_data_and_code/methylation/dss_run.sh

#for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/*region`
#do 
#	for file in `ls $directory | grep ".d10D32_" | grep "_inputDSS.cpgf$"|head -n1`
#	do
#		echo /faststorage/project/spider2/assay_study_anne/methylation/$directory/$file $file
#		#sbatch dss_run.sh /faststorage/project/spider2/assay_study_anne/methylation/$directory/ $file	
#	done
#done


for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/gene_region/C[A-Z][A-Z]`
do 
	for file in `ls $directory | grep ".d10D32_" | grep "_inputDSS.cpgf$\|_inputDSS.chhf$\|_inputDSS.chgf$" |head -n1`
	do
		echo $directory  $file
		sbatch assay_new_setup/derived_data_and_code/methylation/dss_run.sh $directory $file	
	done
done

mv dss_run.sh*out assay_new_setup/derived_data_and_code/methylation/smps/new_files_output/outfiles/



# conclusion:
	# CHG has some genes (4) with population effect on methylation level, but none with temperature effect.
	# CHH has no genes with effect of either population or temperature.
	# CPG has 1537 genes with population effect on methylation level, but none with temperature effect.
		
		# NEW: deduplicated:
		#	type	Temperature	population
		#	CHG		1	68
		#	CHH		0	1
		#	CPG		2	1537









############-------------------- alternative analysis in DSS:-----------########################

	# This code was run as a test to check whether results were consistent. The results were not used in publication, since the same pattern was observed.
	# run the analysis on all SMPs, then intersect to gene afterwards.
	# counts from files: /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/*poptemp*.cpg
	
## Cpg files used here was extracted from the .cov files exported by bismarkmethylationextractor, and filtered for a coverage above 5. So still pretty low depth/cov filtering. It is not too important, since DSS will also take that into account.
#cd assay_new_setup/derived_data_and_code/methylation/smps/second_analysis_DSS/
#directoryis="assay_new_setup/derived_data_and_code/methylation/smps/second_analysis_DSS/"

#for file in `ls /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/*.cpg`
#do
#	filenewname=`basename $file`
#	awk '{ print $1,$2,$5+$6,$5; }' $file >  $directoryis${filenewname/.cpg/_dssINPUT2.cpg}
#done

###--------- Filter for 10 cov.

#for file in `ls $directoryis/*.cpg`
#do
#	filenewname=`basename $file`
#	awk ' $3 >= 10' $file >  $directoryis${filenewname/.cpg/_10cov.txt} ##NOT good
#done


#--------- run DSS ---------------
#directory=assay_new_setup/derived_data_and_code/methylation/smps/second_analysis_DSS/
#chmod a+x $directory/../../dss_run_ver2.sh
#for file in `ls $directory | grep "_dssINPUT2_10cov.txt$"|head -n1`
#do
#	echo $directory/$file
#	sbatch $directory/../../dss_run_ver2.sh $directory/ $file	
#done
	
#	# this will produce files with results for each SMP, afterwards, intersect to get only SMPs in gene regions.

##--------------- Run bedtools intersect to get DSS results for gene regions: ----------------------------#

 ##Obs, I may need to modify output from DSS to have two position columns?
#file=`grep "DSS_result_all_" $directory`
#awk '{print $1,$2,$2,$3,$4,$5,$6,$7,$8}' $file > ${file/.txt/_mod.txt}

##### Gene region:

#cd /faststorage/project/spider2/assay_study_anne/methylation/
#fileis=${file/.txt/_mod.txt}

#qx --no-scratch -m 50G "bedtools intersect -loj -a $fileis -b /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/gene_pos_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/gene_region/${fileis/.txt/_gene.txt}"
		# could be a "-wo" or a "wa" instead?


# Make R script to devide regions into categories:
	# Effect of population
	# Effect of Temperature
	# Effect of both





#--------------- Pick out additional upstream TSS regions

	# pick one 1000bp upsteam and another 1000-2000 bp upstream
	
anno="/faststorage/project/spider2/dumicola_genome/annotation/v2.gff" #annotation

# Non overlapping tss regions.
	# First pick out whether it is a "tss" and what strand it is on (for/rev) --> define region as tss +/- upstream downstream
	#actually, I don't care if they overlap. I guess I only care if they are identical.
upstream=1000
downstream=0
awk -v upstream=$upstream -v downstream=$downstream 'BEGIN{FS="\t";OFS="\t"}\
$3=="tss" && $7=="+"{print $1,$4,$5,$9,$4-upstream,$5+downstream,$7}\
$3=="tss" && $7=="-"{print $1,$4,$5,$9,$5-upstream,$4+downstream,$7}' $anno\
 | sort -s -n -k 3,3 | sort -s -n -k 2,2 | sort -s -k 1,1 | awk '\
 BEGIN{FS="\t";OFS="\t";getline;SEQ=$1;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7}\
 RST<0{RST=1}\
 REN<0{REN=1}\
 $1!=SEQ{print SEQ,ST,EN,RST,REN,STRA,GENE;SEQ=$1;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7;next}\
 $2>EN{print SEQ,ST,EN,RST,REN,STRA,GENE;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7;next}\
 $3>EN{EN=$3}\
 END{print SEQ,ST,EN,RST,REN,STRA,GENE}' > assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_pos_region.txt
	# If SEQ is different on this and next line, print and next
	# if startposition is larger than end position on former line, print and next (ie, if next start is NOT witin the range of current sequence)
	# if current end position is larger than former end, make the new end position the variable.

#Check that a region is not reported twice
cut assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_pos_region.txt -f 1,4,5|uniq|wc -l
# 34442
cut assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_pos_region.txt -f 1,4,5|wc -l
# 34442

#cut out relevant fields to new file (seq, st_reg, en_reg,strand,gene name)
cut assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_pos_region.txt -d$'\t' -f 1,4,5,6,7 >assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_region.txt
rm assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region/tss1000_pos_region.txt

	
	# pick one 1000bp upsteam and another 1000-2000 bp upstream
	
	
cd /faststorage/project/spider2/assay_study_anne/methylation/
anno="/faststorage/project/spider2/dumicola_genome/annotation/v2.gff" #annotation

# Non overlapping tss regions.
	# First pick out whether it is a "tss" and what strand it is on (for/rev) --> define region as tss +/- upstream downstream
	#actually, I don't care if they overlap. I guess I only care if they are identical.
upstream=2000
downstream=-1000
awk -v upstream=$upstream -v downstream=$downstream 'BEGIN{FS="\t";OFS="\t"}\
$3=="tss" && $7=="+"{print $1,$4,$5,$9,$4-upstream,$5+downstream,$7}\
$3=="tss" && $7=="-"{print $1,$4,$5,$9,$5-upstream,$4+downstream,$7}' $anno\
 | sort -s -n -k 3,3 | sort -s -n -k 2,2 | sort -s -k 1,1 | awk '\
 BEGIN{FS="\t";OFS="\t";getline;SEQ=$1;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7}\
 RST<0{RST=1}\
 REN<0{REN=1}\
 $1!=SEQ{print SEQ,ST,EN,RST,REN,STRA,GENE;SEQ=$1;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7;next}\
 $2>EN{print SEQ,ST,EN,RST,REN,STRA,GENE;ST=$2;EN=$3;RST=$5;REN=$6;GENE=$4;STRA=$7;next}\
 $3>EN{EN=$3}\
 END{print SEQ,ST,EN,RST,REN,STRA,GENE}' > assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_pos_region.txt
	# If SEQ is different on this and next line, print and next
	# if startposition is larger than end position on former line, print and next (ie, if next start is NOT witin the range of current sequence)
	# if current end position is larger than former end, make the new end position the variable.

#Check that a region is not reported twice
cut assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_pos_region.txt -f 1,4,5|uniq|wc -l
# 34440
cut assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_pos_region.txt -f 1,4,5|wc -l
# 34442
################## OBS!!!! FIND overlap (similar regions)

#cut out relevant fields to new file (seq, st_reg, en_reg,strand,gene name)
cut assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_pos_region.txt -d$'\t' -f 1,4,5,6,7 >assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_region.txt
rm assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/tss2000_pos_region.txt


#Intesecting regions

##### tss region:

cd /faststorage/project/spider2/assay_study_anne/methylation/
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d5D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b //faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/tss_region/tss_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/tss_region/${file/.cov/_tss.cov}"
done
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d10D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b //faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/tss_region/tss_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/tss_region/${file/.cov/_tss.cov}"
done


##### tss 1000 region:
cd /faststorage/project/spider2/assay_study_anne/methylation/
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d5D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region//tss1000_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region//${file/.cov/_tss1000.cov}"
done
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d10D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region//tss1000_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/tss1000_region//${file/.cov/_tss1000.cov}"
done

##### tss 2000 region:
cd /faststorage/project/spider2/assay_study_anne/methylation/
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d5D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b assay_new_setup/derived_data_and_code/methylation/smps//tss2000_region/tss2000_region.txt> assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/${file/.cov/_tss2000.cov}"
done
for file in `ls assay_new_setup/derived_data_and_code/methylation/smps/ | grep "d10D32.cov$"`
do
	qx --no-scratch -m 50G "bedtools intersect -wo -a assay_new_setup/derived_data_and_code/methylation/smps/$file -b assay_new_setup/derived_data_and_code/methylation/smps//tss2000_region/tss2000_region.txt > assay_new_setup/derived_data_and_code/methylation/smps/tss2000_region/${file/.cov/_tss2000.cov}"
done

