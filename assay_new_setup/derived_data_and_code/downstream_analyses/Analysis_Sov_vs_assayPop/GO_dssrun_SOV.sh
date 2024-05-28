
#conda activate assay

## -------------------------------- Make files ready for input into DSS for modelling

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

#cd /faststorage/project/spider2/assay_study_anne/methylation/
chmod a+x make_input_for_dss.sh

directory_datafiles=../Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region

outdirectory=assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/dss_files_SOV
working_script_directory=assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop

mkdir -p $outdirectory
for file in `ls $directory_datafiles | grep ".d10D30" | grep ".gene.cpg$"`
do
    echo $directory_datafiles/$file
    #echo $file
    sbatch $working_script_directory/make_input_for_dss_SOV.sh $directory_datafiles $file $outdirectory	
done
	# Outfile (_inputDSS.cpg) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads
	# Outfile (_inputDSS.cpg.counts) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads	count_of_Cs_in_region
	

# Run binominal model in R (DSS)

#srun --mem=10g --pty /bin/bash
#conda activate assay; R


mv make_input_for_dss.sh*out outfiles/make_input_outs/

# ------------- R code to run DSS for methylation data per region


#chmod a+x dss_run.sh

#for directory in `ls -d assay_new_setup/derived_data_and_code/methylation/smps/*region`
#do 
	#for file in `ls $directory | grep ".d5D32_" | grep "_inputDSS.cpgf$"|head -n1`
	#do
	#	echo /faststorage/project/spider2/assay_study_anne/methylation/$directory/$file $file
	#	sbatch dss_run.sh /faststorage/project/spider2/assay_study_anne/methylation/$directory/ $file
	#done


# I cannot get dss to run on four populations with no replicates
    # It will not fit a linear model due to lacking degrees of freedom (which makes sence) --> thus I cannot make similar analyses across assay and SOV
   
    # Or - I should take all correlations to climatic parameters (before I only took temperature related parameters)
    # or - do SOV correlations to PCs for both assay (each temp) and SOV data subset.
        # Check overlap
        # However, PCs are proxies for population differences. It is not a great analysis.
        # Also turned out to not reveal any that were shared between any Assay correlations and SOV correlations using B,K,O,S only and 99.99% threshold.
        # (see GO_mantel_tests_ASSAY.sh)
    # BETTER (and USED IDEA):
        # Run DSS using SOV as timepoint 1 (or year/samplingtime 15), and use Assay25 as timepoint 2 (or year/samplingtime 17)
        # This should run with 8 samples (4 pops per 2 samplingtimes)
            # Here we assume that populations are comparable, as nests sampled are not the same.
                # This seems like a valid assumption, assuming that some amount of meth follows genetic signal, and gen div is very low.
        # This will give us an idea of whether methylation pattern can change over time (two generations) in the wild.
        # If that is the case, methylation patterns are plastic.   
            # Just ASSAY temperatures and time frame is not enough to push it.
        
#conda activate assay
directory_assay_files=assay_new_setup/derived_data_and_code/methylation/smps/gene_region/CPG

for file in `ls $outdirectory | grep ".d10D30." | grep "_inputDSS.cpgf$"|head -n1`
do
    echo $outdirectory $file
    sbatch $working_script_directory/dss_run_SOV.sh $outdirectory $file	$directory_assay_files
done


mv dss_run.sh*out outfiles/dss_run_outs/

# Conclusion
    # ~17000 genes change according to sampling time (additive model of pop+samplingtime, below FDR=.1 and pval=0.05)
    # ~500 genes different between populations
        # Same picture whether I use Assay 15, 19, 23, 25 or 29 together with SOV data.
    # What could explain this apparent plastic signal? (data contrary to ASSAY temperature acclimations, that does not provide plastic signal)
        # Genes methylation pattern being plastic in response to all sorts of things, including strong developmental plasticity?
            # Just not across the time frame and temperature pressures experiences in common garden in adult subadult-adultspiders.
        # Methylations vary a lot between individuals or nests?
            # still, it seems that 17000 genes is a lot to be explained just by this
            # could be tested on the level of nests with David Fishers data?
                # data from two time points in time series, time point one and time point just before they die off.
                # As far as we remember, 5 spiders pooled?

                # datafiles should look like one of the following formats: preferably be bismark bedgraph coverage files
                    # see Lin's script /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/run.sh
                '
                /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/B29.bismark.cov
                SEQ_00001       13      13      0       0       1
                SEQ_00001       28      28      0       0       2
                SEQ_00001       34      34      0       0       2
                SEQ_00001       45      45      0       0       3
                SEQ_00001       49      49      0       0       4
                
                /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/B29.cov
                SEQ_00001       13      13      0       0       1       +       CHG
                SEQ_00001       28      28      0       0       2       +       CHH
                SEQ_00001       34      34      0       0       2       +       CHH
                SEQ_00001       45      45      0       0       3       +       CHH
                SEQ_00001       49      49      0       0       4       +       CHG   
                
                /faststorage/project/spider2/Methylation_data/Lin_20180608_bismark_v2_assay/SMP/B29.cpg
                SEQ_00001       1159    1159    71.4285714285714        5       2       -       CpG
                SEQ_00001       1176    1176    62.5    5       3       -       CpG
                SEQ_00001       1283    1283    37.5    3       5       -       CpG
                SEQ_00001       2243    2243    85.7142857142857        6       1       +       CpG
                SEQ_00001       2324    2324    76.9230769230769        10      3       +       CpG 
                '
# Modify below analysis to analyse David Fishers methylation data from Time series dumicola data

## ----------------------- nest replicability of mehylations on gene level? -------------------##

    # 5 nests (T52-T56) from the Otavi population were sampled at two time points ( an earlier time point (maybe t2) and the last timepoint before the nest dies off )
    # Pools of 3 spiders per sample
    # How deeply was it sequenced? coverage of x?

    ## Get coverage files and filter to coverage ten and maybe a high threshold?
directory_samples="/home/anneaa/spider2/faststorage/David/BACKUP/clean/"
    # Cpg files used here was extracted from the .cov files exported by bismarkmethylationextractor.
        # depth is not too important, since DSS will also take that into account.
        # still I will filter on cov 10
        # files are called:
            # $directory_samples/*/bedgraph_output.gz.bismark.cov.gz

##--------- Filter for 10 cov. per cpg site
    
    # Coverage file looks like this:  <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>

outdirectory="assay_new_setup/derived_data_and_code/downstream_analyses/Analysis_Sov_vs_assayPop/estimate_nest_methylation_variability/cov_files/"
for directo in  `ls -d $directory_samples/*`
do 
    sample_id=`echo $directo | rev | cut -f1 -d'/'|rev`
    file_for_covfilter=$directo/"bedgraph_output.gz.bismark.cov.gz"
        # no dependencies:
     copy_id=$(sbatch --parsable --wrap="cp $file_for_covfilter $outdirectory/$sample_id'_bedgraph_output.gz.bismark.cov.gz'")
     unzip_id=$(sbatch --parsable --dependency=afterany:$copy_id --wrap="gunzip $outdirectory/$sample_id'_bedgraph_output.gz.bismark.cov.gz'")
        # no dep:
    #unzip_id=$(sbatch --parsable --wrap="gunzip $outdirectory/$sample_id'_bedgraph_output.gz.bismark.cov.gz'")
    
    file_for_covfilter1=$outdirectory/$sample_id"_bedgraph_output.gz.bismark.cov"
        # dependencies:
     cov_process_id=$(sbatch --parsable --dependency=afterany:$copy_id:$unzip_id --wrap="awk '\$5+\$6 >= 10' $file_for_covfilter1 > ${file_for_covfilter1/.cov/_d10.cov}")
     sbatch --dependency=afterany:$copy_id:$unzip_id:$cov_process_id --wrap="gzip -9 $file_for_covfilter1"
        # fewer dependencies:
    #cov_process_id=$(sbatch --parsable --dependency=afterany:$unzip_id --wrap="awk '\$5+\$6 >= 10' $file_for_covfilter1 > ${file_for_covfilter1/.cov/_d10.cov}")
    #sbatch --dependency=afterany:$unzip_id:$cov_process_id --wrap="gzip -9 $file_for_covfilter1"
done
    # loop:
        # copies coverage file from David Fishers directory
        # Unzips it
        # filter on coverage ( escape doller signs inside awk command when submitting as one liner )
        # compress file again in new directory
    # submits to cluster with dependencies for earlier process to finish


    #---------- Intersect to gene region


##### Gene region:
    # make sure to use correct genone annotation (so posititons are correct for genes)
    # make new position file to intersect with in case of new genome. 
        # He did, according to the scaffold names (hicScaffolds) in the coverage files
        # See /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/GO_overview.sh
        # /faststorage/project/spider2/Methylation_data/AAL_20200207_fishing_regions_SOV/gene_region/GO_get_gene_region_sunrem.sh
conda activate assay_copy
anno="/home/anneaa/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/annotations/DUM_Ste_All_gene.gff3"
    # old: "/faststorage/project/spider2/dumicola_genome/annotation/v2.gff" #annotation

# Non overlapping gene regions.
        # First pick out whether it is a "gene" and what strand it is on (for/rev) --> define region as gene +/- upstream downstream
        #actually, I don't care if they overlap. I guess I only care if they are identical.
awk 'BEGIN{FS="\t";OFS="\t"}\
$3=="gene" && $7=="+"{print $1,$4,$5,$9,$7}\
$3=="gene" && $7=="-"{print $1,$4,$5,$9,$7}' $anno\
 | sort -s -n -k 3,3 | sort -s -n -k 2,2 | sort -s -k 1,1 | awk '\
 BEGIN{FS="\t";OFS="\t";getline;SEQ=$1;ST=$2;EN=$3;GENE=$4;STRA=$5}\
 $1!=SEQ{print SEQ,ST,EN,STRA,GENE;SEQ=$1;ST=$2;EN=$3;GENE=$4;STRA=$5;next}\
 $2>EN{print SEQ,ST,EN,STRA,GENE;ST=$2;EN=$3;GENE=$4;STRA=$5;next}\
 $3>EN{EN=$3}\
 END{print SEQ,ST,EN,STRA,GENE}' > $outdirectory/gene_pos_region_hic_genome.txt
        # If SEQ is different on this and next line, print and next
        # if startposition is larger than end position on former line, print and next (ie, if next start is NOT witin the range of current sequence)
        # if current end position is larger than former end, make the new end position the variable.

#Check that a region is not reported twice
cut $outdirectory/gene_pos_region_hic_genome.txt -f 1,2,3|uniq|wc -l
# 33567
cut $outdirectory/gene_pos_region_hic_genome.txt -f 1,2,3|wc -l
# 33567

##### intersect to Gene region:

for file in `ls $outdirectory/*d10.cov`
do
	sbatch --mem=16G --wrap="bedtools intersect -wo -a $file -b $outdirectory/gene_pos_region_hic_genome.txt > ${file/_bedgraph_output.gz.bismark_d10.cov/_d10_gene.cov}"
done
    # writes overlapping stuff from a, followed by b, followed by the number of bp overlapping (not that I need last part)

## -------------------------- Calculate weighted methylation level for all the picked out regions
			# also check number of C's in each region retained after filtering
					
# To use for calculating methylation level difference per gene among populations.
# Using methylation coverage ratio.
    # Gene, coverage of methy, coverage of unmethy.

for file in `ls $outdirectory/*d10_gene.cov`
do
    echo $file
    sbatch $working_script_directory/calculates_sums_TIMESER.sh $file	
done
	# Outfile loooks like this:     Gene	methylated_reads		unmethylated_reads		count_of_Cs_in_region
	# m = methylated , u = unmethylated

## Employ the scripts that calculate methylation level

for file in `ls $outdirectory/*d10_gene.cov.sum | head -n1`
do
    echo $file
    sbatch $working_script_directory/cal_weighted_mlvl_TIMESER.sh $outdirectory/ $file
done

## ---------------------- make DSS input to analyse differences among nests
    # nests T52-T56
    # time points _1 and _2

    # In DSS test effect of nest

## -------------------------------- Make files ready for input into DSS for modelling

# Files should look like this: CHR	Pos	TOTAL_READS	METHYLATED_reads

#cd /faststorage/project/spider2/assay_study_anne/methylation/

for file in `ls $outdirectory/*d10_gene.cov`
do
    echo $file
    #echo $file
    sbatch $working_script_directory/make_input_for_dss_TIMESER.sh $file
done
	# Outfile (_inputDSS.cpg) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads
	# Outfile (_inputDSS.cpg.counts) loooks like this: Gene	Pos_first_cpg	total_reads		methylated_reads	count_of_Cs_in_region
	

# ------------- R code to run DSS for methylation data per region

#conda activate assay_copy

for file in `ls $outdirectory | grep ".d10D30." | grep "_inputDSS.cpgf$"|head -n1`
do
    echo $outdirectory $file
    sbatch $working_script_directory/dss_run_TIMESER.sh $outdirectory $file	$directory_assay_files
done
