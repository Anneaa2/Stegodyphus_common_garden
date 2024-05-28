
 ## For SOV results and temperature
 ## /faststorage/project/spider2/AAL_20190722_SOV_correlations/correlation_results/only_temp_prec/genewise_SMPS_res/
	## temp_prec_partial_mantel_res_genes_smps_refpop_ALL.csv

# cd /faststorage/project/spider2/assay_study_anne/phenotypes/data/ ; conda activate assay_copy; R
# setwd("/faststorage/project/spider2/assay_study_anne/phenotypes/data/") 
# Rscript
library(scales)
library(readxl)
library(reshape2)
outplace="assay_new_setup/derived_data_and_code/phenotypes/results/"


#output_type = "png"
output_type = "pdf"

if (	output_type == "pdf"	)	{ 	
	pdf(paste0(outplace,"time_to_temp_plots.pdf"),width=7,height=5) 		
}else if ( 	output_type == "png"	){
	png(paste0(outplace,"time_to_temp_plots.png"),width=7,height=5,units="in",res=400) 	
}
#par(mfrow=c(1,2)).
layout(matrix(c(1,1,1,2,2),nrow=1))

# CTmaxx
dat_ctmax <- as.data.frame(read_excel("assay_new_setup/derived_data_and_code/phenotypes/data/time_to_temperature.xlsx",sheet=1))
dat_ctmax <- dat_ctmax[1:25,c(1:5,7:8)]			; dat_ctmax[,6] <- as.numeric(dat_ctmax[,6])
xmax_is <- max(dat_ctmax[,1],na.rm=T)
ymax_is <- max(apply(dat_ctmax[,-1],1,max,na.rm=T),na.rm=T)
avg_is <- rowMeans(dat_ctmax[,-1],na.rm=T)
sd_is <- apply(dat_ctmax[,-1], 1, sd,na.rm=T)
N_per_timepoint <- apply(dat_ctmax[,-1], 1,function(i) NROW(which(!is.na(i))))
names(N_per_timepoint) <- dat_ctmax[,1]; N_per_timepoint
write.csv(N_per_timepoint, paste0(outplace,"time_to_temp_ctmax_N.csv")) 	

plot(dat_ctmax[,1],type="n", xlab="Time (min)",ylab=expression("Temperature ("*~degree*C*")"), xlim=c(0, xmax_is+1),ylim=c(42,ymax_is+1),main="CTmax", yaxt='n')
apply(dat_ctmax[,-1],2,points,x=dat_ctmax[,1],col=alpha("red",.4))
axis(2, las=2)
points(dat_ctmax[,1],avg_is)
arrows(dat_ctmax[,1], avg_is-sd_is, dat_ctmax[,1], avg_is+sd_is, length=0.05, angle=90, code=3)
abline(lm(avg_is~dat_ctmax[,1]))
mtext("a", line=1, adj=-0.1)

lm(avg_is~dat_ctmax[,1])
'Call:
lm(formula = avg_is ~ dat_ctmax[, 1])

Coefficients:
   (Intercept)  dat_ctmax[, 1]
       42.9916          0.0979
'

# ccr:
dat_ccr <- as.data.frame(read_excel("assay_new_setup/derived_data_and_code/phenotypes/data/time_to_temperature.xlsx",sheet=2))
xmax_is <- max(dat_ccr[,1],na.rm=T)
ymax_is <- max(apply(dat_ccr[,-1],1,max,na.rm=T),na.rm=T)
avg_is <- rowMeans(dat_ccr[,-1],na.rm=T)
sd_is <- apply(dat_ccr[,-1], 1, sd,na.rm=T)
N_per_timepoint <- apply(dat_ccr[,-1], 1,function(i) NROW(which(!is.na(i))))
names(N_per_timepoint) <- dat_ccr[,1]; N_per_timepoint
write.csv(N_per_timepoint, paste0(outplace,"time_to_temp_ccr_N.csv")) 	

plot(dat_ccr[,1],type="n", xlab="Time (min)",ylab=expression("Temperature ("*~degree*C*")"), xlim=c(0, xmax_is+1),ylim=c(0,ymax_is+1), main="CCRTemp", yaxt='n')
apply(dat_ccr[,-1],2,points,x=dat_ccr[,1],col=alpha("blue",.1))
axis(2, las=2)
points(dat_ccr[,1],avg_is)
arrows(dat_ccr[,1], avg_is-sd_is, dat_ccr[,1], avg_is+sd_is, length=0.05, angle=90, code=3)
abline(lm(avg_is~dat_ccr[,1]))
mtext("b", line=1, adj=-0.2)

lm(avg_is~dat_ccr[,1])
'Call:
lm(formula = avg_is ~ dat_ccr[, 1])

Coefficients:
 (Intercept)  dat_ccr[, 1]
      0.2537        0.5001
'
dev.off()


