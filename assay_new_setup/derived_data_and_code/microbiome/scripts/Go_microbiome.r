## Auther: Anne Aagaard
## File for analyzing microbiome data:
# import data.
# filter according to filter used in Sources_of_variation (prevalence threshold of 25%) 
# analyze using anova to find temperature, population and interaction effects

#conda activate assay_copy; R			

# setwd("/faststorage/project/spider2/assay_study_anne/microbiome/")
setwd("/faststorage/project/spider2/assay_study_anne/")
library(agricolae)
library(scales)
library(vegan)
#Plot is further below


micro_dat = read.csv(file="assay_new_setup/derived_data_and_code/microbiome/data_input/Assay-tax-asvdata-fraction.csv",sep=",", row.names=1)
meta_dat = read.csv(file="assay_new_setup/derived_data_and_code/microbiome/data_input/assayf-metadata.csv",sep=",")

microdat_newcolnames <- micro_dat[,1]
taxa_dat <- micro_dat[,1:8]
micro_dat <- t(micro_dat[,9:ncol(micro_dat)])
colnames(micro_dat) <- microdat_newcolnames
	#merge with meta data.
micro_dat_all <- merge(meta_dat,micro_dat,by.x=1,by.y=0)
micro_dat <- micro_dat_all[,8:ncol(micro_dat_all)]	;	rownames(micro_dat) <- micro_dat_all$X
pop_nest_vec <- as.factor(paste0(micro_dat_all[,2],"_",micro_dat_all[,3]))

	# count presence in Population/nest groups
total_count_PopNest <- aggregate(x=micro_dat[,1:ncol(micro_dat)],by=list(pop_nest_vec),FUN=length)
presence_count_PopNest <- aggregate(x=micro_dat[,1:ncol(micro_dat)],by=list(pop_nest_vec),FUN=function(c)sum(c!=0))
	
	# count presence in 25% of nests within a population (aggregate on previous aggregate, using new pop_vector)
pop_vector <- sapply(strsplit(as.character(total_count_PopNest[,1]),"_"),"[[",1)
total_count_PopNest_Pop <- aggregate(x=presence_count_PopNest[,2:ncol(presence_count_PopNest)],by=list(pop_vector),FUN=length)
presence_count_PopNest_Pop <- aggregate(x=presence_count_PopNest[,2:ncol(presence_count_PopNest)],by=list(pop_vector),FUN=function(c)sum(c!=0))
	
	# calculate prevalence of nests in a population containing said column
prevalence <- presence_count_PopNest_Pop[2:ncol(presence_count_PopNest_Pop)]/total_count_PopNest_Pop[2:ncol(total_count_PopNest_Pop)]
ASVs_above <- which( colSums(prevalence > .25) >= 1)
NROW(ASVs_above)
	# 78

	# Make new datafranme with only selected ASVs
df_new <- merge(cbind(ASVs_above, NA),t(micro_dat), by=0,all.x=T,all.y=F)
rownames(df_new) <- df_new[,1]			;		df_new <- df_new[,4:ncol(df_new)]
new_df_all <- merge(meta_dat,t(df_new),by.x=1,by.y=0, all.x=F, all.y=T)
new_df <- new_df_all[,8:ncol(new_df_all)]

# Then do anova? with temp + pop + interaction  (should colony be tested in a model?)


# Function for running ANOVA and extracting results. Adapted from metaboanalystR
aov_between <- function(x){
  unlist(summary(aov(x ~ aov_facA*aov_facB)), use.names=F)[c(13,14,15,17,18,19)]
}
# run function in apply
aov_facA <- new_df_all$Population	;	aov_facA_name <- "Population" # pop vector 
aov_facB <- new_df_all$Temp 		;	aov_facB_name <- "Temperature"# temp vector
aov_mat <- t(apply(as.matrix(new_df), 2, aov_between))
	# make pvalue correction
p_correction <- "fdr"

aov_mat2 <- cbind (aov_mat, p.adjust(aov_mat[,4], p_correction),
							p.adjust(aov_mat[,5], p_correction),
							p.adjust(aov_mat[,6], p_correction))
						 
colnames(aov_mat2) <- c(paste(aov_facA_name, "(F.val)", sep = ""), 
                        paste(aov_facB_name, "(F.val)", sep = ""),
                        paste("Interaction", "(F.val)", sep = ""),
                        paste(aov_facA_name, "(raw.p)", sep = ""), 
                        paste(aov_facB_name, "(raw.p)", sep = ""),
                        paste("Interaction", "(raw.p)", sep = ""), 
                        paste(aov_facA_name, "(adj.p)", sep = ""), 
                        paste(aov_facB_name, "(adj.p)", sep = ""), 
                        paste("Interaction", "(adj.p)", sep = ""))

aov_mat2 <- as.data.frame(aov_mat2[, c(1,4,7,2,5,8,3,6,9),drop=F] )


write.csv(file="assay_new_setup/derived_data_and_code/microbiome/data_output/aov_betweenSbjs_all_ASVs.csv",aov_mat2)

	## Now filter on significant etc.
pop_sig <- rownames(aov_mat2)[which(aov_mat2$"Population(adj.p)" < .05)]
temp_sig <- rownames(aov_mat2)[which(aov_mat2$"Temperature(adj.p)" < .05)]
pop__temp_sig <- rownames(aov_mat2)[which(aov_mat2$"Population(adj.p)" < .05 & aov_mat2$"Temperature(adj.p)" < .05 )]
int_sig <- rownames(aov_mat2)[which(aov_mat2$"Interaction(adj.p)" < .05)]

count_sign_df <- data.frame(NROW(pop_sig),NROW(temp_sig),NROW(pop__temp_sig),NROW(int_sig))
colnames(count_sign_df) <- c("Pop_effect","Temp_effect","Pop&Temp_effect","Interaction")

write.csv(file="assay_new_setup/derived_data_and_code/microbiome/data_output/Count_sign_ASVs.csv",count_sign_df)


### Which ASV's is it?

microbe_df <- taxa_dat[0,]
for(ASV_is in int_sig){
	print(taxa_dat[grep(paste0(ASV_is,"$"), taxa_dat[,1])		,		] )#add line ending symbol
}

microbe_df <- taxa_dat[0,]
for(ASV_is in pop_sig){
	microbe_df[nrow(microbe_df)+1,] <- taxa_dat[grep(paste0(ASV_is,"$"), taxa_dat[,1])		,		] #add line ending symbol
}
microbe_df$"study_effects" <- c(rep("pop",4),"pop + int","pop")

write.csv(file="assay_new_setup/derived_data_and_code/microbiome/data_output/ASV_with_effects_taxa.csv",microbe_df)


## make pca of taxa with effects

new_df_pca <- new_df_all[,c(1,2,8:ncol(new_df_all))]
pop_column <- new_df_all[,2]	;	temp_column <- new_df_all[,4]	
rownames(new_df_pca) <- new_df_all[,1]	; new_df_pca <- t(new_df_pca[,-c(1,2)])
new_df_pca <- new_df_pca[ which(rownames(new_df_pca) %in% pop_sig), ]



#pca_popExp_genes <- prcomp(t(na.omit(vst_counts_pca)))
pca_micro <- prcomp(t(new_df_pca))
summa <- summary(pca_micro)
summa
'Importance of components:
                          PC1    PC2    PC3    PC4     PC5     PC6
Standard deviation     0.4709 0.3332 0.2155 0.1473 0.11299 0.01481
Proportion of Variance 0.5357 0.2683 0.1122 0.0524 0.03085 0.00053
Cumulative Proportion  0.5357 0.8040 0.9162 0.9686 0.99947 1.00000
'

popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")
pop_order <- vector()
popcols <- popcols[ grep("Ndumo", names(popcols), invert=T) ]
for( each in pop_column ){	pop_order <- append(pop_order,grep(each,names(popcols)))	}
popcols2 <- popcols[pop_order]



output_type = "pdf" #png #pdf
#output_type = "png" #png #pdf

if (	output_type == "pdf"	)	{ 	
	pdf("assay_new_setup/derived_data_and_code/microbiome/data_output/pca_ASV_pop_plots.pdf",width=7,height=7) 		
}else if ( 	output_type == "png"	){
	png("assay_new_setup/derived_data_and_code/microbiome/data_output/pca_ASV_pop_plots.png",width=7,height=7,units="in",res=400) 	
}
par(mar=c(4.7,4.2,3,1.5), mfrow=c(2,2),cex=.7 , oma=c(0,0,1,0))
plottit <- ordiplot(pca_micro,choices=c(1,2),type="n", 
			xlim=c(-1,1),ylim=c(-1,1), cex=0.8, main="Microbiome, population",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC2 (",round(summa$importance[2,2],2)*100,"%)"))
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))
legend("topright",legend=sort(names(popcols)),col=popcols[sort(names(popcols))],pch=1, 
			inset=c(-(.39),-.21), xpd=NA, ncol=2, box.col="grey90")

plottit <- ordiplot(pca_micro,choices=c(2,3),type="n", 
			xlim=c(-1,1),ylim=c(-1,1), cex=0.8, main="Microbiome, population",
			xlab=paste0("PC1 (",round(summa$importance[2,1],2)*100,"%)"), ylab=paste0("PC3 (",round(summa$importance[2,3],2)*100,"%)"))
points(plottit, "sites",col=popcols2, pch=1)#
ordiellipse(plottit,popcols2,conf=0.80, col=levels(as.factor(popcols2)))

screeplot(pca_micro, type="l", main="Screeplot Microbiome")
biplot(pca_micro, main="Biplot Microbiome")
dev.off()


#Plot the asvs individually
popcols3 <- factor(popcols2,levels=unique(popcols2))
pop_column1 <- factor(pop_column,levels=unique(pop_column))
jit_col <- pop_column1	;	levels(jit_col) <- seq(-.5,.5,length.out=4)	; jit_col <-  as.numeric(as.character(jit_col))


if (	output_type == "pdf"	)	{ 	
	pdf("assay_new_setup/derived_data_and_code/microbiome/data_output/trends_popASVs.pdf",width=7,height=4) 		
}else if ( 	output_type == "png"	){
	png("assay_new_setup/derived_data_and_code/microbiome/data_output/trends_popASVs.png",width=7,height=4,units="in",res=400) 	
}
par(mfrow=c(2,3),oma=c(2,1,2,1), cex=.7,cex.lab=.7, cex.main=.7,cex.axis=.7, mar=c(2,2.5,3,0))
for ( rowss in 1:nrow(new_df_pca) ){
	grepped_name <- microbe_df[	grep(paste0(rownames(new_df_pca)[rowss],"$"), microbe_df$Row.names),	which(colnames(microbe_df)=="Species")	]
	print(grepped_name)
	plot(NULL,NULL,"n",xlim=c(17,30),ylim=c(0,1),xaxt="n",
		main=paste0(rownames(new_df_pca)[rowss],"\n",grepped_name),
		ylab=NA, xlab=NA)
	axis(1,at=seq(19,29,by=2),labels=NA, col.ticks="grey70")
	axis(1,at=c(19,29),labels=c(19,29))
	points(temp_column+jit_col,new_df_pca[rowss,],col=popcols2, cex=.7)
	for (pop in levels(pop_column1)){
		abline(lm(new_df_pca[rowss, which(pop_column1==pop) ] ~ temp_column[ which(pop_column1==pop) ]),col=levels(popcols3)[ which(levels(pop_column1)==pop) ]	)
	}
}
legend(-2.5,3.5,legend=levels(pop_column1), pch=1,col=levels(popcols3), xpd=NA, horiz=T,bty="n",cex=.85)
mtext("Relative abundance", side=2, line=0, outer=T,xpd=NA,cex=.7)
mtext(expression(paste("Acclimation Temperature (", degree,"C)")),side=1,line=.5,outer=T,xpd=NA,cex=.7)
dev.off()
