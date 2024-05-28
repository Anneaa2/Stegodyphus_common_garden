#/faststorage/project/spider2/assay_study_anne/phenotypes/data/

setwd("/faststorage/project/spider2/assay_study_anne/phenotypes/")

ToD_raw <- read.table("data/Complete_data_file_ToD_270218.txt", header = TRUE, sep = "", dec = ',')

ToD <- ToD_raw[!grepl("E001", ToD_raw$Nest),]   ## Why do we remove this nest??
#Nest E001 has very low values in all acclimation temperatures, and is removed here########

w_tod <-  data.frame(ToD$Population,as.factor(ToD$Treatment),ToD$Mass,ToD$ToD)
colnames(w_tod) <- c("Population","Treatment","Mass","ToD")
	# remove NA
w_tod <-  na.omit(w_tod)


pdf(file="results/temp_of_death_vs_weight.pdf")
plot(w_tod$Mass,w_tod$ToD,	xlab="Mass", ylab="CTmax (ToD)")
dev.off()

# pop_colors <- c("red", "darkgoldenrod1", "green3", "steelblue1", "blue", "darkorchid4") # pops <- c("Otavi", "Windhoek", "Stampriet", "Betta", "Karasburg", "Botswana")
pop_colors <- c("red", "steelblue1", "green3", "blue") # otavi, betta, stampriet, karasburg

#subset - ToD - population

pdf(file="results/temp_of_death_vs_weight_pops.pdf")
par(mfrow=c(2,2))
plot(main="Otavi",	w_tod[which(w_tod$Population == "Otavi"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Population == "Otavi"),which(colnames(w_tod) =="ToD")],	col=pop_colors[1],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="Betta",	w_tod[which(w_tod$Population == "Betta"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Population == "Betta"),which(colnames(w_tod) =="ToD")],	col=pop_colors[2],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="Stampriet",	w_tod[which(w_tod$Population == "Stampriet"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Population == "Stampriet"),which(colnames(w_tod) =="ToD")],	col=pop_colors[3],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="Karasburg",	w_tod[which(w_tod$Population == "Karasburg"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Population == "Karasburg"),which(colnames(w_tod) =="ToD")],	col=pop_colors[4],	xlab="Mass", ylab="CTmax (ToD)")
dev.off()


#subset - ToD - treatment

accli_colors <- c("red", "orange", "yellow", "green", "blue") # otavi, betta, stampriet, karasburg


pdf(file="results/temp_of_death_vs_weight_accli.pdf")
par(mfrow=c(3,2))
plot(main="15",	w_tod[which(w_tod$Treatment == "15"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Treatment == "15"),which(colnames(w_tod) =="ToD")],	col=accli_colors[1],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="19",	w_tod[which(w_tod$Treatment == "19"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Treatment == "19"),which(colnames(w_tod) =="ToD")],	col=accli_colors[2],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="23",	w_tod[which(w_tod$Treatment == "23"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Treatment == "23"),which(colnames(w_tod) =="ToD")],	col=accli_colors[3],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="25",	w_tod[which(w_tod$Treatment == "25"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Treatment == "25"),which(colnames(w_tod) =="ToD")],	col=accli_colors[4],	xlab="Mass", ylab="CTmax (ToD)")
plot(main="29",	w_tod[which(w_tod$Treatment == "29"),which(colnames(w_tod) =="Mass")]	,	w_tod[which(w_tod$Treatment == "29"),which(colnames(w_tod) =="ToD")],	col=accli_colors[5],	xlab="Mass", ylab="CTmax (ToD)")
dev.off()

