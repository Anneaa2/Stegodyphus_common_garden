#!/usr/bin/env Rscript
# Script for plotting PCA's for the metabolites

# plots a panel with 4 plots:
  # 1 PC1 and PC2
  # 2 PC2 and PC3
  # 3 Screeplot
  # 4 Biplot

# input should be
  # 1 PCA variable
  # 2 color factor (f.ex. pop colors)
'
args = commandArgs(trailingOnly=TRUE)

if (length(args) !=2 ) {
  stop("Two arguments must be supplied (input PCA element and color factor/vector).n", call.=FALSE)
} else if (length(args) == 2) {
  # default output file
  PCA_element = args[1]
  color_vector = args[2]
}
'
#PCA_element=pca_pop_features; factor_to_use=factor_is; color_vector=colvec1
pca_plot_metabo <- function( PCA_element, factor_to_use, color_vector){
	#class(PCA_element) <- "prcomp"
	par(mfrow=c(2,2), mar=c(4.7,4.2,4,2.5),oma=c(0,0,0,.5), cex=.5, cex.lab=.9, cex.axis=.9)
	summa <- summary(PCA_element)
	uni_cols <- unique(color_vector)
	names(uni_cols) <- unique(names(color_vector))
	#print(names(uni_cols) )
	#pc1 pc2
	ordii<-ordiplot(PCA_element, type="n",choices=c(1,2), main=datname)
	points(ordii, "sites",col=color_vector,pch=1)	; plot_xrange <- max(ordii$sites[,1])-min(ordii$sites[,1])
	ordiellipse(PCA_element,factor_is,conf=0.85,w=NULL, col=uni_cols[sort(names(uni_cols))],choices=c(1,2))#,xlim=xlimsare,ylim=xlimsare)
	legend("topright",legend=sort(names(uni_cols)),col=uni_cols[sort(names(uni_cols))],pch=1, inset=c(-(.25),0),xpd=NA)
	#pc2 pc3
	ordii<-ordiplot(PCA_element, type="n",choices=c(2,3), main=datname)
	points(ordii, "sites",col=color_vector,pch=1)	; plot_xrange <- max(ordii$sites[,1])-min(ordii$sites[,1])
	ordiellipse(PCA_element,factor_is,conf=0.85,w=NULL, col=uni_cols[sort(names(uni_cols))],choices=c(2,3))#,xlim=xlimsare,ylim=xlimsare)
	#legend("topright",legend=sort(names(uni_cols)),col=uni_cols[sort(names(uni_cols))],pch=1, inset=c(-(.15),0),xpd=NA,bty="n")

	#scree
	screeplot(PCA_element, type="l", main=paste0("Screeplot ",datname))
	#biplot
	biplot(PCA_element$x[,1:2],PCA_element$rotation[,1:2], main=paste0("Biplot ",datname))
}


## modify PLDA plot function:

# note from Xia: internal variables and functions not to be modified by users
# note from Xia: This is only for web version
.on.public.web <- FALSE; # only TRUE when on metaboanalyst web server

# note from Xia: note, this is usually used at the end of a function
# note from Xia: for local, return itself; for web, push to global environment
.set.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    mSet <<- mSetObj;
    return (1);
  }
  return(mSetObj);
}

.get.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    return(mSet)
  }else{
    return(mSetObj);
  }
}
 # actual function to modify 

### Fix the levels issue - figure out some other way. STILL SOMETHING OFF WITH COLORS/PCH ELLIPSES???

PlotPLS2DScore_mod <- function(mSetObj=NA, imgName, factor_to_use, col_vector=NULL, format="png", dpi=72, width=NA, inx1, inx2, reg=0.95, show=1, grey.scale=0, use.sparse=FALSE){

	mSetObj <- .get.mSet(mSetObj);

	imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
	if(is.na(width)){
		w <- 9;
	}else if(width == 0){
		w <- 7.2;
	}else{
		w <- width;
	}
	h <- w;

	if( factor_to_use == "pop" ){ # set factor to show
		cls1 <- as.factor(substr(mSetObj$dataSet$cls,1,1))
		cls.type <- mSetObj$dataSet$cls.type
	}else if( factor_to_use == "temp" ){
		cls1 <- mSetObj$dataSet$facB
		cls.type <- mSetObj$dataSet$cls.type
	}

	if( is.vector(col_vector) ){
		if(length(col_vector) == length(levels(cls1))){ #chech if lenghts are equal and define colVec
			#colVec <- col_vector[order(names(col_vector))] #gets embedded into getcolorscheama function
			colVec <- col_vector #[order(names(col_vector))] #gets embedded into getcolorscheama function
				#colvec1 = factor_is = mSet$dataSet$cls
					#col_vector1 = popcols1[as.character(colvec1)]
			###print(colVec)
		}else{
			message("error - color vetor length not equal to the number of levels of data factor.\nTo determine colors yourself, make them equal length.")
		}	
	}

	mSetObj$imgSet$pls.score2d <- imgName;

	lv1 <- mSetObj$analSet$plsr$scores[,inx1];
	lv2 <- mSetObj$analSet$plsr$scores[,inx2];
	xlabel <- paste("Component", inx1, "(", round(100*mSetObj$analSet$plsr$Xvar[inx1]/mSetObj$analSet$plsr$Xtotvar,1), "%)");
	ylabel <- paste("Component", inx2, "(", round(100*mSetObj$analSet$plsr$Xvar[inx2]/mSetObj$analSet$plsr$Xtotvar,1), "%)");

	Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
	par(mar=c(5,5,3,3));
	text.lbls <- substr(rownames(mSetObj$dataSet$norm),1,12) # some names may be too long

	# obtain ellipse points to the scatter plot for each category

	if(cls.type=="integer"){
		cls <- as.factor(as.numeric(levels(cls1))[cls1]);
	}else{
		cls <- cls1;
	}

	lvs <- levels(cls);
	pts.array <- array(0, dim=c(100,2,length(lvs)));
	for(i in 1:length(lvs)){
		inx <- cls1 == lvs[i];
		groupVar <- var(cbind(lv1[inx],lv2[inx]), na.rm=T);
		groupMean <- cbind(mean(lv1[inx], na.rm=T),mean(lv2[inx], na.rm=T));
		pts.array[,,i] <- ellipse::ellipse(groupVar, centre = groupMean, level = reg, npoints=100);
	}

	xrg <- range(lv1, pts.array[,1,]);
	yrg <- range(lv2, pts.array[,2,]);
	x.ext<-(xrg[2]-xrg[1])/12;
	y.ext<-(yrg[2]-yrg[1])/12;
	xlims<-c(xrg[1]-x.ext, xrg[2]+x.ext);
	ylims<-c(yrg[1]-y.ext, yrg[2]+y.ext);

	## cols = as.numeric(dataSet$cls)+1;
	cols <- GetColorSchema_mod(cls, colVec=colVec, grey.scale==1);
	uniq.cols <- unique(cols);
	#print(uniq.cols)

	plot(lv1, lv2, xlab=xlabel, xlim=xlims, ylim=ylims, ylab=ylabel, type='n', main="Scores Plot");
	#plot(lv1, lv2, xlab=xlabel, xlim=xlims, ylim=ylims, ylab=ylabel, main="Scores Plot");
	grid(col = "lightgray", lty = "dotted", lwd = 1);

	# make sure name and number of the same order DO NOT USE levels, which may be different
	###legend.nm <- unique(as.character(sort(cls)));
	legend.nm <- unique(as.character((cls)));
	## uniq.cols <- unique(cols);
	#print(legend.nm)

	pchs <- GetShapeSchema_mod(mSetObj, cls, show=show, grey.scale=grey.scale) #uses facB section

	#print(pchs)
	uniq.pchs <- unique(pchs);
	#print(uniq.pchs)
	names(uniq.pchs) <- legend.nm;
	uniq.pchs <- uniq.pchs[order(names(uniq.pchs))]

	## BHAN: when same color is choosen for black/white; it makes an error
	# names(uniq.cols) <- legend.nm;
	if (length(uniq.cols) > 1) {
		names(uniq.cols) <- legend.nm;
		###print(uniq.cols)
		uniq.cols <- uniq.cols[order(names(uniq.cols))]
		###print(uniq.cols)
	}

	legend.nm <- legend.nm[order(legend.nm)]
	###print(legend.nm)

	# draw ellipse
	for(i in 1:length(lvs)){
		if ( length(uniq.cols) > 1) {
			polygon(pts.array[,,i], col=adjustcolor(uniq.cols[lvs[i]], alpha=0.2), border=NA);
		} else {
			polygon(pts.array[,,i], col=adjustcolor(uniq.cols, alpha=0.2), border=NA);
		}
		if(grey.scale) {
			lines(pts.array[,,i], col=adjustcolor("black", alpha=0.5), lty=2);
		}
	}

	if(grey.scale) {
		cols <- rep("black", length(cols));
	}
	if(show==1){ # display sample name set on
	text(lv1, lv2, label=text.lbls, pos=4, xpd=T, cex=0.75);
	points(lv1, lv2, pch=pchs, col=cols);
	}else{
		if (length(uniq.cols) == 1) {
			points(lv1, lv2, pch=pchs, col=cols, cex=1.0);
		} else {
			if(grey.scale == 1 | (exists("shapeVec") && all(shapeVec>=0))){
				my.cols <- adjustcolor(cols, alpha.f = 0.4);
				my.cols[pchs == 21] <- "black";
				points(lv1, lv2, pch=pchs, col=my.cols, bg=adjustcolor(cols, alpha.f = 0.4), cex=1.8);
			}else{
				points(lv1, lv2, pch=21, bg=adjustcolor(cols, alpha.f = 0.4), cex=2);
			}
		}
	}


	if(grey.scale) {
		uniq.cols <- "black";
	}

	if(length(uniq.cols) != length(levels(cls))){

		if(cls.type=="integer"){
			names(cols) <- as.numeric(cls)
		}else{
			names(cols) <- as.character(cls)
		}

		match.inx <- match(levels(cls), names(cols))
		uniq.cols <- cols[match.inx]
	}

	if(length(uniq.pchs) != length(levels(cls))){
		#if(mSetObj$dataSet$type.cls.lbl=="integer"){
		#  names(pchs) <- as.numeric(cls)
		#}else{
		names(pchs) <- as.character(cls)
		#}

		match.inx <- match(levels(cls), names(pchs))
		uniq.pchs <- pchs[match.inx]
	}

	legend("topright", legend = legend.nm, pch=uniq.pchs, col=uniq.cols);

	dev.off();
	return(.set.mSet(mSetObj));
}





# same as in metaboanalystr - but it would not load?

GetShapeSchema_mod <- function(mSetObj=NA, my.cls, shapeVec=NULL, show.name, grey.scale){
	mSetObj <- .get.mSet(mSetObj);
	if(exists("shapeVec") && all(shapeVec >= 0) && length(shapeVec) == length(levels(my.cls))){
		#sps <- rep(0, length=length(my.cls));
		#clsVec <- as.character(my.cls)
		#grpnms <- names(shapeVec);
		#for(i in 1:length(grpnms)){
		#  sps[clsVec == grpnms[i]] <- shapeVec[i];
		#}
		sps = shapeVec[as.character(my.cls)] # creates named vector of colors - use as input.
		shapes <- sps;
	}else{
		#message("entered else")
		if(show.name | grey.scale){
			shapes <- as.numeric(my.cls)+1;
		}else{
			shapes <- rep(21, length(my.cls));
		}
	}
	return(shapes);
}
pal_18 <- c("#e6194B", "#3cb44b", "#4363d8", "#42d4f4", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45",
                  "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075");

GetColorSchema_mod <- function(my.cls, colVec=NULL, grayscale=F){

	lvs <- levels(my.cls);
	grp.num <- length(lvs);
	if(grayscale){
		dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
	}else if(exists("colVec") && !any(colVec =="#NA") && length(colVec) == length(levels(my.cls))){
		###print(colVec)
		dist.cols <- colVec;
		###print(dist.cols)
	}else{
		if(grp.num <= 18){ # update color and respect default
			dist.cols <- pal_18[1:grp.num];
		}else{
			dist.cols <- colorRampPalette(pal_18)(grp.num);
		}
	}

	#colors_are <- vector(mode="character", length=length(my.cls));
	#for(i in 1:length(lvs)){
	#	colors_are[my.cls == lvs[i]] <- dist.cols[i];
	#}
	col_vector_is = dist.cols[as.character(my.cls)] # creates named vector of colors - use as input.
	#print(colors_are)
	#print(col_vector_is)
	#print(my.cls)
	return (col_vector_is);
}
