Sys.setenv(R_LIBS_USER="~/myRlib_3_5_0")
.libPaths("~/myRlib_3_5_0")
	
	# Modified from Lins script
library(Biostrings)
#f.ref<-"../Lin_20180104_bismark_v2/ref_v2/ref_v2.fasta"
f.ref<-"../Jesper/bismark/genome/Stegodyphus_dumicola_genome_v2.fasta" 
ref<-readDNAStringSet(f.ref)

input<-commandArgs(trailingOnly=T)
f.in<-input[1]
f.out<-input[2]

dict.strand<-c("+","-")
names(dict.strand)<-c("C","G")
dict.context<-c(
	"CpG","CpG","CpG","CpG",
	"CHG","CHG","CHG",
	"CHH","CHH","CHH",
	"CHH","CHH","CHH",
	"CHH","CHH","CHH")
names(dict.context)<-c(
	"CGA","CGT","CGG","CGC",
	"CAG","CTG","CCG",
	"CAA","CAT","CAC",
	"CTA","CTT","CTC",
	"CCA","CCT","CCC")

c.in<-file(f.in,"r")
repeat
{
	a<-try(read.table(c.in,stringsAsFactors=F,nrows=1e6),silent=T)
	if(class(a)=="try-error")break
	x<-as.character(subseq(ref[a[,1]],a[,2],a[,3]))
	strand<-dict.strand[x]
	i<-x=="C"
	x[i]<-as.character(subseq(ref[a[,1]][i],a[i,2],a[i,3]+2))
	x[!i]<-as.character(reverseComplement(subseq(ref[a[,1]][!i],a[!i,2]-2,a[!i,3])))
	context<-dict.context[x]
	a<-cbind(a,strand,context)
	write.table(a,f.out,append=T,sep="\t",row.names=F,col.names=F,quote=F)
}
close(c.in)
