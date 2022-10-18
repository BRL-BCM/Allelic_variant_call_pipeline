setwd("/Volumes/WD_BLACK/GG/ThesisProj/exRNA/refGenome/masked_hg38_WcmmV")
tab.0=read.table("ERCC.cov.hg38.hg38_gnomad30_genome.common005.txt")
popFreq=0.05
tab.0=tab.0[,1:6]

temp=tab.0
ind=NULL
temp=temp[(order(temp$V1,temp$V2)),]
posDistance=temp$V2-c(0,temp$V2[-dim(temp)[1]])
while (any(posDistance<=50 & posDistance>=0)) {
  ind=which(posDistance<=50 & posDistance>=0)
  ind.rm=c()
  
  ind.rm=sapply(ind, function(i) {ifelse(temp$V6[i]>temp$V6[i-1],i-1,i)})

  # for (i in ind){
  #   i.rm=ifelse(temp$V6[i]>temp$V6[i-1],i-1,i)
  #   ind.rm=c(ind.rm,i.rm)
  # }
  if(length(ind.rm)>0) temp=temp[-unique(ind.rm),]
  posDistance=temp$V2-c(0,temp$V2[-dim(temp)[1]])
  ind=NULL
}

tab=temp
dim(tab)
save(tab,file = sprintf("hg38_gnomad.common%s.ERCC.cov.hg38.50bp.RData",popFreq))
rm(list=ls())




## reformat and adding header ##
popFreq=0.05
load(sprintf("hg38_gnomad.common%s.ERCC.cov.hg38.50bp.RData",popFreq))


tab$V3=NULL
tab$V1=sapply(tab$V1, function(x) sub("chr","",x))
tab$V6=paste("NS=1:DP=50:","AF=",tab$V6,sep="")
colnames(tab)=c("#CHROM","POS","REF","ALT","INFO")
tab$FORMAT="GT:GQ:DP:HQ"
tab$GNOMAD="0|1:49:50:60,60"
tab$FILTER="PASS"
tab$QUAL="67"
tab$ID="."
tab=tab[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GNOMAD")]
tab=unique(tab)
dim(tab)

# adding base for deletion/insersion #
require(seqinr)
require(Biostrings)

ref=readDNAStringSet("../masked_hg38/broad-hg38.masked.fa")
del.ind=which(tab$ALT=="-")
# tab[del.ind[1:10],"REF"]
# del.base=
#   sapply(del.ind[1:10], #round(length(del.ind)/10)
#          function(x)
#            as.character(subseq(ref[[paste("chr", tab$`#CHROM`[x], sep="")]],
#            start=(tab$POS[x]-1),
#            width = 10)))
tab=tab[-del.ind,]
temp=tab
# temp$POS[del.ind]=temp$POS[del.ind]-1
# temp$REF[del.ind]=paste(del.base,temp$REF[del.ind],sep = "")
# temp$ALT[del.ind]=del.base

ins.ind=which(tab$REF=="-")
tab=tab[-ins.ind,]

save(tab,file = sprintf("hg38_gnomad.common%s.ERCC.cov.hg38.50bp.reformed.RData",popFreq))


fileDate=paste0(unlist(strsplit(as.character(Sys.Date()),"-")),collapse = "")
fileSource=sprintf("hg38_gnomad.common%s.ERCC.cov.hg38.50bp.reformed",popFreq)

vcfHeader=sprintf("##fileformat=VCFv4.0
##fileDate=%s
##source=%s
##reference=broad-hg38.masked.fa
##phasing=complete
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
",fileDate,fileSource)
cat(vcfHeader,file="vcfHeader.txt")

load(sprintf("hg38_gnomad.common%s.ERCC.cov.hg38.50bp.reformed.RData",popFreq))
write.table(tab,file="vcfBody.txt",quote = F,sep = "\t",row.names = F)
system(sprintf("cat vcfHeader.txt vcfBody.txt > gnomad%s.50bp.reformed.vcf",popFreq))
system("rm vcfHeader.txt vcfBody.txt")
rm(list=ls())
