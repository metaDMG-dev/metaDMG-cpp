library(data.table)
#[fvr124@dandycomp03fl metaPlots]$ cut -f3 -d_ Cambodia_IDs.txt |md5sum
#050133146ac7206c4fab0376184b8a17  -
#[fvr124@dandycomp03fl metaPlots]$ cut -f1 metadata.tsv |cut -f3 -d_|sed 1d|md5sum
#050133146ac7206c4fab0376184b8a17  -
ids_txt <- "Cambodia_IDs.txt"
mdata <- "metadata.tsv"
plot_outdir <- "res0/"


plotter <- function(x){
  est<-x[1:21]
#  print(est[1])
  d<-matrix(x[-c(1:21)],5)
  plot(d[3,],main=mean(d[2,]),pch=1,type='b',col=1,lwd=2)
  legend("top",paste(names(est),est,sep="="))
  lines(1:60,d[4,],pch=2,type='b',lwd=1,col=2)
}



fun <- function(fname,doplot=T){
##    df <- fread("/projects/caeg/data/production/LV30/0527/LV3005272684/20230927_HNLW5DSX5_LV7009026484/results/metadmg/dfit/LV3005272684_LV7009026484_collapsed.dfit.gz")
    df <- fread(fname)
    a<-as.matrix(df)
    a<-a[order(a[,1]),]
    naln <- apply(a,1,function(x) mean(matrix(x[22:321],5)[2,]))
    if(doplot)
      apply(a[naln>20,],1,plotter)
    invisible(a)
}


fns<-read.table(ids_txt,as.is=T)
fns<-matrix(list.files(sapply(fns,function(x) paste("/projects/caeg/data/production/",sep='/',x,"/results/metadmg/dfit/")),full.names=T),2)[2,]

##plot for each library
for(i in fns){
      onam <- paste0(plot_outdir,basename(i),".pdf")
   pdf(onam)
fun(i)
dev.off()
}


alldata <- sapply(fns,fun,doplot=F)
##> table(unlist(lapply(alldata,function(x) x[1,1])))
#
# 1
#43

#> mean(matrix(alldata[[1]][1,-c(1:21)],5)[2,])
#[1] 21251.07

dep<-as.numeric(read.table(mdata,as.is=T,he=T)[,2])
nal <- unlist(lapply(alldata,function(x) mean(matrix(x[1,-c(1:21)],5)[2,])))
pdf(paste0(plot_outdir,"nreads_dep.pdf"))
plot(dep,nal,main="nr of classifying reads over depths",xlab="Depth",ylab="Nr of alignments")
dev.off()


##get mapstatss
fns<-read.table(ids_txt,as.is=T)
fns<-list.files(sapply(fns,function(x) paste("/projects/caeg/data/production/",sep='/',x,"/stats/metadmg/aggregate/")),full.names=T)

d<-lapply(fns,function(x) fread(x,nrows=2)[1,4:5])
d<-matrix(unlist(d),2)


fns<-read.table(ids_txt,as.is=T)
fns<-list.files(sapply(fns,function(x) paste("/projects/caeg/data/production/",sep='/',x,"/stats/align/samtools_stats/")),full.names=T)

a <- sapply(fns,function(x) paste("grep mapped",x,"|head -n1|cut -f3"))
a<-as.numeric(unlist(lapply(a,system,intern=T)))

res<-rbind(a,d)
colnames(res)<-dep
pdf(paste0(plot_outdir,"summary.pdf"),width=21)
barplot(res,col=1:3,legend=c("Mapped","classified alignments","classified reads"))
dev.off()



##read nr of seqs sequenced
fns<-read.table(ids_txt,as.is=T)

oo <-c()
for(i in fns[,1]){
  bname <- list.files(sapply(i,function(x) paste("/projects/caeg/data/production/",sep='/',x,"/stats/reads/fastqc_trim/")),full.names=T)
val <- 0
for(ii in grep("*collapsed_fastqc.zip",bname,val=T)){
       cmd <- paste("unzip -p ",ii, "\"*.txt\"|grep \"Total Seq\"|cut -f2")
val <- val + as.numeric(system(cmd,intern=TRUE))
}
oo <- c(oo,val)
}

oo2 <-c()
for(i in fns[,1]){
  bname <- list.files(sapply(i,function(x) paste("/projects/caeg/data/production/",sep='/',x,"/stats/reads/fastqc_raw/")),full.names=T)
val <- 0
for(ii in grep("*R1_fastqc.zip",bname,val=T)){
       cmd <- paste("unzip -p ",ii, "\"*.txt\"|grep \"Total Seq\"|cut -f2")
val <- val + as.numeric(system(cmd,intern=TRUE))
}
oo2 <- c(oo2,val)
}

oo3<-rbind(oo,oo2)
colnames(oo3)<-sapply(strsplit(fns[,1],"_"),function(x) x[3])
pdf(paste0(plot_outdir,"seqstats_v2.pdf"),width=14)
par(mar=c(11,4,4,4))
 barplot(oo3,las=3,legend=c("Reads retained","Reads sequenced "),col=1:2)
par(mar=c(11,4,4,4))
barplot(oo3[1,]/oo3[2,],las=3,col=3,main="Fraction of retained sequences")


##get all
totalres<-rbind(res,oo3)
barplot(totalres[3,]/totalres[5,],las=3,main="Fraction of classified reads",col=4)
dev.off()
