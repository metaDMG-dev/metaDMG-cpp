########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l)
    unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
    if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
        cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
        q("no")
    }
    arguments<-list()
    for(a in names(args))
        arguments[[a]]<-getArgs(a,l)
    
    if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
        cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
        q("no")
    }
    for(a in names(args))
        if(is.null(arguments[[a]]))
            arguments[[a]]<-args[[match(a,names(args))]]
    arguments
}

print.args<-function(args,des){
    if(missing(des)){
        des<-as.list(rep("",length(args)))
        names(des)<-names(args)
    }
    cat("->  Needed arguments:\n")
    mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
    cat("->  Optional arguments (defaults):\n")
    mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
    q("no")
}

###### ####### ###### ###### ###### #######
# choose your parameters and defaults
# NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(file=NULL,outfile=NA,type=1)
#if no argument are given prints the need arguments and the optional ones with default
des<-list(file=" The damage file, typicaly called meta.res.gz",outfile="name of output file",type="type=1 -r 1,type=1 -r 0 in ./metadamage getdamage [maybe not used here yet]")

######################################
#######get arguments and add to workspace
### do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}

cat("\t-> Using: file: ",file," outfile: ",outfile)
if(is.na(outfile))
    outfile <- paste0(file,".pdf")
cat("\t will write: ",outfile,"\n")
pdf(outfile,width=14)
 par(mfrow=c(1,2))
norm  <- function(x) x/sum(x)

getdam <- function(x){
    as.numeric(apply(matrix(x,4),2,function(x) norm(as.numeric(x))))[-c(1,6,11,16)]
}

  cols <- c(1:9,1:3)
    cols <- sample(cols)
    cols <- rev(cols)
   
d<-read.table(file,as.is=T)

for(j in 1:nrow(d)){
    if(nrow(d)>1)
        dd<-d[j,-c(1,2)]
    else
        dd<-d[j,-c(1)]
    howmany<-length(dd)/16/2
    
    
    dam5 <- apply(matrix(dd[1:(howmany*16)],16),2,getdam)
    dam3 <- apply(matrix(dd[-c(1:(howmany*16))],16),2,getdam)
    dam3 <- dam3[,ncol(dam3):1]
    nam <- c("A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G")
    rownames(dam5)<-nam
    rownames(dam3)<-nam
    pchs<-rep(0:1,6)
    
    
    
  
    cap <- ""
    if(nrow(d)>1)
        cap <- d[j,1]
    cap <- paste0(cap,": ","5'")
    if(any(is.na(dam5))|any(is.na(dam3)))
        next
    plot(1:ncol(dam5),dam5[1,],ylim=c(0,max(c(dam5,dam3))),type='b',lwd=2,col=cols[1],xlab="Position",ylab="Misincorporation",pch=pchs[1],main=cap)
    for(i in 2:nrow(dam5))
        lines(1:ncol(dam5),dam5[i,],col=cols[i],lwd=2,pch=pchs[i],type='b')
    
    legend("topright",nam,pch=pchs,col=1:9)
    cap <- ""
    if(nrow(d)>1)
        cap <- d[j,1]
    cap <- paste0(cap,": ","3'")
    plot(1:ncol(dam3),dam3[1,],ylim=c(0,max(c(dam5,dam3))),type='b',lwd=2,col=cols[1],xlab="Position",ylab="Misincorporation",pch=pchs[1],main=cap)
    for(i in 2:nrow(dam3))
        lines(1:ncol(dam3),dam3[i,],col=cols[i],lwd=2,pch=pchs[i],type='b')
    legend("topleft",nam,pch=pchs,col=1:9)
}
dev.off()
