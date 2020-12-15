pdf("formikkel.pdf")
fil <- "../data/KapK-12-1-35-Ext-12-Lib-12-Index2.col.sorted.sam.gz.species.bdamage.gz.counts.txt"
d<-read.table(fil,as.is=T,he=T,comment.char='!')

fw <- c()
for(i in 0:14)
    fw <- rbind(fw,colSums(d[d$Direction=="5'" & d$Pos==i,-c(1:4)]))

getmismatch <- function(x)
    c(x[1:4]/x[1],x[5:8]/x[6],x[9:12]/x[11],x[13:16]/x[16])

res <- t(apply(fw,1,getmismatch))

stand <- c(1,6,11,16)
res2 <- res[,-stand]

pchs<-rep(0:1,6)
cols <- c(1:9,1:3)
cols <- sample(cols)
cols <- rev(cols)

plot(1:nrow(res2),res2[,1],ylim=range(res2),type='b',pch=pchs[1],col=cols[1],main=fil)
for(i in 2:ncol(res2))
    lines(1:nrow(res2),res2[,i],type='b',pch=pchs[i],col=cols[i])
legend("topright",colnames(res2),pch=pchs,col=cols)


fil <- "../data/KapK-12-1-35-Ext-12-Lib-12-Index2.col.sorted.sam.gz.genus.bdamage.gz.counts.txt"
d<-read.table(fil,as.is=T,he=T,comment.char='!')

fw <- c()
for(i in 0:14)
    fw <- rbind(fw,colSums(d[d$Direction=="5'" & d$Pos==i,-c(1:4)]))

getmismatch <- function(x)
    c(x[1:4]/x[1],x[5:8]/x[6],x[9:12]/x[11],x[13:16]/x[16])

res <- t(apply(fw,1,getmismatch))

stand <- c(1,6,11,16)
res2 <- res[,-stand]

pchs<-rep(0:1,6)
cols <- c(1:9,1:3)
cols <- sample(cols)
cols <- rev(cols)

plot(1:nrow(res2),res2[,1],ylim=range(res2),type='b',pch=pchs[1],col=cols[1],main=fil)
for(i in 2:ncol(res2))
    lines(1:nrow(res2),res2[,i],type='b',pch=pchs[i],col=cols[i])
legend("topright",colnames(res2),pch=pchs,col=cols)


dev.off()
