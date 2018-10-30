######################################################################################
### 2. Data processing using ERCC spike ins ###
######################################################################################
#ENDOGENOUS GENES
# created by Nikolay Oskolkov, nikolay.oskolkov@scilifelabs.se, 2017
Path_Main<-"/CAF_SCS/"
plate1_raw<-read.delim(paste(Path_Main,"/SS2_15_0048/counts.tab",sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate1_raw$gene<-make.unique(as.character(plate1_raw$gene))
colnames(plate1_raw)[2:length(colnames(plate1_raw))]<-paste("SS2_15_0048_",colnames(plate1_raw)[2:length(colnames(plate1_raw))],sep="")
plate2_raw<-read.delim(paste(Path_Main,"/SS2_15_0049/counts.tab",sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate2_raw$gene<-make.unique(as.character(plate2_raw$gene))
colnames(plate2_raw)[2:length(colnames(plate2_raw))]<-paste("SS2_15_0049_",colnames(plate2_raw)[2:length(colnames(plate2_raw))],sep="")
expr_raw<-merge(plate1_raw,plate2_raw,by="gene",all=TRUE)
rownames(expr_raw)<-as.character(expr_raw$gene)
expr_raw$gene<-NULL
sum(expr_raw==0)/(dim(expr_raw)[1]*dim(expr_raw)[2])
#SPIKES
plate1_raw_ercc<-read.delim(paste(Path_Main,"/SS2_15_0048/counts-ercc.tab",sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate1_raw_ercc$gene<-make.unique(as.character(plate1_raw_ercc$gene))
colnames(plate1_raw_ercc)[2:length(colnames(plate1_raw_ercc))]<-paste("SS2_15_0048_",colnames(plate1_raw_ercc)[2:length(colnames(plate1_raw_ercc))],sep="")
plate2_raw_ercc<-read.delim(paste(Path_Main,"/SS2_15_0049/counts-ercc.tab",sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate2_raw_ercc$gene<-make.unique(as.character(plate2_raw_ercc$gene))
colnames(plate2_raw_ercc)[2:length(colnames(plate2_raw_ercc))]<-paste("SS2_15_0049_",colnames(plate2_raw_ercc)[2:length(colnames(plate2_raw_ercc))],sep="")
expr_raw_ercc<-merge(plate1_raw_ercc,plate2_raw_ercc,by="gene",all=TRUE)
rownames(expr_raw_ercc)<-as.character(expr_raw_ercc$gene)
expr_raw_ercc$gene<-NULL
sum(expr_raw_ercc==0)/(dim(expr_raw_ercc)[1]*dim(expr_raw_ercc)[2])
barplot(sort(as.numeric(colSums(expr_raw_ercc)),decreasing=TRUE),ylab="SPIKE LIBRARY SIZE",xlab="CELL INDEX")
hist(as.numeric(colSums(expr_raw_ercc)),col="brown",main="Distribution of Spike Library Sizes",xlab="Spike Library Size",breaks=20)
#MERGING ENDOGENOUS GENES WITH SPIKES
print(paste0("There are ",dim(expr_raw)[1]," endogenous genes"))
print(paste0("There are ",dim(expr_raw_ercc)[1]," spikes"))
all.counts.raw<-rbind(expr_raw,expr_raw_ercc)
sum(all.counts.raw==0)/(dim(all.counts.raw)[1]*dim(all.counts.raw)[2])
dim(all.counts.raw[rowSums(all.counts.raw)==0,])

cell_QC<-read.delim(paste(Path_Main,"/qc/qc_2plates.filtered_cells.txt",sep=""),row.names=1,header=TRUE,sep="\t")
rownames(cell_QC)<-gsub("__","_",rownames(cell_QC))
all.counts.raw<-subset(all.counts.raw,select=colnames(all.counts.raw)[!colnames(all.counts.raw)%in%rownames(cell_QC)])
expr_raw<-subset(expr_raw,select=colnames(expr_raw)[!colnames(expr_raw)%in%rownames(cell_QC)])
expr_raw_ercc<-subset(expr_raw_ercc,select=colnames(expr_raw_ercc)[!colnames(expr_raw_ercc)%in%rownames(cell_QC)])
#FILTER OUT GENES HAVING LESS THAN 1 COUNT IN AVERAGE OVER ALL CELLS
all.counts.raw<-all.counts.raw[rowMeans(all.counts.raw)>0,]
expr_raw<-expr_raw[rowMeans(expr_raw)>=1,]
expr_raw_ercc<-expr_raw_ercc[rowMeans(expr_raw_ercc)>0,]
#PLOT CV^2 vs. MEAN RAW COUNT
library("matrixStats")
mean_expr_raw<-as.numeric(rowMeans(expr_raw,na.rm=TRUE))
sd_expr_raw<-rowSds(as.matrix(expr_raw),na.rm=TRUE)
cv_squared_expr_raw<-(sd_expr_raw/mean_expr_raw)^2
plot(log10(cv_squared_expr_raw)~log10(mean_expr_raw),pch=20,cex=0.5,xlab="log10 ( mean raw count )",ylab="log10 ( CV² )",main="RAW COUNTS")
mean_expr_raw_ercc<-as.numeric(rowMeans(expr_raw_ercc,na.rm=TRUE))
sd_expr_raw_ercc<-rowSds(as.matrix(expr_raw_ercc),na.rm=TRUE)
cv_squared_expr_raw_ercc<-(sd_expr_raw_ercc/mean_expr_raw_ercc)^2
points(log10(cv_squared_expr_raw_ercc)~log10(mean_expr_raw_ercc),col="red",pch=20,cex=1.5)
#FIT SPIKEINS WITH A CURVE
fit_expr_raw_ercc<-loess(log10(cv_squared_expr_raw_ercc)[is.finite(log10(mean_expr_raw_ercc))]~log10(mean_expr_raw_ercc)[is.finite(log10(mean_expr_raw_ercc))],span=1)
j<-order(log10(mean_expr_raw_ercc)[is.finite(log10(mean_expr_raw_ercc))])
lines(fit_expr_raw_ercc$fitted[j]~log10(mean_expr_raw_ercc)[is.finite(log10(mean_expr_raw_ercc))][j],col="red",lwd=3)
pred_expr_raw<-predict(fit_expr_raw_ercc,log10(mean_expr_raw))
#DETERMINE VARIABLE GENES THAT ARE ABOVE THE SPIKEINS CURVE
filtered_expr_raw<-expr_raw[log10(cv_squared_expr_raw)>=pred_expr_raw,]
filtered_expr_raw<-filtered_expr_raw[grepl("NA",rownames(filtered_expr_raw))==FALSE,]



plate1_rpkm<-read.delim(paste(Path_Main,"/SS2_15_0048/rpkms.tab", sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate1_rpkm$gene<-make.unique(as.character(plate1_rpkm$gene))
colnames(plate1_rpkm)[2:length(colnames(plate1_rpkm))]<-paste("SS2_15_0048_",colnames(plate1_rpkm)[2:length(colnames(plate1_rpkm))],sep="")
plate2_rpkm<-read.delim(paste(Path_Main,"/SS2_15_0049/rpkms.tab", sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate2_rpkm$gene<-make.unique(as.character(plate2_rpkm$gene))
colnames(plate2_rpkm)[2:length(colnames(plate2_rpkm))]<-paste("SS2_15_0049_",colnames(plate2_rpkm)[2:length(colnames(plate2_rpkm))],sep="")
expr_rpkm<-merge(plate1_rpkm,plate2_rpkm,by="gene",all=TRUE)
rownames(expr_rpkm)<-as.character(expr_rpkm$gene)
expr_rpkm$gene<-NULL
sum(expr_rpkm==0)/(dim(expr_rpkm)[1]*dim(expr_rpkm)[2])
#SPIKES
plate1_rpkm_ercc<-read.delim(paste(Path_Main,"/SS2_15_0048/rpkms-ercc.tab", sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate1_rpkm_ercc$gene<-make.unique(as.character(plate1_rpkm_ercc$gene))
colnames(plate1_rpkm_ercc)[2:length(colnames(plate1_rpkm_ercc))]<-paste("SS2_15_0048_",colnames(plate1_rpkm_ercc)[2:length(colnames(plate1_rpkm_ercc))],sep="")
plate2_rpkm_ercc<-read.delim(paste(Path_Main,"/SS2_15_0049/rpkms-ercc.tab", sep=""),header=TRUE,check.names=FALSE,sep="\t")
plate2_rpkm_ercc$gene<-make.unique(as.character(plate2_rpkm_ercc$gene))
colnames(plate2_rpkm_ercc)[2:length(colnames(plate2_rpkm_ercc))]<-paste("SS2_15_0049_",colnames(plate2_rpkm_ercc)[2:length(colnames(plate2_rpkm_ercc))],sep="")
expr_rpkm_ercc<-merge(plate1_rpkm_ercc,plate2_rpkm_ercc,by="gene",all=TRUE)
rownames(expr_rpkm_ercc)<-as.character(expr_rpkm_ercc$gene)
expr_rpkm_ercc$gene<-NULL
#MERGING ENDOGENOUS GENES WITH SPIKES
print(paste0("There are ",dim(expr_rpkm)[1]," endogenous genes"))
print(paste0("There are ",dim(expr_rpkm_ercc)[1]," spikes"))
all.counts.rpkm<-rbind(expr_rpkm,expr_rpkm_ercc)
sum(all.counts.rpkm==0)/(dim(all.counts.rpkm)[1]*dim(all.counts.rpkm)[2])
dim(all.counts.rpkm[rowSums(all.counts.rpkm)==0,])
#PLOT BARPLOT OF TOP EXPRESSED GENES
mean_rpkms<-data.frame(Gene=rownames(all.counts.rpkm),Mean_rpkm=rowMedians(as.matrix(all.counts.rpkm)))
mean_rpkms<-mean_rpkms[order(-mean_rpkms$Mean_rpkm),]
head(mean_rpkms,20)
par(mfrow=c(1,1))
barplot(log10(as.numeric(mean_rpkms$Mean_rpkm)[1:50]),names.arg=mean_rpkms$Gene[1:50],las=2,ylim=c(0,5),cex.names=0.8,ylab="Median RPKM")

#FILTER POOR CELLS OUT
cell_QC<-read.delim(paste(Path_Main,"/qc/qc_2plates.filtered_cells.txt",sep=""),row.names=1,header=TRUE,sep="\t")
rownames(cell_QC)<-gsub("__","_",rownames(cell_QC))
all.counts.rpkm<-subset(all.counts.rpkm,select=colnames(all.counts.rpkm)[!colnames(all.counts.rpkm)%in%rownames(cell_QC)])
expr_rpkm<-subset(expr_rpkm,select=colnames(expr_rpkm)[!colnames(expr_rpkm)%in%rownames(cell_QC)])
expr_rpkm_ercc<-subset(expr_rpkm_ercc,select=colnames(expr_rpkm_ercc)[!colnames(expr_rpkm_ercc)%in%rownames(cell_QC)])

#FILTER OUT GENES HAVING LESS THAN 1 RPKM COUNT IN AVERAGE OVER ALL CELLS
RPKM.full<-all.counts.rpkm[rowSums(all.counts.rpkm)>0,]
all.counts.rpkm<-all.counts.rpkm[rowMeans(all.counts.rpkm)>=1,]
expr_rpkm<-expr_rpkm[rowMeans(expr_rpkm)>=1,]
expr_rpkm_ercc<-expr_rpkm_ercc[rowMeans(expr_rpkm_ercc)>=1,]


#PLOT CV^2 vs. MEAN RPKM
library("matrixStats")
mean_expr_rpkm<-as.numeric(rowMeans(expr_rpkm,na.rm=TRUE))
sd_expr_rpkm<-rowSds(as.matrix(expr_rpkm),na.rm=TRUE)
cv_squared_expr_rpkm<-(sd_expr_rpkm/mean_expr_rpkm)^2
plot(log10(cv_squared_expr_rpkm)~log10(mean_expr_rpkm),pch=20,cex=0.5,xlab="log10 ( mean rpkm count )",ylab="log10 ( CV² )",main="RPKM COUNTS")
mean_expr_rpkm_ercc<-as.numeric(rowMeans(expr_rpkm_ercc,na.rm=TRUE))
sd_expr_rpkm_ercc<-rowSds(as.matrix(expr_rpkm_ercc),na.rm=TRUE)
cv_squared_expr_rpkm_ercc<-(sd_expr_rpkm_ercc/mean_expr_rpkm_ercc)^2
points(log10(cv_squared_expr_rpkm_ercc)~log10(mean_expr_rpkm_ercc),col="red",pch=20,cex=1.5)

#FIT SPIKEINS WITH A CURVE
fit_expr_rpkm_ercc<-loess(log10(cv_squared_expr_rpkm_ercc)[is.finite(log10(mean_expr_rpkm_ercc))]~log10(mean_expr_rpkm_ercc)[is.finite(log10(mean_expr_rpkm_ercc))],span=1)
j<-order(log10(mean_expr_rpkm_ercc)[is.finite(log10(mean_expr_rpkm_ercc))])
lines(fit_expr_rpkm_ercc$fitted[j]~log10(mean_expr_rpkm_ercc)[is.finite(log10(mean_expr_rpkm_ercc))][j],col="red",lwd=3)
pred_expr_rpkm<-predict(fit_expr_rpkm_ercc,log10(mean_expr_rpkm))
#DETERMINE VARIABLE GENES THAT ARE ABOVE THE SPIKEINS CURVE
filtered_expr_rpkm<-expr_rpkm[log10(cv_squared_expr_rpkm)>=pred_expr_rpkm,]
RPKM<-filtered_expr_rpkm[grepl("NA",rownames(filtered_expr_rpkm))==FALSE,]
sum(RPKM==0)/(dim(RPKM)[1]*dim(RPKM)[2])
