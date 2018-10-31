######################################################################################
### 4. ROTS analysis ###
######################################################################################
#Differentually expressed genes are detected using Reproducibility-optimized test statistic (ROTS), for each subgroup compared to the other subgroups.
#Michael Bartoschek, michael.bartoschek@med.lu.se, and Nikolay Oskolkov, nikolay.oskolkov@scilifelab.se
library(ROTS)
library(plyr)
groups<-CAFgroups
groups[groups!=1]<-234

ROTS_input<-RPKM.full[rowMeans(RPKM.full)>=1,]
ROTS_input<-as.matrix(log2(ROTS_input+1))
results_pop1 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop1<-data.frame(summary(results_pop1, fdr=1))

groups<-CAFgroups
groups[groups!=2]<-134
results_pop2 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop2<-data.frame(summary(results_pop2, fdr=1))


groups<-CAFgroups
groups[groups!=3]<-124
results_pop3 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop3<-data.frame(summary(results_pop3, fdr=1))


groups<-CAFgroups
groups[groups!=4]<-123
results_pop4 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop4<-data.frame(summary(results_pop4, fdr=1))


####### DESeq2 ###########
library("scran")
library("limSolve")
library(scater)
library(DESeq2)
ann<-data.frame(Plate = factor(unlist(lapply(strsplit(colnames(RPKM.full),"_"),function(x) x[3]))), Population = factor(gsub("(3|4)","2",as.character(CAFgroups)),levels=c("1","2")))
ann<-data.frame(Population = factor(gsub("(3|4)","2",as.character(CAFgroups)),levels=c("1","2")))
rownames(ann)<-colnames(RPKM.full)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = all.counts.raw[rownames(RPKM.full),],
  colData = ann,
  design = ~ Population)

ddsFullCountTable<-DESeq(ddsFullCountTable)
DESeq_result<-results(ddsFullCountTable)
DESeq_result<-DESeq_result[order(DESeq_result$padj, DESeq_result$pvalue),]
head(DESeq_result,30)
write.table(DESeq_result[grep("ERCC", rownames(DESeq_result), invert=TRUE),], "DESeq_result.txt", col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")

                                             
####### EdgeR ############
library("edgeR")
edgeR_Data<-DGEList(counts=all.counts.raw[rownames(RPKM.full),], group=ann$Population)
edgeR_Data<-estimateCommonDisp(edgeR_Data)
edgeR_Data<-estimateTagwiseDisp(edgeR_Data)
edgeR_result<-exactTest(edgeR_Data)
edgeR_result_table<-edgeR_result$table
edgeR_result_table<-edgeR_result_table[order(edgeR_result_table$PValue),]
edgeR_result_table<-edgeR_result_table[grep("ERCC", rownames(edgeR_result_table), invert=TRUE),]
write.table(edgeR_result_table, "EdgeR_result.txt", col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")

###### Wilcox ##########

NumPerm<-1000
POP1_expr<-subset(RPKM.full,select=rownames(ann)[ann==1])
POP2_expr<-subset(RPKM.full,select=rownames(ann)[ann==2])
p_wilcox<-vector()
p_t<-vector()
p_perm<-vector()
statistics<-vector()
median_POP1_expr<-vector()
median_POP2_expr<-vector()
a<-seq(from=0,to=length(rownames(RPKM.full)),by=1000)



print("START DIFFERENTIAL GENE EXPRESSION BETWEEN POP1 AND POP2")
for(i in 1:length(rownames(RPKM.full)))
{
  p_wilcox<-append(p_wilcox,wilcox.test(as.numeric(POP1_expr[rownames(RPKM.full)[i],]),as.numeric(POP2_expr[rownames(RPKM.full)[i],]))$p.value)
  statistics<-append(statistics,wilcox.test(as.numeric(POP1_expr[rownames(RPKM.full)[i],]),as.numeric(POP2_expr[rownames(RPKM.full)[i],]))$statistic)
  p_t<-append(p_t,t.test(as.numeric(POP1_expr[rownames(RPKM.full)[i],]),as.numeric(POP2_expr[rownames(RPKM.full)[i],]))$p.value)
  p_perm<-append(p_perm,PermTest_Median(as.numeric(POP1_expr[rownames(RPKM.full)[i],]),as.numeric(POP2_expr[rownames(RPKM.full)[i],]),NumPerm))
  median_POP1_expr<-append(median_POP1_expr,median(as.numeric(POP1_expr[rownames(RPKM.full)[i],])))
  median_POP2_expr<-append(median_POP2_expr,median(as.numeric(POP2_expr[rownames(RPKM.full)[i],])))
  if(i%in%a){print(paste("FINISHED ",i," GENES",sep=""))}
}
fold_change<-median_POP1_expr/median_POP2_expr
log2_fold_change<-log2(fold_change)
p_adj<-p.adjust(p_wilcox,method="fdr")
output_wilcox<-data.frame(GENE=rownames(RPKM.full),POP1_EXPR=median_POP1_expr,POP234_EXPR=median_POP2_expr,FOLD_CHANGE=fold_change,LOG2FC=log2_fold_change,WILCOX_STAT=statistics,P_T_TEST=p_t,P_PERM=p_perm,P_WILCOX=p_wilcox,FDR=p_adj)
output_wilcox<-output_wilcox[order(output_wilcox$P_PERM,output_wilcox$FDR,output_wilcox$P_WILCOX,output_wilcox$P_T_TEST,-abs(output_wilcox$LOG2FC)),]
print(head(output_wilcox,20))
write.table(output_wilcox,file="Wilcox_Perm_de_results.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
                                             
#SCDE WORKFLOW
library("scde")
# factor determining cell types
sg<-factor(as.numeric(CAFgroups))
# the group factor should be named accordingly
names(sg)<-colnames(expr_raw)
table(sg)
# define two groups of cells
#groups<-sg
groups <- factor(gsub("(3|4)","2",as.character(sg)),levels=c("1","2"))
table(groups)
# calculate models

cd<-apply(expr_raw,2,function(x) {storage.mode(x) <- 'integer'; x})
colnames(cd)<-colnames(expr_raw)
o.ifm<-scde.error.models(counts=cd,groups=groups,n.cores=4,threshold.segmentation=TRUE,save.crossfit.plots=FALSE,save.model.plots=FALSE,verbose=1)
print(head(o.ifm))
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells<-o.ifm$corr.a > 0
table(valid.cells)
o.ifm<-o.ifm[valid.cells, ]
# estimate gene expression prior
o.prior<-scde.expression.prior(models=o.ifm,counts=cd,length.out=400,show.plot=FALSE)
# run differential expression tests on all genes.
ediff<-scde.expression.difference(o.ifm,cd,o.prior,groups=groups,n.randomizations=100,n.cores=4,verbose=1) #batch=batch
# top upregulated genes
ediff_order<-ediff[order(abs(ediff$Z),decreasing=TRUE), ]
head(ediff_order,20)                                           
write.table(ediff_order,file="scde_de_results_1_vs_234.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
