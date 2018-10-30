######################################################################################
### 5. Plotting ###
######################################################################################
#Michael Bartoschek, michael.bartoschek@med.lu.se, 2018
#The following functions were used to generate the plots for the manuscript. For connected panels the function multiplot was used.

#scatterplots were generated with the standard plot function in R

#feature plots represent gene expression for each cell on the t-SNE display. It requires the name of the gene as a string, the output of Rtsne and the expression matrix with rownames representing the gene names.
plot.feature2<-function(gene, tsne.output=tsne.out, DATAuse=DATA){
  plot.frame<-data.frame(x=tsne.output$Y[,1], y=tsne.output$Y[,2], log2expr=as.numeric(log2(DATAuse[gene,]+1)))
  
  
  p<-ggplot(plot.frame,aes(x=x, y=y, col=log2expr))+
    geom_point(size=1) +
    labs(title=paste(gene))+
    theme_classic()+
    scale_color_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"), limits=c(0,14))+
    theme(axis.title = element_blank())+
    theme(axis.text = element_blank())+
    theme(axis.line = element_blank())+
    theme(axis.ticks = element_blank())+
    theme(plot.title = element_text(size=20,face="italic"))+
    theme(legend.title  = element_blank())+
    
    
    theme(legend.position = "none")
  
  
  
  
  return(p)
}

plot.feature2("Pdgfra", opt_tsne_full, RPKM.full)




#violin plots represent gene expression for each subpopulation. The color of each violin represents the mean gene expression after log2 transformation.
#gene: Gene name of interest as string. DATAuse: Gene expression matrix with rownames containing gene names. tsne.popus = dbscan output, axis= if FALSE no axis is printed. legend_position= default "none" indicates where legend is plotted. gene_name = if FALSE gene name will not be included in the plot.
plot.violin2 <- function(gene, DATAuse, tsne.popus, axis=FALSE, legend_position="none", gene_name=FALSE){
  
  testframe<-data.frame(expression=as.numeric(DATAuse[paste(gene),]), Population=tsne.popus$cluster)
  testframe$Population <- as.factor(testframe$Population)
  colnames(testframe)<-c("expression", "Population")
  
  col.mean<-vector()
  for(i in levels(testframe$Population)){
    col.mean<-c(col.mean,mean(testframe$expression[which(testframe$Population ==i)]))
  }
  col.mean<-log2(col.mean+1)
  
  col.means<-vector()
  
  for(i in testframe$Population){
    col.means<-c(col.means,col.mean[as.numeric(i)])
  }
  testframe$Mean<-col.means
  testframe$expression<-log2(testframe$expression+1)
  
  
  
  p <- ggplot(testframe, aes(x=Population, y=expression, fill= Mean, color=Mean))+
    geom_violin(scale="width") +
    labs(title=paste(gene), y ="log2(expression)", x="Population")+
    theme_classic() +
    
    
    
    
    scale_color_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"), limits=c(0,14))+
    scale_fill_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"), limits=c(0,14))+
    
    theme(axis.title.y =  element_blank())+
    theme(axis.ticks.y =  element_blank())+
    theme(axis.line.y =   element_blank())+
    theme(axis.text.y =   element_blank())+
    theme(axis.title.x = element_blank())+
    
    
    theme(legend.position=legend_position )
  
  if(axis == FALSE){
    p<-p+
      theme( axis.line.x=element_blank(),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())
    
  }
  
  if(gene_name == FALSE){
    p<-p+  theme(plot.title = element_blank())   
  }else{ p<-p + theme(plot.title = element_text(size=10,face="bold"))}
  
  p
  
}

plot.violin2(gene = "Pdgfra", DATAuse = RPKM.full, tsne.popus = CAFgroups_full)



population_subset<-c(rownames(summary_pop1[summary_pop1$ROTS.statistic<0,])[1:18],rownames(summary_pop2[summary_pop2$ROTS.statistic<0,])[1:18],rownames(summary_pop3[summary_pop3$ROTS.statistic<0,])[1:18],rownames(summary_pop4[summary_pop4$ROTS.statistic<0,])[1:18])
RPKM_heatmap<-RPKM.full[population_subset,]

RPKM_heatmap<-RPKM_heatmap[,order(CAFgroups_full$cluster)]
RPKM_heatmap<-log2(RPKM_heatmap+1)

popul.col<-sort(CAFgroups_full$cluster)
popul.col<-replace(popul.col, popul.col==1,"#1C86EE" )
popul.col<-replace(popul.col, popul.col==2,"#00EE00" )
popul.col<-replace(popul.col, popul.col==3,"#FF9912" )
popul.col<-replace(popul.col, popul.col==4,"#FF3E96" )
library(gplots)

pdf("heatmap_genes_population.pdf")
heatmap.2(as.matrix(RPKM_heatmap),ColSideColors = as.character(popul.col), tracecol = NA, dendrogram = "none",col=bluered, labCol = FALSE, scale="none", key = TRUE, symkey = F, symm=F,  key.xlab = "", key.ylab = "", density.info = "density", key.title = "log2(RPKM+1)", keysize = 1.2, denscol="black", Colv=FALSE)
dev.off()

####### SC3 Heatmaps #########
library(scater)
library(SC3)

anno<-data.frame(Plate = unlist(lapply(strsplit(colnames(RPKM),"_"),function(x) x[3])), Population = CAFgroups_full$cluster)
rownames(anno)<-colnames(RPKM)
pd <- new("AnnotatedDataFrame", data = anno)
sceset <-  SingleCellExperiment(assays = list(normcounts = as.matrix(RPKM)) ,colData =anno)

counts(sceset)<-normcounts(sceset)
logcounts(sceset)<- log2(normcounts(sceset)+1)
rowData(sceset)$feature_symbol <- rownames(sceset)
sceset<-sceset[!duplicated(rownames(sceset)),]
sceset <- sc3(sceset, ks = 3:5, biology = TRUE)

pdf("sc3_heatmap.pdf", useDingbats = FALSE, width=6, height = 4.5)

sc3_plot_consensus(
  sceset, k=4,
  show_pdata = c(
    "Plate",
    "Population"
  )
)
dev.off()



matrisome<-as.character(read.table("GitHub/170210_naba_matrisome.txt")[,1])
matrisome<-lookuptable[lookuptable[,1] %in% matrisome,4]

sce.matrisome <- SingleCellExperiment(assays = list(normcounts = as.matrix(RPKM.full[matrisome,])), colData = anno)
counts(sce.matrisome)<-normcounts(sce.matrisome)
logcounts(sce.matrisome)<- log2(normcounts(sce.matrisome)+1)
rowData(sce.matrisome)$feature_symbol <- rownames(sce.matrisome)
sceset<-sce.matrisome[!duplicated(rownames(sce.matrisome)),]


sce.matrisome <- sc3(sce.matrisome, ks = 4, biology = TRUE)

pdf("sce_matrisome_expression.pdf",useDingbats = FALSE, width=6, height = 4.5)
sc3_plot_expression(
  sce.matrisome, k=4,
  show_pdata = c(
    "Plate",
    "Population"
  )
)
dev.off()
