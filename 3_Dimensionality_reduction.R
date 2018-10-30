######################################################################################
### 3. PCA and t-SNE ###
######################################################################################

#created by Michael Bartoschek, michael.bartoschek@med.lu.se, 2018
#run t-SNE 50 times and select the optimal run based on the itercosts parameter
#subgroups of cells are assigned by dbscan
library(Rtsne)
N_tsne <- 50
tsne_out <- list(length = N_tsne)
KL <- vector(length = N_tsne)
set.seed(1234)
for(k in 1:N_tsne)
{
  tsne_out[[k]]<-Rtsne(t(log10(RPKM+1)),initial_dims=30,verbose=FALSE,check_duplicates=FALSE,
                       perplexity=27, dims=2,max_iter=5000)
  KL[k]<-tail(tsne_out[[k]]$itercosts,1)
  print(paste0("FINISHED ",k," TSNE ITERATION"))
}
names(KL) <- c(1:N_tsne)
opt_tsne <- tsne_out[[as.numeric(names(KL)[KL==min(KL)])]]$Y
opt_tsne_full<-tsne_out[[as.numeric(names(KL)[KL==min(KL)])]]


library(dbscan)
plot(opt_tsne,  col=dbscan(opt_tsne,eps=3.1)$cluster, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
CAFgroups<-dbscan(opt_tsne,eps=3.1)$cluster
CAFgroups_full<-dbscan(opt_tsne,eps=3.1)
CAFgroups[CAFgroups==0]<-1
CAFgroups_full$cluster[CAFgroups_full$cluster==0]<-1
plot(opt_tsne,  col=CAFgroups, pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")

RPKM.PCA<-prcomp(log2(t(RPKM)+1), center=TRUE)
plot(RPKM.PCA$x,main="first PCA", pch=19, col=CAFgroups)