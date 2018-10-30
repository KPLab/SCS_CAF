######################################################################################
### 4. ROTS analysis ###
######################################################################################
#Differentually expressed genes are detected using Reproducibility-optimized test statistic (ROTS), for each subgroup compared to the other subgroups.
#Michael Bartoschek, michael.bartoschek@med.lu.se
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


