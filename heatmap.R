install.packages("gplots")
install.packages("RColorBrewer")
install.packages("proxy")
source("http://bioconductor.org/biocLite.R")
biocLite(c("RColorBrewer", "pheatmap"))

setwd("/media/anavarro/Data/Backup-Jena/Aura_Navarro-Quezada/Documents/")
library(gplots)
library(grid)
library(RColorBrewer)
library(pheatmap)
library(proxy)



#function to fill the self-self comparison
FillSim<-function(data=NULL){
  data<-as.matrix(data);
  for (i in 1:length(data[,1])){
    data[i,i]<-1;
  }
  return(data);
}

#functions to do IQR filtering
IQRFilter<-function(data=NULL,cutoff=0.1){
  list.rm<-vector();
  for (i in 1:length(row.names(data))){
    if (IQR(as.numeric(data[i,])) <= cutoff  ){
      list.rm<-c(list.rm,i);
    }
  }
  data<-data[-list.rm,];
  return(data);
} 


library(proxy);

data.all=read.table("cdiffExpr.leaf.P0.01_C0.585.matrix",header=TRUE,row=1)
mat = as.matrix(data.all)
head(mat)
## filter data with IQR
data.all.clean<-IQRFilter(data.all, cutoff=10);
#Scritpt for calculating similarity matrix among different collumns

data.sim<-simil(t(data.all.clean),pairwise=TRUE,upper=TRUE,by_rows = TRUE);
mat=as.matrix(data.sim);
mat<-FillSim(mat);

#write.table(mat,file = 'C:\\Zhihao_Linux\\R_script\\distribution\\PSI_tfpkm\\symmetrical_heatmap\\psi_simi.txt',row.names=F,col.names=F,quote=F,sep="\t")

#png(filename="heatmap/all_ngfpkm.png",width=1000,height=1000)
  pheatmap(mat,
         scale="none",
         border_color = NA,#NA,         
         main = "Isoform percentage",cex.main = 0.1,main_height = 5,#ylab=NULL,#xlab="Tissues",
         dendrogram = "both",
         trace = "none",
         legend = T,
         #legend_breaks = 1:2,
         #legend_labels = c("1e-4", "1e-3","1"),
         cellwidth = 25, cellheight = 0.01, fontsize=10, fontsize_row=10,fontsize_column=13, #adjust the cell parameter
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         #annotation = c("1","2","1","2"),
         treeheight_row = 40,
         treeheight_col = 40,
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         #color = colorRampPalette(c("red","yellow","green"))(30),
         show_rownames = T,
         show_colnames = T,)         
dev.off()


## plot heatmap
#rc <- rainbow(nrow(mat), start=0, end=.5)
#cc <- rainbow(ncol(mat), start=0, end=.5)
#heatmap.2(mat, col = bluered(75),vline=F,hline=F,density.info="non",trace="none",
#          dendrogram="both", scale="none",RowSideColors=rc,ColSideColors=cc,keysize=1.2,
#          distfun = dist,
#          hclustfun = hclust,margins = c(4,4));



#filename ="sysheat_psi.pdf",width = 3, height =4,) 

