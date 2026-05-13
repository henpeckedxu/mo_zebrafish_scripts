library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(DESeq2)
library(SeuratWrappers)
library(reticulate)
library(scater)
library(readxl)
library(scCustomize)
library("grid") 
library("gridExtra") 
library("cowplot")
`%notin%` <- Negate(`%in%`)
gex.combined.hicat<-readRDS("./processed_data/snRNA_Analysis/gex.combined.hicat.rds")

#color list for classes
color_list<-c("#DB00FA","#FA0098","#fa0041", "#FF5884",
              "#FF8A00","#ffe000",
              "#116653","#6DCF23","#23CF4E","#23cfa7","#ACF1E1",
              "#6330FF","#318DFF","#31d3ff", "#988298")

#Figure 3a. Dimplot of neuron vs non-neuron
p <- DimPlot(
  subset(gex.combined.hicat, subset=class %notin% c('NA', 'doublet')),
  reduction = "umap.harmony",
  group.by = c("class"),
  label.size = 3,label = FALSE,
  cols = c(rep("#FF5884", 14), "#23cfa7")
)+ggtitle('')+theme(axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.text=element_blank())
p<-p+guides(color = guide_legend(override.aes=list(shape = 15,size=5)))
p<-p+theme(plot.margin = unit(c(0,0,0,0),'mm'))+NoLegend()
p

#Figure 3b. cell proportions of each class
df_plot <- as.data.frame(table(subset(gex.combined.hicat, subset=class %notin% c('NA', 'doublet') & category=='SDA' &orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07"))$class))
colnames(df_plot) <- c('class', 'sda')
df_plot$vda <- as.data.frame(table(subset(gex.combined.hicat, subset=class %notin% c('NA', 'doublet') & category=='VDA')$class))$Freq
df_plot$color <- color_list
df_plot$vda<-df_plot$vda/sum(df_plot$vda)
df_plot$sda<-df_plot$sda/sum(df_plot$sda)
df_plot <- gather(df_plot, category, cell_no, sda:vda, factor_key=TRUE)
plot_ly(df_plot, x = ~category, y = ~cell_no, type = 'bar', width=80,height=880,
        marker=list(color = df_plot$color)) %>% 
  layout(margin = list(l = 0,r = 0,b = 0,t = 0),
         yaxis = list(title = '',showticklabels=FALSE), xaxis = list(categoryarray=~category,title='',showticklabels=FALSE),
         barmode = 'stack')

#Figure 3c. UMAP representation of all cells colored by class
p <- DimPlot(
  subset(gex.combined.hicat, subset=class %notin% c('NA', 'doublet')),
  reduction = "umap.harmony",
  group.by = c("class"),
  label.size = 3,label = FALSE,
  cols = color_list
)+ggtitle('')+theme(axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.text=element_blank())
p<-p+guides(color = guide_legend(override.aes=list(shape = 15,size=5)))
p<-p+theme(plot.margin = unit(c(0,0,0,0),'mm'))+NoLegend()
p
legend <- get_legend(p)
grid.newpage()

##Figure 3d. UMAP representation of all cells colored by no. of DEGs
temp_list<-list()
for (cluster_id in unique(cell_annotation$harmony_clusters_pca_2)){
  temp<-gex.combined.hicat@meta.data[colnames(subset(gex.combined.hicat,subset=harmony_clusters_pca_2==cluster_id[[1]])),]
  if (cluster_id %in%colnames(de_log2fc)){
    temp$deg_no<-colSums(de_log2fc!=0)[[cluster_id]]
    temp$deg_up_no<-colSums(de_log2fc>0)[[cluster_id]]
    temp$deg_down_no<-colSums(de_log2fc<0)[[cluster_id]]
  }
  else{
    temp$deg_no<-0
    temp$deg_up_no<-0
    temp$deg_down_no<-0
  }
  temp_list<-rbind(temp_list,temp)
}
gex.combined.hicat<-AddMetaData(gex.combined.hicat, temp_list)


p <- DimPlot(
  subset(gex.combined.hicat,subset=class %notin% c('NA', 'doublet')),
  reduction = "umap.harmony",
  group.by = c("deg_no"),
  label.size = 5, label = FALSE,
  cols = paletteer_c("grDevices::Viridis", 100)[c(1,sort(unique(subset(gex.combined.hicat,subset=class %notin% c('NA', 'doublet'))$deg_no))[-1]+28)]
)+theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_blank(),
        plot.margin = unit(c(0,0,0,0),'mm'))+NoLegend()+ggtitle("")
p

##ED Figure 7b. UMAP representation of all cells colored by nt type
p <- DimPlot(
  subset(gex.combined.hicat,subset=class %notin% c('doublet')),
  reduction = "umap.harmony",
  group.by = c("nt_type_new"),
  label.size = 5, label = FALSE,
  cols = c(paletteer_c("ggthemes::Classic Red-Blue", 30)[18],paletteer_c("grDevices::Reds 3", 30)[1],
           paletteer_c("grDevices::heat.colors", 30)[seq(1,30,5)],
           paletteer_c("grDevices::Green-Yellow", 30)[seq(1,16,6)], paletteer_c("ggthemes::Gold-Purple Diverging", 30)[c(30)],
           '#CCCCCCFF')
)+ggtitle('')
p+theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

##ED Figure 7c. UMAP representation of all cells colored by anatomy
p <- DimPlot(
  subset(gex.combined.hicat,subset=class %notin% c('doublet')),
  reduction = "umap.harmony",
  group.by = c("anatomy"),
  label.size = 5, label = FALSE,
  cols = c(paletteer_c("ggthemes::Classic Red-Blue", 30)[c(18,21,24)],paletteer_c("ggthemes::Classic Red-Blue", 30)[seq(1,11,2)],
           paletteer_c("grDevices::Green-Yellow", 30)[seq(1,15,3)], paletteer_c("grDevices::Red-Yellow", 30)[c(15,20)], '#CCCCCCFF', 
           paletteer_c("ggthemes::Gold-Purple Diverging", 30)[c(25,30)])
)+ggtitle('')
p+theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p

