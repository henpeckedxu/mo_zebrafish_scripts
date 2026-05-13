library(sceasy)
library(reticulate)
library(S4Vectors)
library(Matrix)
library(Matrix.utils)
library(pheatmap)
library(tibble)
library(tidyverse)
library(readxl)
library(scales)
library(paletteer)
library(dittoSeq)
library(DESeq2)
library(plotly)
#####################################################################
####################### Pre-test for DEG analysis ###################
#####################################################################
##Step 1. convert Seurat object to Single Cell Experiment object (sce)
gex.combined.hicat <- readRDS('./processed_data/gex.combined.hicat.rds')
sce<-subset(gex.combined.hicat,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07"))
sce<-subset(gex.temp)

sce$sample.category<-paste(sce$category, sce$orig.ident, sep = "_")
sce<-convertFormat(sce, from = "seurat", to="sce")

##Step 2. 
###1. extract count matrix with each rownames as gene names and colnames as cell barcode
cluster_level = 'harmony_clusters_pca_2'
cluster_level = 'temp_class'
cluster_names <- levels(droplevels(as.factor(colData(sce)[[cluster_level]])))
sample_names<-levels(as.factor(colData(sce)$sample.category))
groups <- colData(sce)[, c(cluster_level, "sample.category")]
groups[[cluster_level]]<-as.factor(groups[[cluster_level]])
groups$sample.category<-as.factor(groups$sample.category)

###2. generate an aggregated matrix that sums cells from a same sample and a same cluster
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
aggr_counts <- t(aggr_counts)
counts_ls<-list()
for (i in 1:length(cluster_names)){
  counts_ls[[i]]<-aggr_counts[,substr(colnames(aggr_counts), 1, nchar(colnames(aggr_counts))-12)==cluster_names[i]]
  names(counts_ls)[i]<-cluster_names[i]
}

##Step 3. prepare metadata that matches the aggregated matrix generated in step 2
metadata<-colData(sce)%>%
  as.data.frame()%>%
  dplyr::select(sample.category,category)
metadata<-metadata[!duplicated(metadata), ]
rownames(metadata) <- metadata$sample.category
t<-table(colData(sce)$sample.category,
         colData(sce)[[cluster_level]])

metadata_ls<-list()
for (i in 1:length(counts_ls)) {
  df<-data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  df[[cluster_level]]<-substr(df$cluster_sample_id,1,nchar(df$cluster_sample_id)-12)
  df$sample.category<-substr(df$cluster_sample_id,nchar(df$cluster_sample_id)-10,nchar(df$cluster_sample_id))
  cell_counts<-t[,colnames(t)==unique(df[[cluster_level]])]
  cell_counts <- cell_counts[cell_counts > 0]
  sample_order <- match(df$sample.category, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  df$cell_count <- cell_counts
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  rownames(df) <- df$cluster_sample_id
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df[[cluster_level]])
}

##Step 4. perform DESeq for each cluster
###create a function that takes cluster id as an input###
DESeq_analysis<-function(cluster){
  idx <- which(names(counts_ls) == cluster)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ category)
  rld <- rlog(dds, blind=TRUE)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds,name = "category_VDA_vs_SDA",alpha = 0.05)
  res <- lfcShrink(dds, 
                   res=res,
                   coef = "category_VDA_vs_SDA",
                   type = "apeglm")
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj)
  sig_res <- dplyr::filter(res_tbl, padj < 0.05) %>%
    dplyr::arrange(padj)
  sig_res$cluster<-cluster
  return (sig_res)
}

DEseq_res_ls<-list()
for(i in 1:length(cluster_names)){
  DEseq_res_ls[[i]]<-DESeq_analysis(cluster_names[i])
  print(i)
}
sv_deseq2_6s <- bind_rows(DEseq_res_ls)
colnames(sv_deseq2_6s)[7]=cluster_level
sv_deseq2_6s <- merge(sv_deseq2_6s, cell_annotation,by=cluster_level)
sv_deseq2_6s_flt<-sv_deseq2_6s[abs(sv_deseq2_6s$log2FoldChange)>1,]
gwas_gene<-read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Crispr/data/psy_genes_for_KO/20240617_GWAS_loci_and_candidate_protein-coding_genes.xlsx', 
                      sheet = 'Candidate_gene_human_ortho', range = "A5:AD342",col_names=FALSE)
sv_deseq2_6s_flt$gwas<-ifelse(paste0("ENSDARG", str_split_i(sv_deseq2_6s_flt$gene, "-ENSDARG",2)) %in% gwas_gene[[1]], 'y', 'n')

saveRDS(sv_deseq2_6s, file = "./DESeq2/sv_deseq2_6s_raw_20250501.rds")
saveRDS(sv_deseq2_6s_flt, file = "./DESeq2/sv_deseq2_6s_flt_20250501.rds")


#####################################################################
####################### Part 1. DEG analysis #######################
#####################################################################
####develop a function to perform DEG analysis based on DESeq2####
DESeq2_custom <- function(seurat_obj, cluster_level){
  print('seurat object imported as below: ')
  print(seurat_obj)
  ##Step 1. 
  ##exclude clusters that are not represented by every sample
  bad_cl <- c()
  for (cl in names(table(seurat_obj[[cluster_level]]))){
    if (length(unique(subset(seurat_obj, subset=!!sym(cluster_level)==cl)$orig.ident))<6){
      bad_cl <- c(bad_cl, cl)
    }
  }
  print(bad_cl)
  seurat_obj <- subset(seurat_obj, subset = !!sym(cluster_level) %notin% bad_cl)
  seurat_obj$sample.category<-paste(seurat_obj$category, seurat_obj$orig.ident, sep = "_")
  
  sce<-convertFormat(seurat_obj, from = "seurat", to="sce")
  print('object created as below')
  print(seurat_obj)
  ##Step 2. 
  ###1. extract count matrix with each rownames as gene names and colnames as cell barcode
  
  cluster_level = cluster_level
  cluster_names <- levels(droplevels(as.factor(colData(sce)[[cluster_level]])))
  sample_names<-levels(as.factor(colData(sce)$sample.category))
  groups <- colData(sce)[, c(cluster_level, "sample.category")]
  groups[[cluster_level]]<-as.factor(groups[[cluster_level]])
  groups$sample.category<-as.factor(groups$sample.category)
  
  
  ###2. generate an aggregated matrix that sums cells from a same sample and a same cluster
  aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
  aggr_counts <- t(aggr_counts)
  counts_ls<-list()
  for (i in 1:length(cluster_names)){
    counts_ls[[i]]<-aggr_counts[,substr(colnames(aggr_counts), 1, nchar(colnames(aggr_counts))-12)==cluster_names[i]]
    names(counts_ls)[i]<-cluster_names[i]
  }
  
  ##Step 3. prepare metadata that matches the aggregated matrix generated in step 2
  metadata<-colData(sce)%>%
    as.data.frame()%>%
    dplyr::select(sample.category,category)
  metadata<-metadata[!duplicated(metadata), ]
  rownames(metadata) <- metadata$sample.category
  t<-table(colData(sce)$sample.category,
           colData(sce)[[cluster_level]])
  
  metadata_ls<-list()
  for (i in 1:length(counts_ls)) {
    df<-data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    df[[cluster_level]]<-substr(df$cluster_sample_id,1,nchar(df$cluster_sample_id)-12)
    df$sample.category<-substr(df$cluster_sample_id,nchar(df$cluster_sample_id)-10,nchar(df$cluster_sample_id))
    cell_counts<-t[,colnames(t)==unique(df[[cluster_level]])]
    cell_counts <- cell_counts[cell_counts > 0]
    sample_order <- match(df$sample.category, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    df$cell_count <- cell_counts
    df <- plyr::join(df, metadata, 
                     by = intersect(names(df), names(metadata)))
    rownames(df) <- df$cluster_sample_id
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df[[cluster_level]])
  }
  
  ##Step 4. perform DESeq for each cluster
  ###create a function that takes cluster id as an input###
  DESeq_analysis<-function(cluster){
    idx <- which(names(counts_ls) == cluster)
    cluster_counts <- counts_ls[[idx]]
    cluster_metadata <- metadata_ls[[idx]]
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ category)
    rld <- rlog(dds, blind=TRUE)
    
    dds <- DESeq(dds)
    resultsNames(dds)
    res <- results(dds,name = "category_VDA_vs_SDA",alpha = 0.05)
    res <- lfcShrink(dds, 
                     res=res,
                     coef = "category_VDA_vs_SDA",
                     type = "apeglm")
    res_tbl <- res %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(padj)
    sig_res <- dplyr::filter(res_tbl, padj < 0.05) %>%
      dplyr::arrange(padj)
    sig_res$cluster<-cluster
    return (sig_res)
  }
  
  DEseq_res_ls<-list()
 
  for(i in 1:length(cluster_names)){
    DEseq_res_ls[[i]]<-DESeq_analysis(cluster_names[i])
    print(i)
  }
  
  res_raw <- bind_rows(DEseq_res_ls)
  colnames(res_raw)[7]=cluster_level
  #if (length(cluster_names) >1){
    #res_raw <- merge(res_raw, cell_annotation,by=cluster_level)
  #}
  gwas_gene<-read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Crispr/data/psy_genes_for_KO/20240617_GWAS_loci_and_candidate_protein-coding_genes.xlsx', 
                        sheet = 'Candidate_gene_human_ortho', range = "A5:AD342",col_names=FALSE)
  res_raw$gwas<-ifelse(paste0("ENSDARG", str_split_i(res_raw$gene, "-ENSDARG",2)) %in% gwas_gene[[1]], 'y', 'n')
  res_flt<-res_raw[abs(res_raw$log2FoldChange)>1,]
  return(list(res_flt, res_raw))
}

####apply the function of DESeq2 to a seurat object####
df_deg_res <- DESeq2_custom(subset(gex.combined.hicat,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07")), 
                      'harmony_clusters_pca_2')

####perform DESeq2 with customized clusters####
##e.g. perform DEG analysis between SDA and VDA for cells from cluster 126 and 127 of the "hicat" class
gex.combined.hicat <- readRDS('./processed_data/gex.combined.hicat.rds')
gex.temp <- subset(gex.combined.hicat, subset=hicat %in% c(126, 127))
gex.temp <- subset(gex.temp,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07"))
gex.temp[['temp_class']] <- 'selected'

subset_sce <- readRDS('./processed_data/subset_sce_001004037.rds')
gex.temp <- subset(subset_sce,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07"))
DESeq2_temp_class <- DESeq2_custom(gex.temp, 'rc_1')

#####################################################################
####################### Part 2. Visualize DEG analysis #######################
#####################################################################
##2.1 general stats of raw results
#make a histogram of log2FoldChange across all differential expression events
df_deg_raw<-do.call(rbind.data.frame, df_deg_res[2])
#write.csv(df_deg_raw, "./processed_data/deg_raw_v1.csv", row.names = FALSE)
#saveRDS(df_deg_raw, './processed_data/df_deg_raw_v1.rds')
df_deg_raw <- readRDS('./processed_data/df_deg_raw_v1.rds')

ggplot(data = df_deg_raw, aes(x = log2FoldChange)) +
  geom_histogram(bins = 20,color = "white", fill='#808080')+theme_classic()+
  geom_segment(aes(x=-1,y=0,xend=-1,yend=400),color='red',linewidth = 0.5)+geom_segment(aes(x=1,y=0,xend=1,yend=400),color='red',linewidth = 0.5)+
  theme(axis.title=element_blank(),axis.text=element_text(size=30))

##2.2 general stats of filtered results
#filtered DEG by excluding mt genes
df_deg<-do.call(rbind.data.frame, df_deg_res[1])
#write.csv(df_deg, "./processed_data/deg_v1.csv", row.names = FALSE)
#saveRDS(df_deg, './processed_data/df_deg_v1.rds')
df_deg <- readRDS('./processed_data/df_deg_v1.rds')
df_deg<-df_deg[df_deg$class!='doublet',]

#make a pie chart of different types of genes with differential expression
gene_cat<-c("mitochondrial", "non-coding", "conserved_protein-coding","non-conserved_protein-coding")
counts<-c(count(grepl('mt-',unique(df_deg$gene))),
          count(grepl('CU|BX|AL|:|CR',unique(df_deg$gene))),
          length(unique(df_deg[str_split_i(df_deg$gene, "-",-1)%in%bg_genes$zebrafish_ens,]$gene)),
          count(!grepl('CU|BX|AL|:|mt-|CR',unique(df_deg$gene)))-length(unique(df_deg[str_split_i(df_deg$gene, "-",-1)%in%bg_genes$zebrafish_ens,]$gene)))
df_plot<-data.frame(cat=gene_cat, pct = counts/sum(counts))
df_plot$labels<-as.factor(scales::percent(df_plot$pct))
ggplot(df_plot, aes(x = "", y = pct, fill = gene_cat))+
  geom_col(color="black")+
  #geom_text(aes(label = labels),color=c("#f2f2f2","#f2f2f2","#000000"), 
  #position = position_stack(vjust = 0.5),
  #show.legend = FALSE) +
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c("#60ddc8", "#9fc5e8",
                               "#EACF65", "#dd6075"))+
  guides(fill = guide_legend(title = ""))+theme(legend.position = "")

###visualized general stats of DEG###
##make a bar chart of number of related clusters per DEG
df_deg<-df_deg[!grepl('mt-', df_deg$gene),]
df_plot_all<-data.frame(rbind(table(table(df_deg$gene)==1), 
                              table((table(df_deg$gene)>1)&(table(df_deg$gene)<5)),
                              table((table(df_deg$gene)>=5)&(table(df_deg$gene)<10)),
                              table(table(df_deg$gene)>=10)))

df_plot_conserved<-df_deg[str_split_i(df_deg$gene,"-",-1)%in%bg_genes$zebrafish_ens,]
df_plot_conserved<-data.frame(rbind(table(table(df_plot_conserved$gene)==1), 
                                    table((table(df_plot_conserved$gene)>1)&(table(df_plot_conserved$gene)<5)),
                                    table((table(df_plot_conserved$gene)>=5)&(table(df_plot_conserved$gene)<10)),
                                    table(table(df_plot_conserved$gene)>=10)))
df_plot_conserved$cluter_no<-c("1", "2-5", "6-10", ">10")

df_plot_nonconserved<-df_deg[str_split_i(df_deg$gene,"-",-1)%notin%bg_genes$zebrafish_ens,]
df_plot_nonconserved<-data.frame(rbind(table(table(df_plot_nonconserved$gene)==1), 
                                       table((table(df_plot_nonconserved$gene)>1)&(table(df_plot_nonconserved$gene)<5)),
                                       table((table(df_plot_nonconserved$gene)>=5)&(table(df_plot_nonconserved$gene)<10)),
                                       table(table(df_plot_nonconserved$gene)>=10)))
df_plot_nonconserved$cluter_no<-c("1", "2-5", "6-10", ">10")

df_plot<-data.frame(cbind(df_plot_conserved[,3],df_plot_conserved[,2],df_plot_nonconserved[,2]))
colnames(df_plot)<-c('category','conserved_deg','nonconserved_deg')
df_plot$conserved_deg<-as.numeric(df_plot$conserved_deg)
df_plot$nonconserved_deg<-as.numeric(df_plot$nonconserved_deg)
fig <- plot_ly(df_plot, x = ~category, y = ~conserved_deg, type = 'bar', name = 'conserved',
               marker=list(color = rep('#AA6926',4))
)
fig <- fig %>% add_trace(y = ~nonconserved_deg, name = 'none-conserved',marker=list(color = rep('#EAD6B8',4)))
fig <- fig %>% layout(yaxis = list(title = '',showgrid=F),
                      xaxis = list(title = '', categoryarray=~category),
                      barmode = 'stack', font=list(size=30), showlegend=F)
fig

##make a bar chart of number of DEGs per cluster
df_plot_conserved<-df_deg[str_split_i(df_deg$gene,"-",-1)%in%bg_genes$zebrafish_ens,]
df_plot_conserved<-data.frame(rbind(table(table(df_plot_conserved$harmony_clusters_pca_2)==1), 
                                    table((table(df_plot_conserved$harmony_clusters_pca_2)>1)&(table(df_plot_conserved$harmony_clusters_pca_2)<5)),
                                    table((table(df_plot_conserved$harmony_clusters_pca_2)>=5)&(table(df_plot_conserved$harmony_clusters_pca_2)<10)),
                                    table(table(df_plot_conserved$harmony_clusters_pca_2)>=10)))
df_plot_conserved$cluter_no<-c("1", "2-5", "6-10", ">10")

df_plot_nonconserved<-df_deg[str_split_i(df_deg$gene,"-",-1)%notin%bg_genes$zebrafish_ens,]
df_plot_nonconserved<-data.frame(rbind(table(table(df_plot_nonconserved$harmony_clusters_pca_2)==1), 
                                       table((table(df_plot_nonconserved$harmony_clusters_pca_2)>1)&(table(df_plot_nonconserved$harmony_clusters_pca_2)<5)),
                                       table((table(df_plot_nonconserved$harmony_clusters_pca_2)>=5)&(table(df_plot_nonconserved$harmony_clusters_pca_2)<10)),
                                       table(table(df_plot_nonconserved$harmony_clusters_pca_2)>=10)))
df_plot_nonconserved$cluter_no<-c("1", "2-5", "6-10", ">10")

df_plot<-data.frame(cbind(df_plot_conserved[,3],df_plot_conserved[,2],df_plot_nonconserved[,2]))
colnames(df_plot)<-c('category','conserved_deg','nonconserved_deg')
df_plot$conserved_deg<-as.numeric(df_plot$conserved_deg)
df_plot$nonconserved_deg<-as.numeric(df_plot$nonconserved_deg)
fig <- plot_ly(df_plot, x = ~category, y = ~conserved_deg, type = 'bar', name = 'conserved',
               marker=list(color = rep('#AA6926',4))
)
fig <- fig %>% add_trace(y = ~nonconserved_deg, name = 'none-conserved',type='bar',marker=list(color = rep('#EAD6B8',4)))
fig <- fig %>% layout(yaxis = list(title = '',showgrid=F,autorange=T),
                      xaxis = list(title = '',categoryarray=~category),
                      font=list(size=30), showlegend=F,barmode = 'stack')
fig

##2.3 visualize expression x cluster matrix
####visualize DEG by cluster###
cell_annotation <- read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/doc/submission_2025/Supplementary_tables/Table S11 Cell annotation for snRNA-seq.xlsx')
cell_annotation <- as.data.frame(cell_annotation)
class_color <- read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/doc/submission_2025/Supplementary_tables/Table S11 Cell annotation for snRNA-seq.xlsx', sheet = 'class_color')
class_color <- as.data.frame(class_color)
de_pval<-DataFrame(matrix(nrow=length(unique(df_deg$gene)),ncol=1))
de_log2fc<-DataFrame(matrix(nrow=length(unique(df_deg$gene)),ncol=1))
rownames(de_pval)<-unique(df_deg$gene)
rownames(de_log2fc)<-unique(df_deg$gene)
for (cluster_id in unique(df_deg$harmony_clusters_pca_2)){
  temp<-DataFrame(matrix(nrow=length(unique(df_deg$gene)),ncol=1))
  rownames(temp)<-unique(df_deg$gene)
  shared_gene<-df_deg[df_deg$harmony_clusters_pca_2==cluster_id,]$gene
  
  temp[shared_gene,1]<-df_deg[df_deg$harmony_clusters_pca_2==cluster_id&df_deg$gene %in% shared_gene,]$padj
  colnames(temp)[1]<-as.character(cluster_id)
  de_pval<-cbind(de_pval,temp)
  
  temp[shared_gene,1]<-df_deg[df_deg$harmony_clusters_pca_2==cluster_id&df_deg$gene %in% shared_gene,]$log2FoldChange
  colnames(temp)[1]<-as.character(cluster_id)
  de_log2fc<-cbind(de_log2fc,temp)
}
de_pval<-de_pval[,-1]
de_pval<-as.matrix(de_pval)
de_pval[is.na(de_pval)]<-1

de_log2fc<-de_log2fc[,-1]
de_log2fc<-as.matrix(de_log2fc)
de_log2fc[is.na(de_log2fc)]<-0
#find genes upper regulated in VDA in all clusters and sort them based on number of related clusters
gene_pos_only<-rownames(de_log2fc)[rowSums(de_log2fc>=0)==dim(de_log2fc)[[2]]] 
gene_pos_only<-gene_pos_only[order(rowSums(de_log2fc[gene_pos_only,]>0),decreasing=T)]
#find genes down regulated in VDA in all clusters and sort them based on number of related clusters
gene_neg_only<-rownames(de_log2fc)[rowSums(de_log2fc<=0)==dim(de_log2fc)[[2]]]
gene_neg_only<-gene_neg_only[order(rowSums(de_log2fc[gene_neg_only,]<0),decreasing=T)]
#find genes regulated in both directions
gene_neg_pos<-setdiff(rownames(de_log2fc), c(gene_pos_only, gene_neg_only))
#order clusters based on its class 
ordered_clusters<-as.character(cell_annotation[order(cell_annotation$class, decreasing = F)&cell_annotation$harmony_clusters_pca_2 %in% colnames(de_pval), ]$harmony_clusters_pca_2)
col_colors<-cell_annotation$color
names(col_colors) <- cell_annotation$class
#col_colors<-merge(cell_annotation[cell_annotation$harmony_clusters_pca_2 %in% ordered_clusters,], class_color,by='class')$color
#names(col_colors)<-merge(cell_annotation[cell_annotation$harmony_clusters_pca_2 %in% ordered_clusters,], class_color,by='class')$class
col_colors<-as.matrix(col_colors)
colnames(col_colors)<-'cell_type'

####heatmap with DEG sorted by cluster specificity and conservation
gene_list_pos<-c()
gene_list_neg<-c()
#gene_list_pos_sc<-c()
#gene_list_pos_mc<-c()
for (cluster in unique(sort(df_deg$cluster))){
  #get positive deg for one cluster
  gene_list_temp_pos<-df_deg[df_deg$cluster==cluster&df_deg$log2FoldChange>0,]$gene
  gene_list_temp_pos<-table(df_deg[df_deg$gene%in%gene_list_temp_pos&df_deg$log2FoldChange>0,]$gene)
  gene_list_temp_neg<-df_deg[df_deg$cluster==cluster&df_deg$log2FoldChange<0,]$gene
  gene_list_temp_neg<-table(df_deg[df_deg$gene%in%gene_list_temp_neg&df_deg$log2FoldChange<0,]$gene)
  #sort deg based on number of related clusters ascendingly
  gene_list_pos<-c(gene_list_pos,names(sort(gene_list_temp_pos)))
  gene_list_neg<-c(gene_list_neg,names(sort(gene_list_temp_neg)))
  #gene_list_temp_pos_sc<-names(gene_list_temp_pos)[gene_list_temp_pos==1]
  #gene_list_temp_pos_mc<-names(gene_list_temp_pos)[gene_list_temp_pos>1]
  #gene_list_temp_pos_mc<-setdiff(gene_list_temp_pos_mc,gene_list_temp_pos_sc)
  #gene_list_pos_sc<-c(gene_list_pos_sc,gene_list_temp_pos_sc)
  #gene_list_pos_mc<-c(gene_list_pos_mc,gene_list_temp_pos_mc)
}

df_deg_spe<-data.frame(table(df_deg$class))
for(gene in unique(df_deg$gene)){
  df_deg_spe<-merge(df_deg_spe,data.frame(table(df_deg[df_deg$gene==gene,]$class)),by="Var1", all.x=T)
  df_deg_spe[is.na(df_deg_spe)]<-0
}
df_deg_spe<-sweep(df_deg_spe[,-1], MARGIN = 2, STATS = colSums(df_deg_spe[,-1]), FUN = "/")
df_deg_spe<-df_deg_spe[,-1]
colnames(df_deg_spe)<-unique(df_deg$gene)
rownames(df_deg_spe)<-data.frame(table(df_deg$class))$Var1
df_deg_spe<-df_deg_spe[,1:(ncol(df_deg_spe)-1)]
df_deg_spe_cat<-data.frame(class=apply(df_deg_spe, MARGIN=2, which.max),
                           pct_max=sapply(df_deg_spe, max, na.rm = TRUE))
df_deg_spe_cat$class<-rownames(df_deg_spe)[df_deg_spe_cat$class]
df_deg_spe_cat$conserved<-'n'
df_deg_spe_cat$dir<-'n'
df_deg_spe_cat$ne_vs_nn<-'n'
df_deg_spe_cat$region_spe<-2

df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe[,df_deg_spe["15_NN",]==0]),]$ne_vs_nn<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe[,df_deg_spe["15_NN",]>0&df_deg_spe["15_NN",]<=0.5]),]$ne_vs_nn<-2
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe[,df_deg_spe["15_NN",]==1]),]$ne_vs_nn<-3

df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_pos)[str_split_i(unique(gene_list_pos),"-",-1)%in% bg_genes$zebrafish_ens],"conserved"]<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_pos)[str_split_i(unique(gene_list_pos),"-",-1)%in% bg_genes$zebrafish_ens],"dir"]<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_pos)[str_split_i(unique(gene_list_pos),"-",-1)%notin% bg_genes$zebrafish_ens],"conserved"]<-2
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_pos)[str_split_i(unique(gene_list_pos),"-",-1)%notin% bg_genes$zebrafish_ens],"dir"]<-1

df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_neg)[str_split_i(unique(gene_list_neg),"-",-1)%in% bg_genes$zebrafish_ens],"conserved"]<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_neg)[str_split_i(unique(gene_list_neg),"-",-1)%in% bg_genes$zebrafish_ens],"dir"]<-2
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_neg)[str_split_i(unique(gene_list_neg),"-",-1)%notin% bg_genes$zebrafish_ens],"conserved"]<-2
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% unique(gene_list_neg)[str_split_i(unique(gene_list_neg),"-",-1)%notin% bg_genes$zebrafish_ens],"dir"]<-2

df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe)[(colSums(df_deg_spe[grepl("_FB",rownames(df_deg_spe)),])>0.8)],]$region_spe<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe)[(colSums(df_deg_spe[grepl("_MB",rownames(df_deg_spe)),])>0.8)],]$region_spe<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe)[(colSums(df_deg_spe[grepl("_HB",rownames(df_deg_spe)),])>0.8)],]$region_spe<-1
df_deg_spe_cat[rownames(df_deg_spe_cat) %in% colnames(df_deg_spe)[(colSums(df_deg_spe[grepl("_OB|_OE",rownames(df_deg_spe)),])>0.8)],]$region_spe<-1

gene_order<-rownames(df_deg_spe_cat[order(df_deg_spe_cat$conserved,
                                          df_deg_spe_cat$dir,
                                          df_deg_spe_cat$ne_vs_nn,
                                          df_deg_spe_cat$region_spe,
                                          df_deg_spe_cat$class,
                                          df_deg_spe_cat$pct_max),])
gene_order<-gene_order[gene_order %in% rownames(df_deg_spe_cat[df_deg_spe_cat$conserved==1,])]
#gene_order<-gene_order[gene_order %in% rownames(df_deg_spe_cat[df_deg_spe_cat$conserved==2,])] #run this to get results of non-conserved genes
library(heatmaply)
p<-heatmaply(de_log2fc[gene_order,ordered_clusters],
             width=1212, height=578,
             dendrogram = 'none',show_dendrogram = c(FALSE, TRUE),
             showticklabels=FALSE,grid_gap = 0.5,
             colors = c(as.vector(paletteer_c("ggthemes::Red-Blue Diverging", 30))[1:3],'#EBE0D9', as.vector(paletteer_c("ggthemes::Red-Blue Diverging", 30))[28:30]),
             limits = c(-8,8),
             margins = c(0,0,60,0),colorbar_xanchor='middle',colorbar_yanchor='top',
             colorbar_xpos=1.0, colorbar_ypos=1.2,plot_method="plotly",colorbar_thickness=20,colorbar_len=0.2
             #ColSideColors = col_colors, 
             #RowSideColors = c(rep('#60ddc8',163),rep('#dd6075',91)),
             #subplot_heights=c(0.02,0.98),subplot_widths=c(0.98,0.02)
)
p$x$layout$showlegend<-FALSE
p$x$layout$xaxis2$ticktext<-""
p$x$layout$yaxis$ticktext<-""
p

##a stacked bar chart to show cell type specificity of each DEG
df_plot<-data.frame(t(df_deg_spe))
df_plot<-df_plot[gene_order,]
df_plot[,'gene']<-rownames(df_plot)
df_plot<-reshape2::melt(df_plot, id.vars='gene')
df_plot$class<-str_split_i(df_plot$variable,"X",-1)
col_colors_temp<-data.frame(col_colors)
col_colors_temp$class<-str_split_i(rownames(col_colors_temp),"X",-1)
df_plot<-merge(df_plot, col_colors_temp, by='class')
df_plot$cell_type<-as.character(df_plot$cell_type)
plot_ly(df_plot, x = ~gene, y = ~value, type = 'bar', width=880,height=80,
        marker=list(color = df_plot$cell_type)) %>% 
  layout(margin = list(l = 0,r = 0,b = 0,t = 0),
         yaxis = list(title = '',showticklabels=FALSE), xaxis = list(categoryarray=~gene,title='',showticklabels=FALSE),
         barmode = 'stack')

##table of up and down regulated DEG number per class
df_temp1 <- df_deg[df_deg$gene %in% rownames(df_deg_spe_cat[df_deg_spe_cat$conserved==1,]),]
df_temp1 <- df_deg[df_deg$log2FoldChange>0, ]
df_temp1 <- table(df_temp1[!duplicated(df_temp1[c("gene", "class")]), ]$class)
df_temp1['08_HB-GABA']<-0
df_temp2 <- df_deg[df_deg$gene %in% rownames(df_deg_spe_cat[df_deg_spe_cat$conserved==1,]),]
df_temp2 <- df_deg[df_deg$log2FoldChange<0, ]
df_temp2 <- table(df_temp2[!duplicated(df_temp2[c("gene", "class")]), ]$class)
df_temp2['08_HB-GABA']<-0
df_temp1[setdiff(names(df_temp2), names(df_temp1))] <- 0
df_temp2[setdiff(names(df_temp1), names(df_temp2))] <- 0
df_temp <- rbind(df_temp1[sort(names(df_temp1))], df_temp2[sort(names(df_temp2))])
df_cellnote <- df_temp
df_temp[2,] <- df_temp[2,]*(-1)
p<-heatmaply(
  t(df_temp),
  width=200, height=500,
  dendrogram = 'none',show_dendrogram = c(FALSE, TRUE),
  cellnote = t(df_cellnote), draw_cellnote = TRUE,cellnote_textposition = 'middle center',cellnote_color='black',cellnote_size = 15,
  showticklabels=FALSE,grid_gap = 0.1,grid_color = 'black',
  colors = c(as.vector(paletteer_c("ggthemes::Red-Blue Diverging", 30))[1:5],'#ffffff',as.vector(paletteer_c("ggthemes::Red-Blue Diverging", 30))[26:30]),
  limit= c(-50, 50), hide_colorbar = TRUE
  )
p$x$layout$showlegend<-FALSE
p$x$layout$xaxis2$ticktext<-""
p$x$layout$yaxis$ticktext<-""
p$x$layout$xaxis$showline<-FALSE
p$x$layout$yaxis$showline<-FALSE
p
##dotplot of nuclei numbers for each cluster
df_plot<-table(subset(gex.combined.hicat,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07"))$harmony_clusters_pca_2)[ordered_clusters]
df_plot<-data.frame(df_plot)
colnames(df_plot)<-c("cluster", "size")
df_plot$color<-sapply(df_plot$cluster, function(x)cell_annotation$color[which(cell_annotation$harmony_clusters_pca_2==x)])
df_plot$xpos<-as.numeric(rownames(df_plot))
df_plot$ypos<-5
colnames(df_plot)[3]<-'color'
g<-plot_ly(df_plot,x=~xpos,y=~ypos,
           width=1300,height=97,type = 'scatter',mode='markers',symbol = 0,symbols=c('circle-dot'),
           size = ~size, marker=list(color = as.character(df_plot$color),
                                     line = list(color = df_plot$color,width = 0.5)
           )
)%>% layout(margin = list(l = 0,r = 0,b = 0,t = 0),
            yaxis = list(title = '',showticklabels=FALSE,showgrid=F,showticklabels=FALSE), 
            xaxis = list(categoryarray=~gene,title='',showgrid=F,showticklabels=FALSE))
g

##count number of clusters per deg###
df_deg_count<-de_log2fc[gene_order,ordered_clusters]
df_deg_count<-data.frame(cluster_no=rowSums(df_deg_count!=0), gene=rownames(df_deg_count))
p_bar<-plot_ly(df_deg_count,width=160,height=800,
               x=~cluster_no,y=~gene,type="bar",
               marker=list(color = c(rep('#6c9a8b',163)))
)
p_bar <- p_bar %>% layout(margin = list(l = 0,r = 10,b = 10,t = 10),
                          xaxis = list(showticklabels=F,showgrid=T,gridwidth=2,autorange='reversed',title = "",tickfont=list(size=25)),
                          yaxis = list(showticklabels=F,showgrid=F, categoryarray=~gene,autorange='reversed',title = ""))
p_bar

##count number of deg per cluster###
df_deg_count<-de_log2fc[gene_order,ordered_clusters]
df_cluster_count<-data.frame(gene_no=colSums(df_deg_count!=0), cluster_id=colnames(df_deg_count))
p_bar<-plot_ly(df_cluster_count,width=1060,height=74,
               y=~gene_no,x=~cluster_id,type="bar",
               marker=list(color = col_colors[,1])
)
p_bar <- p_bar %>% layout(margin = list(l = 0,r = 5,b = 5,t = 5),
                          xaxis = list(showticklabels=FALSE,showgrid=F, categoryarray=~cluster_id,title=F),
                          yaxis = list(showticklabels=T,showgrid=T,title=F,gridwidth=2,title = "",tickfont=list(size=10)))
p_bar
##label genes in enriched GO term
target_term<-'oxidative phosphorylation'
hgenes<-as.list(strsplit(goea_de_hgene_all_pos[goea_de_hgene_all_pos$term_name==target_term, ]$intersection,",")[[1]])
hgenes<-c('ENSG00000110492', 'ENSG00000125534', 'ENSG00000116329', 'ENSG00000033170',
          'ENSG00000168028', 'ENSG00000140612', 'ENSG00000150672', 'ENSG00000151150',
          'ENSG00000129158', 'ENSG00000096433', 'ENSG00000141668')
target_genes<-bg_genes[bg_genes$human_ens%in%hgenes & bg_genes$zebrafish_ens%in%str_split_i(df_deg$gene,"-",-1), ]
de_target_genes<-df_deg[str_split_i(df_deg$gene,"-",-1) %in% target_genes$zebrafish_ens,]

df_plot<-data.frame(t(df_deg_spe))[gene_order,]
df_plot<-rownames_to_column(df_plot, 'gene')
y_cor<-as.numeric(rownames(df_plot[df_plot$gene %in% unique(de_target_genes$gene),]))
y_cor<-nrow(df_plot)-y_cor+1

g<-plot_ly(data.frame(x=rep(1,length(y_cor)), y=y_cor),x=~x,y=~y,
           width=150,height=750,type = 'scatter',mode='markers',symbol = ~x,symbols = c("star-dot"),
           marker = list(size = 5,color = '#16A697',
                         line = list(color = '#16A697',width = 0.5)
           )
) 
g<-g%>%layout(
  margin = list(l = 0,r = 0,b = 0,t = 0),
  xaxis = list(showticklabels=F,showgrid=F, title="",range=c(0,2)),
  yaxis = list(showticklabels=F,showgrid=F,title="", range=c(0, nrow(df_plot)))
)
g
