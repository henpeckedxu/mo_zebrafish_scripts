library(gprofiler2)
library(readr)
library(stringr)
library(heatmaply)

df_deg <- readRDS('./processed_data/df_deg_v1.rds')
df_temp <- bg_genes[!is.na(bg_genes$human_symbol),]

df_query <- df_deg
df_query <- df_query[df_query$gene %in% intersect(df_query$gene, paste0(df_temp$zebrafish_symbol, '-', df_temp$zebrafish_ens)), ]
df_query <- df_query[df_query$log2FoldChange<0,]
goea_res <- goea_deg(df_query, 'all')[[2]]
saveRDS(goea_res, './processed_data/GOEA/goea_de_hgene_all_neg.rds')
write_csv(goea_res, './processed_data/GOEA/goea_de_hgene_all_neg.csv')

df_query <- df_deg
df_query <- df_query[df_query$gene %in% intersect(df_query$gene, paste0(df_temp$zebrafish_symbol, '-', df_temp$zebrafish_ens)), ]
df_query <- df_query[df_query$log2FoldChange>0,]
goea_res <- goea_deg(df_query, 'all')[[2]]
saveRDS(goea_res, './processed_data/GOEA/goea_de_hgene_all_pos.rds')
write_csv(goea_res, './processed_data/GOEA/goea_de_hgene_all_pos.csv')

df_query <- df_deg
df_query <- df_query[df_query$gene %in% intersect(df_query$gene, paste0(df_temp$zebrafish_symbol, '-', df_temp$zebrafish_ens)), ]
goea_res <- goea_deg(df_query, 'all')[[2]]
saveRDS(goea_res, './processed_data/GOEA/goea_de_hgene_all.rds')
write_csv(goea_res, './processed_data/GOEA/goea_de_hgene_all.csv')



goea_de_hgene_all_pos<- readRDS('./processed_data/GOEA/goea_de_hgene_all_pos.rds')
goea_de_hgene_all_pos <- goea_res[goea_res$term_id %in% names(sort(table(goea_res$term_id))),]


#bar plot of enriched GO terms
target_terms<-c('GO:0005184','GO:0007218','GO:0098992','HPA:0260141','KEGG:04080')
df_plot<-goea_de_hgene_all_pos[goea_de_hgene_all_pos$term_id %in% target_terms,]
df_plot$log10_pvalue<-log(df_plot$p_value,base=10)*(-1)
df_plot$color<-c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600")
fig <- plot_ly(width=600, height = 600, 
               df_plot, x = ~term_name, y = ~log10_pvalue, type = 'bar', name = '',
               marker=list(color = df_plot$color))%>%
  layout(margin = list(l = 0,r = 0,b = 0,t = 0),
         yaxis = list(title = '',showticks=T,showticklabels=F,showgrid=F), xaxis = list(categoryarray=~term_name,title='',showticklabels=F)
  )
fig
##find zebrafish genes enriched for one specific term
target_terms<-c('GO:0005184','GO:0007218','GO:0098992','HPA:0260141','KEGG:04080')
hgenes<-unique(unlist(strsplit(goea_de_hgene_all_pos[goea_de_hgene_all_pos$term_id %in% target_terms,]$intersection,","), recursive = FALSE))
target_genes<-bg_genes[bg_genes$human_ens%in%hgenes & bg_genes$zebrafish_ens%in%str_split_i(df_deg$gene,"-",-1), ]
de_target_genes<-df_deg[str_split_i(df_deg$gene,"-",-1) %in% target_genes$zebrafish_ens,]
de_target_genes<-de_target_genes[grepl('-',de_target_genes$gene),]

##find overlap between the gene set of the enriched terms and gene sets for human psychiatric diseases
scz_genes_overlap <- c()
scz_genes_all <- c()
for (cell_type in excel_sheets("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/SCZ_scRNA_PMID38781388_s4.xlsx")){
  scz_genes <- read_excel("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/SCZ_scRNA_PMID38781388_s4.xlsx", sheet =cell_type)
  scz_genes <- as.data.frame(scz_genes[(scz_genes$Meta_adj.P.Val<0.05)&(abs(scz_genes$Meta_logFC)>0.1), ])
  scz_genes$cell_type = cell_type
  scz_genes_all <- rbind(scz_genes_all, scz_genes)
  scz_genes <- scz_genes[scz_genes$gene %in% bg_genes[bg_genes$human_ens%in%hgenes, ]$human_symbol, ]
  if (dim(scz_genes)[[1]]>0){
    scz_genes_overlap <- rbind(scz_genes_overlap, scz_genes)
  }
}

##find overlap between the gene set of the enriched terms and gene sets for human psychiatric diseases
asd_genes1 <- read_excel("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/ASD_scRNA_PMID38781372_s4.xlsx",sheet = 1,skip = 2)
asd_genes <- unique(unlist(strsplit(asd_genes1$`Genes (by rows for category)`, "\r\n")))
asd_genes <- unique(unlist(strsplit(asd_genes, " ")))
asd_genes_overlap <- intersect(asd_genes, bg_genes[bg_genes$human_ens%in%hgenes, ]$human_symbol)
##find overlap between the gene set of the enriched terms and gene sets for human PTSD and MDD
ptsd_genes_overlap <- c()
mdd_genes_overlap <- c()
mast_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_MAST.txt", sep='\t')
wilcox_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_Wilcox.txt", sep='\t')
celltypelist <- c('IN', 'LAMP5','KCNG1', 'VIP', 'SST','PVALB', 'EXN', 'CUX2','RORB','FEZF2','OPRK1')
ptsd_mdd_deg_all <- cbind(mast_res ,wilcox_res[,5:8])
ptsd_mdd_deg_all <- ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist, ]
ptsd_deg_all <- ptsd_mdd_deg_all[(ptsd_mdd_deg_all$PTSD.MAST.FDR<0.05)&(ptsd_mdd_deg_all$PTSD.Wilcox.FDR<0.05),]
ptsd_deg_all <- ptsd_deg_all[mapply(function(a,b)max(abs(a), abs(b))>log2(1.2), ptsd_deg_all$PTSD.MAST.log2FC, ptsd_deg_all$PTSD.Wilcox.log2FC),]
ptsd_deg_all <- ptsd_deg_all[!is.na(ptsd_deg_all),]
ptsd_deg_all <- ptsd_deg_all[, !grepl('MDD.',colnames(ptsd_deg_all))]



ptsd_mast_res <- mast_res[(mast_res$PTSD.MAST.FDR<0.05)&(mast_res$Celltype %in% celltypelist),]
mdd_mast_res <- mast_res[(mast_res$MDD.MAST.FDR<0.05)&(mast_res$Celltype %in% celltypelist),]
ptsd_wilcox_res <- wilcox_res[(wilcox_res$PTSD.Wilcox.FDR<0.05)&(wilcox_res$Celltype %in% celltypelist),]
mdd_wilcox_res <- wilcox_res[(wilcox_res$MDD.Wilcox.FDR<0.05)&(wilcox_res$Celltype %in% celltypelist),]

ptsd_genes <- intersect(ptsd_mast_res$Genename, ptsd_wilcox_res$Genename)
mdd_genes <- intersect(mdd_mast_res$Genename, mdd_wilcox_res$Genename)

ptsd_genes_overlap <- intersect(ptsd_genes, bg_genes[bg_genes$human_ens%in%hgenes, ]$human_symbol)
mdd_genes_overlap <- intersect(mdd_genes, bg_genes[bg_genes$human_ens%in%hgenes, ]$human_symbol)

library(ggvenn)
ggvenn(
  list(
    'ASD' = setdiff(asd_genes_overlap, c('GRP','GHR')),
    'SCZ'= setdiff(unique(scz_genes_overlap$gene), c('GRP','GHR')), 
    'PTSD'= setdiff(ptsd_genes_overlap,c('GRP','GHR')),
    'MDD'=setdiff(mdd_genes_overlap,c('GRP','GHR'))),
  fill_color = c('#a559aa',"#59a89c", "#f0c571", "#e02b35"),
  stroke_size = 1, set_name_size = 10,text_size = 8,
  show_percentage =FALSE
)

library(ggvenn)
ggvenn(
  list(
    'DEG' = unique(bg_genes[bg_genes$human_ens%in%hgenes, ]$human_symbol),
    'PTSD'= ptsd_genes), 
  fill_color = c('#a559aa',"#59a89c", "#f0c571", "#e02b35"),
  stroke_size = 1, set_name_size = 10,text_size = 0,
  show_set_totals = "none",
  show_stats = 'none',
  show_percentage =FALSE,
  auto_scale = TRUE
)
##heatmap of DEGs related to top enriched terms
npp_cluster_order <- unique(de_target_genes[order(de_target_genes$cluster),]$harmony_clusters_pca_2)
npp_gene_order <- unique(de_target_genes[order(de_target_genes$cluster),]$gene)

p<-heatmaply(de_log2fc[npp_gene_order,npp_cluster_order],
             width=300, height=500,
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

##Venn Diagram of DEnpgs expressed in clusters of different neurotransmitter types
library(eulerr)
list_data <- list(
  GABA = unique(de_target_genes$gene[grepl('GABA',de_target_genes$class)]),
  Glu = unique(de_target_genes$gene[grepl('Glu',de_target_genes$class)])
)
fit <- euler(list_data)
plot(fit, fills = c("#fbd501", "#CC79A7"), labels = '', quantities=FALSE)



