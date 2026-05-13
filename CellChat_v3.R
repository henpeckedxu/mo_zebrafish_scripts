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
library(CellChat)
library(Seurat)
library(RColorBrewer)
####This script is to perform cell-cell communcation analysis with zebrafish genes
##############Chapter 1. develop a database for cellchat####
#<notes>:the original database is based on human genes. we need change genes in all layers of the database to zebrafish genes
#Procedures
##1.1. read raw data
##cluster annotation 
cluster_annotation <- read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/doc/submission_2025/Supplementary_tables/Table S11 Cell annotation for snRNA-seq.xlsx')
##human database for cellchat
hcellchat_db <- CellChatDB.human
##human-zebrafish orthologus
bg_genes <- read_excel('/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/ZebrafishGenomeAnnotation/ortho_genes_20250227.xlsx', sheet = 1)
bg_genes <- bg_genes[bg_genes$human_ens!="#N/A",]
##DEnpg data
de_target_genes <- readRDS('./processed_data/GOEA/deg_goea.rds')
##snRNA-seq object
gex.combined.hicat <- readRDS('./processed_data/snRNA_Analysis/gex.combined.hicat.rds')

##1.2. get the four datasets of the database
h_interaction <-hcellchat_db$interaction
h_complex <- hcellchat_db$complex
h_cofactor <- hcellchat_db$cofactor
h_geneinfo <- hcellchat_db$geneInfo

##1.3. convert human gene to zebrafish gene in interaction dataset
###1.3.1. develop a function to convert human symbol to zebrafish symbol
#<notes>:complexes will be converted with all possible combinations of zebrafish orthologs.
h2z <- function(h_symbol){
  if(length(str_split(h_symbol,', ')[[1]])==1){
    #if there is only one gene listed as the ligand or receptor, simply return zebrafish ortholog(s)
    #if there are multiple zebrafish orthologs, return them all separated by ';' 
    #e.g: h2z('CRH') --> 'crhb;crha'
    z_symbol = bg_genes$zebrafish_symbol[which(bg_genes$human_symbol==h_symbol)]
    z_symbol = paste(z_symbol, collapse=';')
  }
  
  else{
    #if there are multiple genes listed as the ligand or receptor
    #return all combos of the zebrafish orthologs seperated by '|'
    #e.g.h2z('CRH, MEF2C') --> 'crhb,mef2cb|crha,mef2cb|crhb,mef2ca|crha,mef2ca'
    z_symbol = str_split(h_symbol,', ')[[1]]
    temp_list = c()
    for(gene in z_symbol){
      temp_element <- bg_genes$zebrafish_symbol[which(bg_genes$human_symbol==gene)]
      temp_element <- list(temp_element)
      temp_list = c(temp_list, temp_element)
    }
    z_symbol = paste(do.call(paste, c(expand.grid(temp_list), sep = ",")), collapse = '|')
  }
  return (z_symbol)
}
###1.3.2. apply the symbol conversion function to two columns of the interaction dataset: 'ligand.symbol' and 'receptor.symbol
#copy human interaction databset to a new variable name
z_interaction <- h_interaction
#apply the function and create two new column after conversion: zreceptor.symbol and zligand.symbol
z_interaction$zreceptor.symbol<-sapply(z_interaction$receptor.symbol, function(x)h2z(x))
z_interaction$zligand.symbol<-sapply(z_interaction$ligand.symbol, function(x)h2z(x))
#separate one row into multiple rows if there are multiple genes returned with the function
z_interaction <- separate_rows(z_interaction, zreceptor.symbol, sep = '\\|')
z_interaction <- separate_rows(z_interaction, zreceptor.symbol, sep = '\\, ')
z_interaction <- separate_rows(z_interaction, zreceptor.symbol, sep = ';')
z_interaction <- separate_rows(z_interaction, zligand.symbol, sep = '\\|')
z_interaction <- separate_rows(z_interaction, zligand.symbol, sep = ';')
#remove rows with empty ligand or receptors after conversion
z_interaction <- z_interaction[!is.na(z_interaction$zligand.symbol)&!is.na(z_interaction$zreceptor.symbol),]
#create two new columns: zligand and zreceptor by connecting multiple genes within each cell with '_' 
z_interaction$zligand <- sapply(z_interaction$zligand.symbol, function(x)paste(str_split(x,',')[[1]],collapse = '_'))
z_interaction$zreceptor <- sapply(z_interaction$zreceptor.symbol, function(x)paste(str_split(x,',')[[1]],collapse = '_'))

##1.4. convert human gene to zebrafish gene in complex dataset
###1.4.1 ligand complex dataset
z_complex_ligand <- z_interaction[z_interaction$ligand %in% rownames(h_complex),c('ligand','zligand.symbol','zligand')]
z_complex_ligand <- z_complex_ligand[z_complex_ligand$zligand.symbol!="",]
z_complex_ligand <- z_complex_ligand[!duplicated(z_complex_ligand$zligand),]
z_complex_ligand <- column_to_rownames(z_complex_ligand,var='zligand')
z_complex_ligand <- separate(z_complex_ligand,zligand.symbol, into=colnames(h_complex), sep = ",")
z_complex_ligand$ligand <- NULL

###1.4.1 receptor complex dataset
z_complex_receptor <- z_interaction[z_interaction$receptor %in% rownames(h_complex),c('receptor','zreceptor.symbol','zreceptor')]
z_complex_receptor <- z_complex_receptor[z_complex_receptor$zreceptor.symbol!="",]
z_complex_receptor <- z_complex_receptor[!duplicated(z_complex_receptor$zreceptor),]
z_complex_receptor <- column_to_rownames(z_complex_receptor,var='zreceptor')
z_complex_receptor <- separate(z_complex_receptor,zreceptor.symbol, into=colnames(h_complex), sep = ",")
z_complex_receptor$receptor <- NULL

###1.4.2 combine the two complex dataset
z_complex <- rbind(z_complex_ligand,z_complex_receptor)

##1.5. convert human gene to zebrafish gene in cofactor dataset
z_cofactor <- as.data.frame(apply(h_cofactor,c(1,2), function(x)h2z(x)))
z_cofactor$h_cofactor <- rownames(z_cofactor)
for(i in c(1:(ncol(z_cofactor)-1))){
  z_cofactor<-separate_rows(z_cofactor, i, sep = ";")
}
z_cofactor$z_cofactor_temp <- z_cofactor$h_cofactor
##zebrafish has multiple paralogs sharing orthology with a same human gene
##therefore, one human cofactor group corresponds to mutiple zebrafish cofactor groups
##we split one human cofactor into multiple rows in this cofactor dataset and assign a unique name using the z_cofator column
dup_id <- ave(z_cofactor$z_cofactor_temp, z_cofactor$z_cofactor_temp, FUN = seq_along)
z_cofactor$z_cofactor <- paste0(z_cofactor$z_cofactor_temp, "_", dup_id)

##1.6 update columns in interaction dataset related to the cofactor dataset
#<notes>: there are four columns related: 'agonist','antagonist','co_I_recepotr' and 'co_A_receptor'
#<notes>: the h_cofactor column in the z_cofactor dataset is used as the key for update
###1.6.1 update
z_interaction$z_agonist <- lapply(z_interaction$agonist, function(x)paste(z_cofactor$z_cofactor[which(z_cofactor$h_cofactor==x)],collapse = ';'))
z_interaction <- separate_rows(z_interaction,z_agonist, sep = ";")
z_interaction$z_antagonist <- lapply(z_interaction$antagonist, function(x)paste(z_cofactor$z_cofactor[which(z_cofactor$h_cofactor==x)],collapse = ';'))
z_interaction <- separate_rows(z_interaction,z_antagonist, sep = ";")
z_interaction$z_co_A_receptor <- lapply(z_interaction$co_A_receptor, function(x)paste(z_cofactor$z_cofactor[which(z_cofactor$h_cofactor==x)],collapse = ';'))
z_interaction <- separate_rows(z_interaction,z_co_A_receptor, sep = ";")
z_interaction$z_co_I_receptor <- lapply(z_interaction$co_I_receptor, function(x)paste(z_cofactor$z_cofactor[which(z_cofactor$h_cofactor==x)],collapse = ';'))
z_interaction <- separate_rows(z_interaction,z_co_I_receptor, sep = ";")
###1.6.2 convert z_cofactor column as row name for z_cofactor dataset
z_cofactor <- column_to_rownames(z_cofactor,var='z_cofactor')
###1.6.3 remove last two columns in the z_cofactor dataset
z_cofactor <- z_cofactor[, c(1:(ncol(z_cofactor)-2))]

##1.7 convert human gene to zebrafish gene in gene_info dataset
z_gene_info <- bg_genes[bg_genes$human_symbol %in% h_geneinfo$Symbol,c(2:6)]
colnames(z_gene_info) <- plyr::mapvalues(colnames(z_gene_info), from = c("zebrafish_symbol"), to = c("Symbol"),warn_missing = TRUE)

##1.8 rename z_interaction column
#<notes>: columns with the human version will be removed
#<notes>: the remaining columns will be renamed and organized in the same way as the origin dataset
#<notes>: row name and interaction name will made based on original column but will add appendix to make it unique 
dup_id <- ave(z_interaction$interaction_name, z_interaction$interaction_name, FUN = seq_along)
z_interaction$interaction_name <- paste0(z_interaction$interaction_name , "_", dup_id)
z_interaction <- column_to_rownames(z_interaction,var='interaction_name')
z_interaction$interaction_name <- rownames(z_interaction)
dup_id <- ave(z_interaction$interaction_name_2, z_interaction$interaction_name_2, FUN = seq_along)
z_interaction$interaction_name_2 <- paste0(z_interaction$interaction_name_2 , "_", dup_id)
z_interaction_new<-z_interaction[,c('interaction_name', 'pathway_name','zligand', 'zreceptor',
                                    'z_agonist','z_antagonist','z_co_A_receptor', 'z_co_I_receptor', 
                                    'annotation','interaction_name_2', 'evidence','is_neurotransmitter',
                                    'zligand.symbol','ligand.family','ligand.location','ligand.keyword','ligand.secreted_type','ligand.transmembrane',
                                    'zreceptor.symbol', 'receptor.family', 'receptor.location', 'receptor.keyword', 'receptor.surfaceome_main','receptor.surfaceome_sub',
                                    'receptor.adhesome', 'receptor.secreted_type', 'receptor.transmembrane', 'version')]
colnames(z_interaction_new)<-colnames(h_interaction)
zcellchat_db <- list()
zcellchat_db$interaction <- z_interaction_new
zcellchat_db$complex <- z_complex
zcellchat_db$cofactor<-z_cofactor
zcellchat_db$geneInfo<-z_gene_info

##1.9 manually add or modify some LR pairs in cellchatDB for specific analysis purpose
##1) add npy8ar as the receptor for npy
df_temp <- zcellchat_db$interaction['NPY_NPY1R_1',]
rownames(df_temp)<-"NPY_NPY8AR_1"
df_temp[] <- lapply(df_temp, function(x) {
  if (is.character(x)) gsub("NPY1R", "npy8ar", x,ignore.case = TRUE)
  else x
})
zcellchat_db$interaction <- rbind(zcellchat_db$interaction, df_temp)
cols <- colnames(zcellchat_db$geneInfo)
values <- list("npy8ar", 'ENSDARG00000017234', 'neuropeptide Y receptor Y8a', '', '')
df_temp <- setNames(as.data.frame(values), cols)
zcellchat_db$geneInfo <- rbind(zcellchat_db$geneInfo, df_temp)

##2) add npbwr2a as receptor for npb
df_temp <- zcellchat_db$interaction['NPB_NPBWR1_1',]
df_temp$receptor <- 'npbwr2a'
df_temp$receptor.symbol <- 'npbwr2a'
zcellchat_db$interaction['NPB_NPBWR1_1',] <- df_temp
saveRDS(zcellchat_db, './processed_data/CellChat/zcellchat_db_v1.rds')

##############Chapter 2. CellChat analysis for SDA vs. VDA##############
##2.1 select LR pairs only related to DEnpgs
denpgs <- unique(de_target_genes$gene)
denpgs <- sapply(denpgs, function(x)paste0('ENSDARG', str_split_i(x, 'ENSDARG',2)))
denpgs <- z_gene_info$Symbol[sapply(denpgs, function(x)which(z_gene_info$zebrafish_ens==x))]
selected_LR_pairs <-zcellchat_db$interaction[(zcellchat_db$interaction$ligand %in% denpgs | zcellchat_db$interaction$receptor %in% denpgs),c('ligand','receptor','interaction_name')]
selected_LR_pairs['NPY_NPY8AR_1',]$interaction_name <- 'NPY_NPY8AR_1'
##2.2 make combination of cluster and DEnpgs 
df_temp <- de_target_genes
df_temp$gene <- sapply(df_temp$gene, function(x)paste0('ENSDARG', str_split_i(x, 'ENSDARG',2)))
df_temp$gene <- sapply(df_temp$gene, function(x)z_gene_info$Symbol[which(z_gene_info$zebrafish_ens==x)]) 
df_temp$cluster<-as.integer(df_temp$cluster)
targets <- paste0(df_temp$cluster[df_temp$cluster<=60],'_', df_temp$gene[df_temp$cluster<=60])

##2.3 cellchat object construct for SDA data
###2.3.1. create a cellchat object using a seurat object
Idents(gex.combined.hicat) <-c('cluster')
cellchat_sda <- createCellChat(object = subset(gex.combined.hicat,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07")&category=='SDA'), group.by = "cluster", assay = "RNA")
cellchat_sda@DB <- zcellchat_db
###2.3.2. change gene id to be consistent with cellchat database
#<note>: gene id in cellchat is in the format "gene_symbol-ensembl_id" as it is inherited from the seurat object
#<note>: the change is to keep gene_symbol only and consistent with that in gene_info dataset in cellchat database
####2.3.2.1 make a function for conversion
gene_id_convert<-function(old_id){
  if(old_id %in% paste0(z_gene_info$Symbol, '-', z_gene_info$zebrafish_ens)){
    ##find genes with the gene symbol and Ensembl id matched in gene_info dataset
    ##change the gene_id to be the symbol corresponding to the Ensembl id based on gene_info dataset
    new_id = z_gene_info$Symbol[which(paste0(z_gene_info$Symbol, '-', z_gene_info$zebrafish_ens)==old_id)][1]
  }else if(paste0('ENSDARG', str_split_i(old_id, 'ENSDARG',2)) %in% z_gene_info$zebrafish_ens){
    ##find genes with Ensembl id mismatching the gene symbol in gene_info dataset
    ##change the gene_id to be the symbol corresponding to the Ensembl id based on gene_info dataset
    new_id = z_gene_info$Symbol[which(z_gene_info$zebrafish_ens==paste0('ENSDARG', str_split_i(old_id, 'ENSDARG',2)))][1]
  }else{
    ##keep all other genes as they are. these are probably non-conserved genes that will not be used later
    new_id=old_id
    }
  return (new_id)
}
####2.3.2.2 apply the function for converions
rownames(cellchat_sda@data) <- unlist(lapply(rownames(cellchat_sda@data), function(x)gene_id_convert(x)))
####2.3.2.3 manually edit the gene npy8ar as it is a non-conserved gene
rownames(cellchat_sda@data)[grepl('npy8ar', rownames(cellchat_sda@data))]<-'npy8ar'
###2.3.3 subset the expression data of signaling genes for saving computation cost 
cellchat_sda <- subsetData(cellchat_sda)
###2.3.4 drop the levels of idents in cellchat object in case some idents disappear after subset the cellchat object
cellchat_sda@idents<-droplevels(cellchat_sda@idents, exclude = setdiff(levels(cellchat_sda@idents),unique(cellchat_sda@idents)))
###2.3.5 identify over-expressed ligands or receptors in each cluster
future::plan("sequential")
cellchat_sda <- identifyOverExpressedGenes(cellchat_sda)
options(future.globals.maxSize = 2000 * 1024^2)
###2.3.6 identify over-expressed LR in each of cluster with over-expressed ligands or receptor
cellchat_sda <- identifyOverExpressedInteractions(cellchat_sda)
options(future.globals.maxSize = 10000 * 1024^2)
###2.3.7 compute significance of each LR pair in each cluster
#<note>: the threshold for the percent of cells expressing the ligand or receptor gene per cluster is 4%
cellchat_sda <- computeCommunProb(cellchat_sda, type = "truncatedMean", trim = 0.04)
###2.3.8 calculates the aggregated cell-cell communication network
cellchat_sda <- aggregateNet(cellchat_sda)
#saveRDS(cellchat_sda, './processed_data/CellChat/cellchat_sda_v1.rds')

###2.3.9 get the dataframe consisting of all the inferred CCCs at the level of ligands/receptors 
df_net_sda <- subsetCommunication(cellchat_sda)
###2.3.10. Subset SDA ccc results for DEnpgs in relevant clusters
cellchat_sda_filter <- cellchat_sda
#<note>: clusters with less than 20 cells were excluded 
cellchat_sda_filter <- filterCommunication(cellchat_sda_filter, min.cells=20)
df_net_sda_denpg <- subsetCommunication(cellchat_sda_filter,pairLR.use=selected_LR_pairs)
df_net_sda_denpg <- df_net_sda_denpg %>% drop_na()
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
df_sda_ccc <-df_net_sda_denpg
df_sda_ccc$selected<-paste0(df_sda_ccc$source, '_', df_sda_ccc$ligand)
df_sda_ccc <- df_sda_ccc[df_sda_ccc$selected%in%targets,]
#<note>: non-neuronal clusters were excluded
df_sda_ccc<-df_sda_ccc[df_sda_ccc$target<=60,]
#saveRDS(df_sda_ccc, './processed_data/CellChat/df_sda_ccc_v1.rds')

#2.3.11 count ccc in SDA and plot in a heatmap 
df_lr_count_sda <- as.data.frame(matrix(0, nrow = 83, ncol = 83))
rownames(df_lr_count_sda) <- c(1:83)
colnames(df_lr_count_sda) <- c(1:83)
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_sda[i,j] <- dim(df_sda_ccc[(df_sda_ccc$source==i)&(df_sda_ccc$target==j),])[1]
  }
}
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
temp_list<- temp_list[temp_list<=60]
rownames(df_lr_count_sda) <- sapply(rownames(df_lr_count_sda), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_sda) <- sapply(colnames(df_lr_count_sda), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
pheatmap(df_lr_count_sda[temp_list,1:60],cluster_rows=F, cluster_cols=F,
         color = c('#000000', colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)),
         breaks = seq(0, 10, length.out = 100),legend = FALSE,
         show_rownames = F, show_colnames = T)

df_lr_count_sda_vdenpg <- as.data.frame(matrix(0, nrow = 83, ncol = 83))
rownames(df_lr_count_sda_vdenpg) <- c(1:83)
colnames(df_lr_count_sda_vdenpg) <- c(1:83)
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_sda_vdenpg[i,j] <- dim(df_sda_ccc[(df_sda_ccc$source==i)&(df_sda_ccc$target==j)&(df_sda_ccc$ligand %in% validated_depngs),])[1]
  }
}
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
temp_list<- temp_list[temp_list<=60]
rownames(df_lr_count_sda_vdenpg) <- sapply(rownames(df_lr_count_sda_vdenpg), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_sda_vdenpg) <- sapply(colnames(df_lr_count_sda_vdenpg), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
pheatmap(df_lr_count_sda_vdenpg[temp_list,1:60],cluster_rows=F, cluster_cols=F,
         color = c('#000000', colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)),
         breaks = seq(0, 10, length.out = 100),legend = FALSE,
         show_rownames = F, show_colnames = T)

##2.4 cellchat object construct for VDA data
###2.4.1. create a cellchat object using a seurat object
cellchat_vda <- createCellChat(object = subset(gex.combined.hicat,subset=orig.ident %in% c("GEX_s01","GEX_s02","GEX_s03","GEX_s05","GEX_s06","GEX_s07")&category=='VDA'), group.by = "cluster", assay = "RNA")
cellchat_vda@DB <- zcellchat_db
###2.4.2. change gene id to be consistent with cellchat database
####2.4.2.1 apply the function for converions
rownames(cellchat_vda@data) <- unlist(lapply(rownames(cellchat_vda@data), function(x)gene_id_convert(x)))
####2.4.2.2 manually edit the gene npy8ar as it is a non-conserved gene
rownames(cellchat_vda@data)[grepl('npy8ar', rownames(cellchat_vda@data))]<-'npy8ar'
###2.4.3 subset the expression data of signaling genes for saving computation cost
cellchat_vda <- subsetData(cellchat_vda)
###2.4.4 drop the levels of idents in cellchat object in case some idents disappear after subset the cellchat object
cellchat_vda@idents<-droplevels(cellchat_vda@idents, exclude = setdiff(levels(cellchat_vda@idents),unique(cellchat_vda@idents)))
###2.4.5 identify over-expressed ligands or receptors in each cluster
future::plan("sequential")
cellchat_vda <- identifyOverExpressedGenes(cellchat_vda)
###2.4.6 identify over-expressed LR in each of cluster with over-expressed ligands or receptor
options(future.globals.maxSize = 2000 * 1024^2)
cellchat_vda <- identifyOverExpressedInteractions(cellchat_vda)
options(future.globals.maxSize = 10000 * 1024^2)
###2.4.7 compute significance of each LR pair in each cluster
#<note>: the threshold for the percent of cells expressing the ligand or receptor gene per cluster is 4%
cellchat_vda <- computeCommunProb(cellchat_vda, type = "truncatedMean", trim = 0.04)
###2.4.8 calculates the aggregated cell-cell communication network
cellchat_vda <- aggregateNet(cellchat_vda)
#saveRDS(cellchat_vda, './processed_data/CellChat/cellchat_vda_v1.rds')
###2.4.9 get the dataframe consisting of all the inferred CCCs at the level of ligands/receptors
df_net_vda <- subsetCommunication(cellchat_vda)
###2.4.10. Subset SDA ccc results for DEnpgs in relevant clusters
cellchat_vda_filter <- cellchat_vda
#<note>: clusters with less than 20 cells were excluded
cellchat_vda_filter <- filterCommunication(cellchat_vda_filter, min.cells=20)
df_net_vda_denpg <- subsetCommunication(cellchat_vda_filter,pairLR.use=selected_LR_pairs)
df_net_vda_denpg <- df_net_vda_denpg %>% drop_na()
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
df_vda_ccc <-df_net_vda_denpg
df_vda_ccc$selected<-paste0(df_vda_ccc$source, '_', df_vda_ccc$ligand)
df_vda_ccc <- df_vda_ccc[df_vda_ccc$selected%in%targets,]
#<note>: non-neuronal clusters were excluded
df_vda_ccc<-df_vda_ccc[df_vda_ccc$target<=60,]
#saveRDS(df_vda_ccc, './processed_data/CellChat/df_vda_ccc_v1.rds')

#2.4.11 count ccc in VDA and plot in a heatmap 
df_lr_count_vda <- as.data.frame(matrix(0, nrow = 83, ncol = 83))
rownames(df_lr_count_vda) <- c(1:83)
colnames(df_lr_count_vda) <- c(1:83)
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_vda[i,j] <- dim(df_vda_ccc[(df_vda_ccc$source==i)&(df_vda_ccc$target==j),])[1]
  }
}
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
temp_list<- temp_list[temp_list<=60]
rownames(df_lr_count_vda) <- sapply(rownames(df_lr_count_vda), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_vda) <- sapply(colnames(df_lr_count_vda), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
pheatmap(df_lr_count_vda[temp_list,1:60],cluster_rows=F, cluster_cols=F,
         color = c('#000000', colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)),
         breaks = seq(0, 10, length.out = 100),legend = FALSE,
         show_rownames = F, show_colnames = T)

df_lr_count_vda_vdenpg <- as.data.frame(matrix(0, nrow = 83, ncol = 83))
rownames(df_lr_count_vda_vdenpg) <- c(1:83)
colnames(df_lr_count_vda_vdenpg) <- c(1:83)
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_vda_vdenpg[i,j] <- dim(df_vda_ccc[(df_vda_ccc$source==i)&(df_vda_ccc$target==j)&(df_vda_ccc$ligand %in% validated_depngs),])[1]
  }
}
temp_list<- sort(as.integer(unique(str_split_i(targets, '_',1))))
temp_list<- temp_list[temp_list<=60]
rownames(df_lr_count_vda_vdenpg) <- sapply(rownames(df_lr_count_vda_vdenpg), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_vda_vdenpg) <- sapply(colnames(df_lr_count_vda_vdenpg), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
pheatmap(df_lr_count_vda_vdenpg[temp_list,1:60],cluster_rows=F, cluster_cols=F,
         color = c('#000000', colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)),
         breaks = seq(0, 10, length.out = 100),legend = FALSE,
         show_rownames = F, show_colnames = T)

##############Chapter 3. Comparative analysis between SDA and VDA####
##3.1 Venn Diagram of LR pairs shared or distinct 
###3.1.1 LR pairs with all DEnpgs
library(eulerr)
list_data <- list(
  VDA = unique(df_vda_ccc$interaction_name),
  SDA = unique(df_sda_ccc$interaction_name)
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', quantities=TRUE)

###3.1.2 LR pairs with validated DEnpgs
validated_depngs <- c('npy','penkb','crhb','pyya','pyyb','sst7')
list_data <- list(
  VDA = unique(paste0(df_vda_ccc$ligand[df_vda_ccc$ligand %in% validated_depngs], '-', df_vda_ccc$receptor[df_vda_ccc$ligand %in% validated_depngs])),
  SDA = unique(paste0(df_sda_ccc$ligand[df_sda_ccc$ligand %in% validated_depngs], '-', df_sda_ccc$receptor[df_sda_ccc$ligand %in% validated_depngs]))
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', 
     quantities=TRUE)

###3.1.3 LR pairs with non-validated DEnpgs
list_data <- list(
  VDA = unique(paste0(df_vda_ccc$ligand[!df_vda_ccc$ligand %in% validated_depngs], '-', df_vda_ccc$receptor[!df_vda_ccc$ligand %in% validated_depngs])),
  SDA = unique(paste0(df_sda_ccc$ligand[!df_sda_ccc$ligand %in% validated_depngs], '-', df_sda_ccc$receptor[!df_sda_ccc$ligand %in% validated_depngs]))
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', 
     quantities=TRUE)


##3.2 Venn Diagram of cluster pairs shared or distinct 
###3.2.1  cluster pairs with all DEnpgs
list_data <- list(
  VDA = unique(paste0(df_vda_ccc$source, '-', df_vda_ccc$target)),
  SDA = unique(paste0(df_sda_ccc$source, '-', df_sda_ccc$target))
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', 
     quantities=TRUE)
###3.2.2  cluster pairs with validated DEnpgs
validated_depngs <- c('npy','penkb','crhb','pyya','pyyb','sst7')
list_data <- list(
  VDA = unique(paste0(df_vda_ccc$source[df_vda_ccc$ligand %in% validated_depngs], '-', df_vda_ccc$target[df_vda_ccc$ligand %in% validated_depngs])),
  SDA = unique(paste0(df_sda_ccc$source[df_sda_ccc$ligand %in% validated_depngs], '-', df_sda_ccc$target[df_sda_ccc$ligand %in% validated_depngs]))
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', 
     quantities=TRUE)

###3.2.3  cluster pairs with non-validated DEnpgs
validated_depngs <- c('npy','penkb','crhb','pyya','pyyb','sst7')
list_data <- list(
  VDA = unique(paste0(df_vda_ccc$source[!df_vda_ccc$ligand %in% validated_depngs], '-', df_vda_ccc$target[!df_vda_ccc$ligand %in% validated_depngs])),
  SDA = unique(paste0(df_sda_ccc$source[!df_sda_ccc$ligand %in% validated_depngs], '-', df_sda_ccc$target[!df_sda_ccc$ligand %in% validated_depngs]))
)
fit <- euler(list_data)
plot(fit, fills = c("#2C5F8A", "#C06C84"), labels = '', 
     quantities=TRUE)

##3.3 heatmap for cluster pairs shared or distinct
##3.3.1 heatmap for cluster pairs using all DEnpgs
df_lr_count_sv <- as.data.frame(matrix(0, nrow = length(temp_list), ncol = dim(df_lr_count_vda[temp_list,1:60])[2]))
rownames(df_lr_count_sv) <- c(1:dim(df_lr_count_sv)[1])
colnames(df_lr_count_sv) <- c(1:dim(df_lr_count_sv)[2])

for(i in temp_list){
  for(j in c(1:dim(df_lr_count_vda[temp_list,1:60])[2])){
    if(df_lr_count_vda[i, j]>0&df_lr_count_sda[i, j]>0){df_lr_count_sv[i,j]=3}
    else if(df_lr_count_vda[i, j]>0&df_lr_count_sda[i, j]==0){df_lr_count_sv[i,j]=1}
    else if(df_lr_count_vda[i, j]==0&df_lr_count_sda[i, j]>0){df_lr_count_sv[i,j]=2}
    else{df_lr_count_sv[i,j]=0}
  }
}
df_lr_count_sv <- df_lr_count_sv[temp_list,]
pheatmap(df_lr_count_sv,cluster_rows=F, cluster_cols=F,
         color = c('white', "#2C5F8A", "#C06C84", "#FFC000"),
         legend = FALSE,
         show_rownames = F, show_colnames = T)

##3.3.2 heatmap for cluster pairs using validated DEnpgs
df_lr_count_sv_vdenpg <- as.data.frame(matrix(0, nrow = length(temp_list), ncol = dim(df_lr_count_vda_vdenpg[temp_list,1:60])[2]))
rownames(df_lr_count_sv_vdenpg) <- c(1:dim(df_lr_count_sv_vdenpg)[1])
colnames(df_lr_count_sv_vdenpg) <- c(1:dim(df_lr_count_sv_vdenpg)[2])
for(i in temp_list){
  for(j in c(1:dim(df_lr_count_sv_vdenpg[temp_list,1:60])[2])){
    if(df_lr_count_vda_vdenpg[i, j]>0&df_lr_count_sda_vdenpg[i, j]>0){df_lr_count_sv_vdenpg[i,j]=3}
    else if(df_lr_count_vda_vdenpg[i, j]>0&df_lr_count_sda_vdenpg[i, j]==0){df_lr_count_sv_vdenpg[i,j]=1}
    else if(df_lr_count_vda_vdenpg[i, j]==0&df_lr_count_sda_vdenpg[i, j]>0){df_lr_count_sv_vdenpg[i,j]=2}
    else{df_lr_count_sv_vdenpg[i,j]=0}
  }
}
df_lr_count_sv_vdenpg <- df_lr_count_sv_vdenpg[temp_list,]
pheatmap(df_lr_count_sv_vdenpg,cluster_rows=F, cluster_cols=F,
         color = c('white', "#2C5F8A", "#C06C84", "#FFC000"),
         legend = FALSE,
         show_rownames = F, show_colnames = T)

##3.4. chord diagram for one selected DEnpg in one selected cluster
##3.4.1 prepare count matrix for VDA
df_lr_count_onegene <- as.data.frame(matrix(NA, nrow = 83, ncol = 83))
rownames(df_lr_count_onegene) <- c(1:83)
colnames(df_lr_count_onegene) <- c(1:83)
target_gene <- c('crhb')
#target_receptor <- c('npy1r', 'npy8ar', 'npy2rl','gpr83')
target_receptor <- c('crhr2')
target_cluster <- de_target_genes$cluster[sapply(de_target_genes$gene, function(x)str_split_i(x, '-ENSDARG',1) %in% target_gene)]
df_temp <- df_vda_ccc[(df_vda_ccc$ligand %in% target_gene)&(df_vda_ccc$receptor %in% target_receptor),]
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_onegene[i,j] <- dim(df_temp[(df_temp$source==i)&(df_temp$target==j),])[1]
  }
}
rownames(df_lr_count_onegene) <- sapply(rownames(df_lr_count_onegene), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_onegene) <- sapply(colnames(df_lr_count_onegene), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
df_temp <-df_lr_count_onegene[target_cluster,]
df_temp <- df_temp[, colSums(df_temp)>0]
df_temp_vda<-df_temp

##3.4.2 prepare count matrix for SDA
df_lr_count_onegene <- as.data.frame(matrix(NA, nrow = 83, ncol = 83))
rownames(df_lr_count_onegene) <- c(1:83)
colnames(df_lr_count_onegene) <- c(1:83)
df_temp <- df_sda_ccc[(df_sda_ccc$ligand %in% target_gene)&(df_sda_ccc$receptor %in% target_receptor),]
for(i in c(1:83)){
  for(j in c(1:83)){
    df_lr_count_onegene[i,j] <- dim(df_temp[(df_temp$source==i)&(df_temp$target==j),])[1]
  }
}
rownames(df_lr_count_onegene) <- sapply(rownames(df_lr_count_onegene), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
colnames(df_lr_count_onegene) <- sapply(colnames(df_lr_count_onegene), function(x)ifelse(as.integer(x)<10, paste0('00', as.character(x)), paste0('0', as.character(x))))
df_temp <-df_lr_count_onegene[target_cluster,]
df_temp <- df_temp[, colSums(df_temp)>0]
df_temp_sda<-df_temp
related_clusters <- sort(unique(c(rownames(df_temp_sda), colnames(df_temp_sda), rownames(df_temp_vda), colnames(df_temp_vda))))

##3.4.3 combined the data for a consensus of clusters used for plotting
###3.4.3.1. reform count matrix for sda
mat_temp <- matrix(0, nrow = length(related_clusters), ncol = length(related_clusters),dimnames = list(related_clusters, related_clusters))
mat_temp[rownames(df_temp_sda['010',]),colnames(df_temp_sda)]<-as.matrix(df_temp_sda['010',])
###3.4.3.2. setup the colors
mycolor <- c()
for(cluster_id in rownames(mat_temp)){
  mycolor[cluster_id] = cluster_annotation[cluster_annotation$cluster==cluster_id, 'color'][[1]]
}
mycolor <- unlist(mycolor)
###3.4.3.3 plot chord diagram
netVisual_circle(mat_temp, weight.scale = T, label.edge= F, title.name = "Number of interactions",
                 color.use = mycolor, alpha.edge = 1)

###3.4.3.4. reform count matrix for vda
mat_temp <- matrix(0, nrow = length(related_clusters), ncol = length(related_clusters),dimnames = list(related_clusters, related_clusters))
mat_temp[rownames(df_temp_vda['010',]),colnames(df_temp_vda)]<-as.matrix(df_temp_vda['010',])
###3.4.3.5. setup the colors
mycolor <- c()
for(cluster_id in rownames(mat_temp)){
  mycolor[cluster_id] = cluster_annotation[cluster_annotation$cluster==cluster_id, 'color'][[1]]
}
mycolor <- unlist(mycolor)
###3.4.3.6 plot chord diagram
netVisual_circle(mat_temp, weight.scale = T, label.edge= F, title.name = "Number of interactions",
                 color.use = mycolor, alpha.edge = 1)
