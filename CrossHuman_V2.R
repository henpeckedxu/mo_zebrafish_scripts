##1. backgroud gene list with all human genes orthologus to zebrafish genes####
bg_genes_flt <- bg_genes[!grepl(":", bg_genes$zebrafish_symbol),]#remove genes without a proper gene name (predicted genes, BAC clone genes, etc.)
bg_genes_flt <- bg_genes_flt[!is.na(bg_genes_flt$zebrafish_symbol),]
total_human_sym <- unique(bg_genes_flt$human_symbol)
total_human_sym <- total_human_sym[!is.na(total_human_sym)]

##2.query gene list of human orthologs for DEnpgs####
goea_de_hgene_all_pos <- readRDS('./processed_data/GOEA/goea_de_hgene_all_pos.rds')
target_terms<-c('GO:0005184','GO:0007218','GO:0098992','HPA:0260141','KEGG:04080')
hgenes<-unique(unlist(strsplit(goea_de_hgene_all_pos[goea_de_hgene_all_pos$term_id %in% target_terms,]$intersection,","), recursive = FALSE))
DEnpg_human_sym <- unique(bg_genes_flt[bg_genes_flt$human_ens%in%hgenes, ]$human_symbol)

##3. find overlap genes in scRNA-seq studies of four psychiatric disorders
##3.1 SCZ
scz_deg_all <- c()
for (cell_type in excel_sheets("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/SCZ_scRNA_PMID38781388_s4.xlsx")){
  scz_genes <- read_excel("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/SCZ_scRNA_PMID38781388_s4.xlsx", sheet =cell_type)
  scz_genes <- as.data.frame(scz_genes[(scz_genes$Meta_adj.P.Val<0.05)&(abs(scz_genes$Meta_logFC)>0.1), ])
  scz_genes$cell_type = cell_type
  scz_deg_all <- rbind(scz_deg_all, scz_genes)
}
scz_deg_all <- scz_deg_all[scz_deg_all$gene %in% total_human_sym, ]#dataset of all SCZ DEG orthologos 
scz_deg_all_in <- setdiff(unique(scz_deg_all[grepl("In-",scz_deg_all$cell_type),]$gene), unique(scz_deg_all[grepl("Ex-",scz_deg_all$cell_type),]$gene))#IN neuron specific SCZ DEG orthologos
scz_deg_all_ex <- setdiff(unique(scz_deg_all[grepl("Ex-",scz_deg_all$cell_type),]$gene), unique(scz_deg_all[grepl("In-",scz_deg_all$cell_type),]$gene))#EX neuron specific SCZ DEG orthologos
scz_deg_all_inex <- intersect(unique(scz_deg_all[grepl("In-",scz_deg_all$cell_type),]$gene), unique(scz_deg_all[grepl("Ex-",scz_deg_all$cell_type),]$gene))# SCZ orthologos with DE in both EX and IN neurons
scz_nondeg_all <- setdiff(total_human_sym, unique(scz_deg_all[grepl("In-|Ex-",scz_deg_all$cell_type),]$gene))#human orthologos without DE in EX or IN for SCZ

scz_DEnpg_all <- scz_deg_all[scz_deg_all$gene %in% DEnpg_human_sym, ]
scz_DEnpg_all_in <- setdiff(scz_DEnpg_all[grepl("In-",scz_DEnpg_all$cell_type),]$gene, scz_DEnpg_all[grepl("Ex-",scz_DEnpg_all$cell_type),]$gene)
scz_DEnpg_all_ex <- setdiff(scz_DEnpg_all[grepl("Ex-",scz_DEnpg_all$cell_type),]$gene, scz_DEnpg_all[grepl("In-",scz_DEnpg_all$cell_type),]$gene)
scz_DEnpg_all_inex <- intersect(scz_DEnpg_all[grepl("Ex-",scz_DEnpg_all$cell_type),]$gene, scz_DEnpg_all[grepl("In-",scz_DEnpg_all$cell_type),]$gene)
scz_nonDEnpg_all <- setdiff(DEnpg_human_sym, c(scz_DEnpg_all_in,scz_DEnpg_all_inex,scz_DEnpg_all_ex))

##3.2 ASD
asd_genes1 <- read_excel("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/ASD_scRNA_PMID38781372_s4.xlsx",sheet = 1,skip = 2)
asd_ex_deg <- asd_genes1[grepl("EXT_", asd_genes1$CELLTYPE),]
asd_ex_deg <- unique(unlist(strsplit(asd_ex_deg$`Genes (by rows for category)`, "\r\n")))
asd_ex_deg <- unique(unlist(strsplit(asd_ex_deg, " ")))
asd_ex_deg <- asd_ex_deg[!is.na(asd_ex_deg)]
asd_ex_deg <- asd_ex_deg[asd_ex_deg %in% total_human_sym]#overlap with all human orthologues of zebrafish genes

asd_in_deg <- asd_genes1[grepl("INT_", asd_genes1$CELLTYPE),]
asd_in_deg <- unique(unlist(strsplit(asd_in_deg$`Genes (by rows for category)`, "\r\n")))
asd_in_deg <- unique(unlist(strsplit(asd_in_deg, " ")))
asd_in_deg <- asd_in_deg[!is.na(asd_in_deg)]
asd_in_deg <- asd_in_deg[asd_in_deg %in% total_human_sym]#overlap with all human orthologues of zebrafish genes


asd_deg_all_in <- setdiff(asd_in_deg, asd_ex_deg)#IN neuron specific ASD DEG orthologos
asd_deg_all_ex <- setdiff(asd_ex_deg, asd_in_deg)#EX neuron specific ASD DEG orthologos
asd_deg_all_inex <- intersect(asd_ex_deg, asd_in_deg) #ASD DEG orthologos show DE in both IN and EX neurons
asd_deg_all <- c(asd_deg_all_in,asd_deg_all_ex,asd_deg_all_inex) # ASD DEG orthologos show DE in IN and/or EX neurons

asd_DEnpg_all <- intersect(asd_deg_all,DEnpg_human_sym)
asd_DEnpg_all_in <- intersect(asd_deg_all_in,DEnpg_human_sym)
asd_DEnpg_all_ex <- intersect(asd_deg_all_ex,DEnpg_human_sym)
asd_DEnpg_all_inex <- intersect(asd_deg_all_inex,DEnpg_human_sym)
asd_nonDEnpg_all <- setdiff(DEnpg_human_sym, c(asd_DEnpg_all_in,asd_DEnpg_all_inex,asd_DEnpg_all_ex))

##3.3 PTSD
mast_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_MAST.txt", sep='\t')##DEG test results using MAST
wilcox_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_Wilcox.txt", sep='\t')##DEG test results using Wilcox
celltypelist <- c('IN', 'LAMP5','KCNG1', 'VIP', 'SST','PVALB', 'EXN', 'CUX2','RORB','FEZF2','OPRK1')##select neuronal cell clusters 
ptsd_mdd_deg_all <- cbind(mast_res ,wilcox_res[,5:8])##combine the two stats test data
ptsd_mdd_deg_all <- ptsd_mdd_deg_all[ptsd_mdd_deg_all$Celltype %in% celltypelist, ]##subset the data for neuronal cell clusters
ptsd_deg_all <- ptsd_mdd_deg_all[(ptsd_mdd_deg_all$PTSD.MAST.FDR<0.05)&(ptsd_mdd_deg_all$PTSD.Wilcox.FDR<0.05),]##select DEGs with FDR<0.05 in both methods
ptsd_deg_all <- ptsd_deg_all[mapply(function(a,b)max(abs(a), abs(b))>log2(1.2), ptsd_deg_all$PTSD.MAST.log2FC, ptsd_deg_all$PTSD.Wilcox.log2FC),]##select DEGs with at least one FC > 1.2
ptsd_deg_all <- ptsd_deg_all[!is.na(ptsd_deg_all),]
ptsd_deg_all <- ptsd_deg_all[, !grepl('MDD.',colnames(ptsd_deg_all))]##drop columns unrelated to PTSD
ptsd_deg_all <- ptsd_deg_all[ptsd_deg_all$Genename %in% total_human_sym,]##select human genes only orthologous to zebrafish

ptsd_deg_all_in <- setdiff(ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[1:6],]$Genename, ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[7:11],]$Genename)#IN neuron specific PTSD DEG orthologos
ptsd_deg_all_ex <- setdiff(ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[7:11],]$Genename, ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[1:6],]$Genename)#Ex neuron specific PTSD DEG orthologos
ptsd_deg_all_inex <- intersect(unique(ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[1:6],]$Genename), unique(ptsd_deg_all[ptsd_deg_all$Celltype %in% celltypelist[7:11],]$Genename))# PTSD orthologos with DE in both EX and IN neurons
ptsd_nondeg_all <- setdiff(total_human_sym, c(ptsd_deg_all_in,ptsd_deg_all_ex, ptsd_deg_all_inex))#human orthologos without DE in EX or IN for PTSD

ptsd_DEnpg_all_in <- intersect(ptsd_deg_all_in,DEnpg_human_sym)
ptsd_DEnpg_all_ex <- intersect(ptsd_deg_all_ex,DEnpg_human_sym)
ptsd_DEnpg_all_inex <- intersect(ptsd_deg_all_inex,DEnpg_human_sym)
ptsd_DEnpg_all <- c(ptsd_DEnpg_all_in,ptsd_DEnpg_all_inex,ptsd_DEnpg_all_ex)
ptsd_nonDEnpg_all <- setdiff(DEnpg_human_sym, ptsd_DEnpg_all)

##3.4 MDD
mast_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_MAST.txt", sep='\t')##DEG test results using MAST
wilcox_res <- read.csv("/Users/jialexu/Desktop/Project2GWAS-BehvaioralGenetics/experiments/Database/hpsy_genes/PTSD_MDD_Wilcox.txt", sep='\t')##DEG test results using Wilcox
celltypelist <- c('IN', 'LAMP5','KCNG1', 'VIP', 'SST','PVALB', 'EXN', 'CUX2','RORB','FEZF2','OPRK1')##select neuronal cell clusters 
ptsd_mdd_deg_all <- cbind(mast_res ,wilcox_res[,5:8])##combine the two stats test data
ptsd_mdd_deg_all <- ptsd_mdd_deg_all[ptsd_mdd_deg_all$Celltype %in% celltypelist, ]##subset the data for neuronal cell clusters
mdd_deg_all <- ptsd_mdd_deg_all[(ptsd_mdd_deg_all$MDD.MAST.FDR<0.05)&(ptsd_mdd_deg_all$MDD.Wilcox.FDR<0.05),]##select DEGs with FDR<0.05 for MDD using both methods
mdd_deg_all <- mdd_deg_all[mapply(function(a,b)max(abs(a), abs(b))>log2(1.2), mdd_deg_all$MDD.MAST.log2FC, mdd_deg_all$MDD.Wilcox.log2FC),]##select DEGs with at least one FC > 1.2
mdd_deg_all <- mdd_deg_all[!is.na(mdd_deg_all),]
mdd_deg_all <- mdd_deg_all[, !grepl('MDD.',colnames(mdd_deg_all))]##drop columns unrelated to PTSD
mdd_deg_all <- mdd_deg_all[mdd_deg_all$Genename %in% total_human_sym,]##select human genes only orthologous to zebrafish

mdd_deg_all_in <- setdiff(mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[1:6],]$Genename, mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[7:11],]$Genename)#IN neuron specific mdd DEG orthologos
mdd_deg_all_ex <- setdiff(mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[7:11],]$Genename, mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[1:6],]$Genename)#Ex neuron specific mdd DEG orthologos
mdd_deg_all_inex <- intersect(unique(mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[1:6],]$Genename), unique(mdd_deg_all[mdd_deg_all$Celltype %in% celltypelist[7:11],]$Genename))# mdd orthologos with DE in both EX and IN neurons
mdd_nondeg_all <- setdiff(total_human_sym, c(mdd_deg_all_in,mdd_deg_all_ex, mdd_deg_all_inex))#human orthologos without DE in EX or IN for mdd

mdd_DEnpg_all_in <- intersect(mdd_deg_all_in,DEnpg_human_sym)
mdd_DEnpg_all_ex <- intersect(mdd_deg_all_ex,DEnpg_human_sym)
mdd_DEnpg_all_inex <- intersect(mdd_deg_all_inex,DEnpg_human_sym)
mdd_DEnpg_all <- c(mdd_DEnpg_all_in,mdd_DEnpg_all_inex,mdd_DEnpg_all_ex)
mdd_nonDEnpg_all <- setdiff(DEnpg_human_sym, mdd_DEnpg_all)


library(ggvenn)
ggvenn(
  list(
    'SCZ'= c(scz_DEnpg_all_in, scz_DEnpg_all_ex, scz_DEnpg_all_inex),
    'ASD' = c(asd_DEnpg_all_in, asd_DEnpg_all_ex, asd_DEnpg_all_inex),
    'PTSD'= c(ptsd_DEnpg_all_in, ptsd_DEnpg_all_ex, ptsd_DEnpg_all_inex), 
    'MDD'= c(mdd_DEnpg_all_in, mdd_DEnpg_all_ex, mdd_DEnpg_all_inex)), 
  fill_color = c('#a559aa',"#59a89c", "#f0c571", "#e02b35"),
  stroke_size = 1, set_name_size = 10,text_size = 8,
  show_percentage =FALSE
)

intersect(intersect(intersect(c(scz_DEnpg_all_in, scz_DEnpg_all_ex, scz_DEnpg_all_inex), c(asd_DEnpg_all_in, asd_DEnpg_all_ex, asd_DEnpg_all_inex)),
                    c(ptsd_DEnpg_all_in, ptsd_DEnpg_all_ex, ptsd_DEnpg_all_inex)),c(mdd_DEnpg_all_in, mdd_DEnpg_all_ex, mdd_DEnpg_all_inex))
