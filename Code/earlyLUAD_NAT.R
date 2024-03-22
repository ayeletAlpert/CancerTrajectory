library(readr)
library(Matrix)
library(Seurat)
require(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)
require(stringr)
require(ggplot2)
library(BayesPrism)
library(devtools)
library(matrixStats)
library(SummarizedExperiment)
library(InstaPrism)
library(survival)
library(ComplexHeatmap)
library(clusterProfiler)
library(TCGAbiolinks)
library(SummarizedExperiment)
# library(biomaRt)
library(e1071)
library(KEGGREST)
library(httr)
library(GSEABase)
require(reshape2)
library(dplyr)
library(glmnet)
library(pROC)
library(umap)
library(GEOquery)
library(qusage)
theme_set(theme_bw())

data_dir <- "D:\\Ayelet\\GE_NAT\\LUAD\\GE_data"
GE_NAT <- read.csv(file = file.path(data_dir, "GSE229705_counts-normalized.csv"), row.names = 1, check.names = F)
GE_NAT <- GE_NAT[,str_detect(colnames(GE_NAT), "\\-N")]
colnames(GE_NAT) <- str_replace(colnames(GE_NAT), "\\-N", "")

# Identify constant columns
constant_columns <- apply(GE_NAT, 2, function(col) all(col == col[1]))

# Remove constant columns
GE_NAT_filtered <- GE_NAT[, , drop = FALSE][, !constant_columns]

#calculate IL-6 activation pathway:
hallmark_genes_path <- "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt"
hallmark_genes <- read.gmt(hallmark_genes_path)
inflammatory_pathways <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                           "HALLMARK_IL2_STAT5_SIGNALING")
tt_progression = do.call('rbind', lapply(inflammatory_pathways, function(inflammatory_pathway){
  print(inflammatory_pathway)
  pathway_genes <- unique(hallmark_genes$gene[hallmark_genes$term  == inflammatory_pathway])
  prcomp_res <- prcomp(t(GE_NAT_filtered[intersect(pathway_genes, rownames(GE_NAT_filtered)),]), scale. = T)
  biplot(prcomp_res)
  pca_res_pathway <- prcomp_res$x[,c("PC1", "PC2")]
  progression_data$pathway <- pca_res_pathway[rownames(progression_data),"PC1"]
  ggplot(progression_data, aes(x = progression, y = pathway)) + geom_boxplot()
  
  ttres <- t.test(progression_data$pathway[progression_data$progression == "Progression"],
                  progression_data$pathway[progression_data$progression == "No Progression (>5y)"])[[3]]
  
  return(data.frame(inflammatory_pathway = inflammatory_pathway, tt_p = ttres))
}))


#extract phenotypic data:
gse <- getGEO("GSE229705")[[1]]
progression_data <- pData(gse)[,str_detect(colnames(pData(gse)), "characteristics_ch1")]
progression_data <- progression_data[progression_data$characteristics_ch1.1 == "tissue: Lung adjacent normal",]
rownames(progression_data) <- str_replace(progression_data$characteristics_ch1, "patient: ", "")
progression_data$progression <- str_replace(progression_data$characteristics_ch1.3, "progression: ", "")
progression_data <- progression_data[intersect(rownames(progression_data), rownames(pca_res_il6)),]


data_proteomics <- read.table(file = "~/MDPhD/oncology/BayesPrismLUAD/onco_age/proteomic_data/pancancer_olink_data_biostudies_v2.txt", header = T)
data_proteomics$protein_name <- mapIds(org.Hs.eg.db, keys = data_proteomics$UniProt, column = "SYMBOL", keytype = "UNIPROT")
data_proteomics <- data_proteomics[complete.cases(data_proteomics),]
data_proteomics_wide <- dcast(formula = Sample_ID + Cancer ~ protein_name, value.var = "NPX", data = data_proteomics)
rownames(data_proteomics_wide) <- data_proteomics_wide[,"Sample_ID"]
data_proteomics_wide <- data_proteomics_wide[complete.cases(data_proteomics_wide),]

# solid_cancers <- c("BRC", "CVX", "CRC", "ENDC","GLIOM","LUNGC","OVC","PRC")
solid_cancers <- c("BRC", "CRC", "LUNGC", "PRC")
solid_cancers <- c("BRC", "CRC", "LUNGC", "PRC", "CVX", "ENDC", "GLIOM", "OVC")
data_proteomics <- data_proteomics[data_proteomics$Cancer %in% solid_cancers,]
protein_name <- c("IFNG", "IL6", "IL2", "TNF", "TGFB1")
ggplot(data_proteomics[data_proteomics$protein_name %in% protein_name,], aes(x = Cancer, y = NPX)) + 
  geom_boxplot() + geom_jitter(alpha = 0.2) + facet_wrap(~protein_name, scales = "free") + coord_flip() + 
  geom_hline(yintercept = 0, color = "red")
