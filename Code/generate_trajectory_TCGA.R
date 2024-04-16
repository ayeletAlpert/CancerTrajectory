# initialize settings
rm(list = ls())

#load the required packaes:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/initial_settings.R")

#source the tailored functions:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/my_functions.R")

data_dir_bulk_tissue <- "/media/chronos/Storage/ayelet/TCGA_bulk_GE"

#start with BLCA:

cancer_types_short <- c("BLCA", "COAD", "READ", "KIRP", "KIRC", "LUAD", "LUSC", "READ")
cancer_type = "BLCA"

cancer_type_GE_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0(cancer_type, "_TCGA.rds")))
rownames(cancer_type_GE_data) <- sub("\\.\\d+$", "", rownames(cancer_type_GE_data))

#extract only those samples that are related to cancer patients
rel_samples <- colnames(cancer_type_GE_data)[cancer_type_GE_data$definition == "Primary solid Tumor" & !is.na(cancer_type_GE_data$ajcc_pathologic_stage)]

#sort the stages alphabetically:
stages <- gsub("([ABC])$", "", colData(cancer_type_GE_data)[rel_samples,"ajcc_pathologic_stage"])

# Create a sorted vector of unique stages
unique_stages <- sort(unique(stages))

#get stage-associated genes:
stage_genes <- readRDS(file = "/media/chronos/Storage/ayelet/saved_data/cancer_trajectory/sig_genes_stage_lm_per_cancer.rds")
cancer_specific_stage_genes <- stage_genes[[cancer_type]]

#generate trajectory:

#PCA:
gene_expression_stage_genes <- assays(cancer_type_GE_data)$tpm_unstrand[intersect(sym2ensembl(cancer_specific_stage_genes), rownames(cancer_type_GE_data)), 
                                                                        rel_samples]
pca_res <- data.frame(prcomp(t(gene_expression_stage_genes), scale. = T)$x, stage = stages)
ggplot(subset(pca_res, PC1 < 50), aes(x = PC1, y = PC2, color = stage)) + geom_point() + theme_minimal() + ggtitle()
ggplot(subset(pca_res, PC1 < 50), aes(x = PC1, y = ..density.., fill = stage)) + geom_density(alpha = 0.2) + theme_minimal()

#show the dynamics of different genes along PC1:
gene <- "ATL2"
pca_res_gene <- data.frame(pca_res, gene_exp = assays(cancer_type_GE_data)$tpm_unstrand[sym2ensembl(gene), rownames(pca_res)])
ggplot(subset(pca_res_gene, PC1 < 50), aes(x = PC1, y = gene_exp)) + geom_point() + geom_smooth() + theme_minimal() + ggtitle(gene)

total_genes <- nrow(cancer_type_GE_data)totacancer_typel_genes <- nrow(cancer_type_GE_data)
p_vals_genes_stages <- do.call('rbind', lapply(1:total_genes, function(i){
  if(sum(assays(cancer_type_GE_data)$tpm_unstrand[i, rel_samples] == 0) == length(rel_samples)){return(NULL)}
  lm_res <- summary(lm(assays(cancer_type_GE_data)$tpm_unstrand[i, rel_samples] ~ numerical_stages))
  # plot(cancer_type_RnaseqSE_stage, assays(exp_data)$tpm_unstrand[gene,rel_samples_healthy])
  
  # Calculate progress percentage
  progress <- round((i / total_genes) * 100)
  if (progress %% progress_interval == 0) {
    cat("\rProgress:", progress, "% ", sep = "")
  }
  
  return(data.frame(gene_ensembl = rownames(cancer_type_GE_data)[i],
                    #gene_symbol = ensembl2sym(rownames(cancer_type_GE_data)[i]),
                    p_val = coef(lm_res)[2,4]))
}))
