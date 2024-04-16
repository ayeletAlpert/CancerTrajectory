# initialize settings
rm(list = ls())

#load the required packaes:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/initial_settings.R")

#source the tailored functions:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/my_functions.R")

data_dir_bulk_tissue <- "/media/chronos/Storage/ayelet/TCGA_bulk_GE"

cancer_types <- list.files(path = data_dir_bulk_tissue)
exclude_cancers <- c("GBM_TCGA.rds", "OV_TCGA.rds", "PRAD_TCGA.rds", "SARC_TCGA.rds", "UCEC_TCGA.rds")
cancer_types <- cancer_types[!(cancer_types %in% exclude_cancers)]

#per cancer type, get the number of samples per stage:
stages_per_cancer <- do.call('rbind', lapply(cancer_types, function(cancer_type){
  message(cancer_type)
  cancer_type_GE_data <- readRDS(file = file.path(data_dir_bulk_tissue, cancer_type))
  rownames(cancer_type_GE_data) <- sub("\\.\\d+$", "", rownames(cancer_type_GE_data))
  
  #extract only those samples that are related to cancer patients
  rel_samples <- colnames(cancer_type_GE_data)[cancer_type_GE_data$definition == "Primary solid Tumor"]
  rel_stages <- colData(cancer_type_GE_data)[rel_samples,"ajcc_pathologic_stage"]
  rel_stages <- rel_stages[!is.na(rel_stages)]
  rel_stages <- gsub("([ABC])$", "", rel_stages)
  return(data.frame(cancer_type = sub("\\_TCGA\\.rds$", "", cancer_type), rel_stages = rel_stages))
}))

# Define colors for each stage
stage_colors <- c("Stage I" = "blue", "Stage II" = "green", "Stage III" = "orange", "Stage IV" = "red")

ggplot(stages_per_cancer, aes(x = rel_stages, fill = rel_stages)) +
  geom_bar(position = "stack") +
  facet_wrap(~ cancer_type, scales = "free") +
  scale_fill_manual(values = stage_colors) +
  labs(x = "Stage", y = "Count", title = "Distribution of Cancer Stages by Cancer Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Identify per cancer the DE genes over stages:
cancer_types_short <- c("BLCA", "COAD", "READ", "KIRP", "KIRC", "LUAD", "LUSC", "READ")
total_cancer_types <- length(cancer_types_short)
progress_interval <- 1  # Adjust this value as needed

DE_genes_Xross_stages <- do.call('rbind', lapply(paste0(cancer_types_short, "_TCGA.rds"), function(cancer_type){
  cat(paste("\nReading data for:", cancer_type, "\n"))
  cancer_type_GE_data <- readRDS(file = file.path(data_dir_bulk_tissue, cancer_type))
  rownames(cancer_type_GE_data) <- sub("\\.\\d+$", "", rownames(cancer_type_GE_data))
  
  #extract only those samples that are related to cancer patients
  rel_samples <- colnames(cancer_type_GE_data)[cancer_type_GE_data$definition == "Primary solid Tumor" & !is.na(cancer_type_GE_data$ajcc_pathologic_stage)]
  
  #sort the stages alphabetically:
  stages <- gsub("([ABC])$", "", colData(cancer_type_GE_data)[rel_samples,"ajcc_pathologic_stage"])
  
  # Create a sorted vector of unique stages
  unique_stages <- sort(unique(stages))
  
  # Create a mapping dictionary
  mapping <- setNames(1:length(unique_stages), unique_stages)
  
  # Convert stages to numerical values
  numerical_stages <- sapply(stages, function(stage) mapping[stage])
  
  total_genes <- nrow(cancer_type_GE_data)
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
  
  return(data.frame(cancer_type = sub("\\_TCGA\\.rds$", "", cancer_type), p_vals_genes_stages))
}))
ggplot(DE_genes_Xross_stages, aes(x = p_val)) + geom_histogram() + facet_wrap(~cancer_type)
DE_genes_Xross_stages$gene_symbol <- ensembl2sym(DE_genes_Xross_stages$gene_ensembl)
cat("\n")  # Add a new line after progress is completed

# DE genes Xross tumor ----------------------------------------------------
DE_genes_Xross_T <- do.call('rbind', lapply(paste0(cancer_types_short, "_TCGA.rds"), function(cancer_type){
  cat(paste("\nReading data for:", cancer_type, "\n"))
  cancer_type_GE_data <- readRDS(file = file.path(data_dir_bulk_tissue, cancer_type))
  rownames(cancer_type_GE_data) <- sub("\\.\\d+$", "", rownames(cancer_type_GE_data))
  
  #extract only those samples that are related to cancer patients
  rel_samples <- colnames(cancer_type_GE_data)[cancer_type_GE_data$definition == "Primary solid Tumor" & !is.na(cancer_type_GE_data$ajcc_pathologic_stage)]
  
  #sort the stages alphabetically:
  stages <- gsub("([abc])$", "", colData(cancer_type_GE_data)[rel_samples,"ajcc_pathologic_t"])
  
  # Create a sorted vector of unique stages
  unique_stages <- sort(unique(stages))
  
  # Create a mapping dictionary
  mapping <- setNames(1:length(unique_stages), unique_stages)
  
  # Convert stages to numerical values
  numerical_stages <- sapply(stages, function(stage) mapping[stage])
  
  total_genes <- nrow(cancer_type_GE_data)
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
  
  return(data.frame(cancer_type = sub("\\_TCGA\\.rds$", "", cancer_type), p_vals_genes_stages))
}))
DE_genes_Xross_T$gene_symbol <- ensembl2sym(DE_genes_Xross_T$gene_ensembl)

#if this is not TRUE, DO NOT PROCEED
if(identical(p_vals_genes_stages$gene_ensembl, DE_genes_Xross_stages$gene_ensembl[DE_genes_Xross_stages$cancer_type == "BLCA"])){
  comp_DE_T_stage <- data.frame(gene_ensemble = p_vals_genes_stages$gene_ensembl, 
                                p_by_stage = DE_genes_Xross_stages$p_val[DE_genes_Xross_stages$cancer_type == "BLCA"], 
                                p_by_T = p_vals_genes_stages$p_val)
}
ggplot(comp_DE_T_stage, aes(x = -log(p_by_stage), y = -log(p_by_T))) + geom_point() + theme_bw()
heatscatter(-log(comp_DE_T_stage$p_by_stage), y = -log(comp_DE_T_stage$p_by_T), main = "-logP stage VS. -logP T")

#gene-selection:
DE_genes_Xross_stages <- DE_genes_Xross_stages[!is.na(DE_genes_Xross_stages$gene_ensembl) & !is.na(DE_genes_Xross_stages$gene_symbol),]
DE_genes_Xross_stages$p_adj <- unlist(lapply(unique(DE_genes_Xross_stages$cancer_type), function(cancer_type){
  return(p.adjust(DE_genes_Xross_stages$p_val[DE_genes_Xross_stages$cancer_type == cancer_type], method = "BH"))
}))

sig_genes_per_cancer <- lapply(unique(DE_genes_Xross_stages$cancer_type), function(cancer_type){
  return(DE_genes_Xross_stages$gene_symbol[DE_genes_Xross_stages$p_adj < 0.05 & DE_genes_Xross_stages$cancer_type == cancer_type])
})
names(sig_genes_per_cancer) <- unique(DE_genes_Xross_stages$cancer_type)
num_sig_genes_per_cancer <- do.call('rbind', lapply(unique(DE_genes_Xross_stages$cancer_type), function(cancer_type){
  return(data.frame(cancer_type = cancer_type, sum_sig_genes = length(sig_genes_per_cancer[[cancer_type]])))}))
ggplot(num_sig_genes_per_cancer, aes(x = cancer_type, y = sum_sig_genes)) + geom_bar(stat = "identity") + theme_minimal() +
  ggtitle('number of significant (BH-corrected P < 0.05) per cancer type')
saveRDS(sig_genes_per_cancer, file = "/media/chronos/Storage/ayelet/saved_data/cancer_trajectory/sig_genes_stage_lm_per_cancer.rds")

#get the UBC longitudinal model to check the seed features change along TCGA disease stages:
model <- readRDS(file = "/media/chronos/Storage/ayelet/saved_data/cancer_trajectory/model_UBC_longitudinal.rds")

seed_genes <- model$seed

BLCA_p_stage_T <- data.frame(gene_symbol = DE_genes_Xross_T$gene_symbol[DE_genes_Xross_T$cancer_type == "BLCA"],
                             gene_in_seed = (DE_genes_Xross_T$gene_symbol[DE_genes_Xross_T$cancer_type == "BLCA"]) %in% seed_genes,
                             p_val_T = DE_genes_Xross_T$p_val[DE_genes_Xross_T$cancer_type == "BLCA"],
                             p_val_stage = DE_genes_Xross_stages$p_val[DE_genes_Xross_stages$cancer_type == "BLCA"])
ggplot(BLCA_p_stage_T, aes(x = -log(p_val_T), y = -log(p_val_stage), color = gene_in_seed)) + geom_point(aes(alpha = 0.2*gene_in_seed)) + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + theme_minimal()

#check the dynamics of a highly significant gene along the TimeAx trajectory:
BLCA_p_stage_T <- BLCA_p_stage_T[!is.na(BLCA_p_stage_T$gene_symbol) & !is.na(BLCA_p_stage_T$p_val_T),]
head(BLCA_p_stage_T[order(BLCA_p_stage_T$p_val_stage),], 10)

