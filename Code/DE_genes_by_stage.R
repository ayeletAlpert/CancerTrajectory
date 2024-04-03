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
DE_genes_Xross_stages <- do.call('rbind', lapply(cancer_types_short, function(cancer_type){
  message(cancer_type)
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
  
  p_vals_genes_stages <- do.call('rbind', lapply(1:nrow(cancer_type_GE_data), function(i){
    message(i)
    if(sum(assays(cancer_type_GE_data)$tpm_unstrand[i, rel_samples] == 0) == length(rel_samples)){return(NULL)}
    lm_res <- summary(lm(assays(cancer_type_GE_data)$tpm_unstrand[i, rel_samples] ~ numerical_stages))
    # plot(cancer_type_RnaseqSE_stage, assays(exp_data)$tpm_unstrand[gene,rel_samples_healthy])
    return(data.frame(gene_ensembl = rownames(cancer_type_GE_data)[i],
                      #gene_symbol = ensembl2sym(rownames(cancer_type_GE_data)[i]),
                      p_val = coef(lm_res)[2,4]))
  }))
  rel_stages <- colData(cancer_type_GE_data)[rel_samples,"ajcc_pathologic_stage"]
  rel_stages <- rel_stages[!is.na(rel_stages)]
  rel_stages <- gsub("([ABC])$", "", rel_stages)
  return(data.frame(cancer_type = sub("\\_TCGA\\.rds$", "", cancer_type), rel_stages = rel_stages))
}))