# initialize settings
rm(list = ls())

#load the required packaes:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/initial_settings.R")

data_dir_bulk_tissue <- "/media/chronos/Storage/ayelet/TCGA_bulk_GE"

# count the number of NAT samples for each cancer type --------------------------------------------------------
cancer_types <- list.dirs(data_dir_bulk_tissue, full.names = FALSE, recursive = FALSE)
normal_tissue_Xross_cancers <- do.call('rbind', lapply(cancer_types, function(cancer_type){
  print(cancer_type)
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, cancer_type, paste0(sub("^TCGA-", "", cancer_type), "_TCGA.rds")))
  # if(sum(exp_data$definition == "Solid Tissue Normal") == 0){return()}
  return(data.frame(cancer_type = sub("^TCGA-", "", cancer_type), 
                    NAT_num = length(unique(exp_data$patient[exp_data$definition == "Solid Tissue Normal"]))))
}))
normal_tissue_Xross_cancers <- normal_tissue_Xross_cancers[order(normal_tissue_Xross_cancers$NAT_num),]
normal_tissue_Xross_cancers$cancer_type <- factor(normal_tissue_Xross_cancers$cancer_type, levels = normal_tissue_Xross_cancers$cancer_type)
ggplot(normal_tissue_Xross_cancers, aes(x = cancer_type, y = NAT_num)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  geom_hline(yintercept = 30)
THRESH_NAT_samples <- 30
rel_cancers_NAT <- as.character(normal_tissue_Xross_cancers$cancer_type[normal_tissue_Xross_cancers$NAT_num >= THRESH_NAT_samples])

# NAT DEG VS stage --------------------------------------------------------

#save the NAT samples per cancer:
data_dir <- "D:/Ayelet/bulk_TCGA"
TCGA_data_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = F)[grep("/TCGA-", list.dirs(data_dir, full.names = TRUE, recursive = F))]
for(TCGA_data_dir in TCGA_data_dirs){
  cancer_type <- gsub("TCGA-", "", tools::file_path_sans_ext(basename(TCGA_data_dir)))
  print(cancer_type)
  TCGA_full_data <- readRDS(file = file.path(TCGA_data_dir, paste0(cancer_type, "_TCGA.rds")))
  if(sum(TCGA_full_data$definition == "Solid Tissue Normal") > 0){
    saveRDS(TCGA_full_data[,TCGA_full_data$definition == "Solid Tissue Normal"], 
            file = file.path(TCGA_data_dir, paste0(cancer_type, "_NAT_TCGA.rds")))
  }
}

#  FUNCTIONS FOR MAPPING GENE SYMBOLS TO EMSEBLE IDS ------------------------
exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, cancer_types[1], paste0(sub("^TCGA-", "", cancer_types[1]), "_TCGA.rds")))

res_symbol <- mapIds(org.Hs.eg.db, keys <- gsub("\\.\\d+", "", rownames(rowData(exp_data))), 
                     column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
ensembl2sym <- function(ensembl){return(res_symbol[ensembl])}
res_op_sym <- names(res_symbol)
names(res_op_sym) <- res_symbol
sym2ensembl <- function(sym){return(res_op_sym[sym])}

#  FUNCTIONS FOR MAPPING ENTREZ_ID TO EMSEBLE IDS ------------------------
res_entrez <- mapIds(org.Hs.eg.db, keys <- gsub("\\.\\d+", "", rownames(rowData(exp_data))), 
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
ensembl2entrez <- function(ensembl){return(res_entrez[ensembl])}
res_op_entrez <- names(res_entrez)
names(res_op_entrez) <- res_entrez
entrez2ensembl <- function(sym){return(res_op_entrez[sym])}

lapply(rel_cancers_NAT, function(cancer_type){
  print(cancer_type)
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_TCGA.rds")))
  
  stage_col_name <- colnames(colData(exp_data))[str_detect(colnames(colData(exp_data)), "ajcc_pathologic_stage|figo_stage|primary_gleason_grade")]
  exp_data$stage_short <- str_extract(colData(exp_data)[,stage_col_name], "[I]+V*|[0-9]+")
  
  relevant_stages <- c("I", "II", "III", "IV")
  if(cancer_type == "PRAD"){relevant_stages <- c("2", "3", "4", "5")}
  
  rel_samples_healthy <- colnames(exp_data)[exp_data$stage_short %in% relevant_stages &
                                              (exp_data$definition == "Solid Tissue Normal")]
  cancer_type_RnaseqSE_stage <- sapply(rel_samples_healthy, function(sample){
    return(which(relevant_stages == colData(exp_data)[sample,'stage_short']))})
  rownames(exp_data) <- str_replace(rownames(exp_data), "\\.[0-9]+","")

  total_genes <- nrow(exp_data)
  interval <- 1000
  cat("Progress:", "0%")
  
  p_per_gene_NAT_stage <- do.call('rbind', lapply(1:nrow(exp_data), function(gene){
    if (gene %% interval == 0) {
      percentage <- round((gene / total_genes) * 100, 2)
      cat("\rProgress:", percentage, "%")
    }
    
    lm_res <- summary(lm(assays(exp_data)$tpm_unstrand[gene,rel_samples_healthy] ~ cancer_type_RnaseqSE_stage))
    # plot(cancer_type_RnaseqSE_stage, assays(exp_data)$tpm_unstrand[gene,rel_samples_healthy])
    return(data.frame(gene_ensembl = rownames(exp_data)[gene], t_res = coefficients(lm_res)[2,3],
                      p_res = coefficients(lm_res)[2,4], coef = coefficients(lm_res)[2,1]))
  }))
  p_per_gene_NAT_stage$gene_symbol <- ensembl2sym(p_per_gene_NAT_stage$gene_ensembl)
  p_per_gene_NAT_stage$gene_entrez <- ensembl2entrez(p_per_gene_NAT_stage$gene_ensembl)

  saveRDS(p_per_gene_NAT_stage, file = file.path("D:\\Ayelet\\saved_data\\NAT_GE_by_stage", paste0(cancer_type, "__per_gene_NAT_stage.rds")))
  
  return(NULL)
})

#functional enrichment analysis across cancers:
data_dir_out_NAT_analysis <- "D:\\Ayelet\\saved_data\\NAT_GE_by_stage"
cancer_types_out_NAT_analysis <- list.files(data_dir_out_NAT_analysis)
par(mfrow = c(3,4))
skewness_p <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  p_vals_gene_NAT_stage <- readRDS(file = file.path(data_dir_out_NAT_analysis, cancer_type_out_NAT_analysis))
  p_vals_gene_NAT_stage <- p_vals_gene_NAT_stage[!is.na(p_vals_gene_NAT_stage$gene_symbol),]
  skewness_value <- skewness(p_vals_gene_NAT_stage$p_res, na.rm = T)
  hist(p_vals_gene_NAT_stage$p_res, 50, main = paste(cancer_type, "skewness", sprintf("%.2f", skewness_value)), xlab = "p values")
  
  return(data.frame(cancer_type = cancer_type, skewness = skewness_value))
}))

#for those cancers with high skewness, calculate the functional enrichment analysis of the significanr genes:
high_skewness_cancers <- skewness_p$cancer_type[skewness_p$skewness > 0.1]
enrich_pwys_per_cancer <- lapply(high_skewness_cancers, function(high_skewness_cancer){
  print(high_skewness_cancer)
  p_vals_gene_NAT_stage <- readRDS(file = file.path(data_dir_out_NAT_analysis, paste0(high_skewness_cancer,"__per_gene_NAT_stage.rds")))
  p_vals_gene_NAT_stage <- p_vals_gene_NAT_stage[!is.na(p_vals_gene_NAT_stage$gene_symbol) & !is.na(p_vals_gene_NAT_stage$t_res),]
  p_vals_gene_NAT_stage$p_adj <- p.adjust(p_vals_gene_NAT_stage$p_res, method = "BH")
  
  P_THRESH <- 0.5
  sig_genes <- p_vals_gene_NAT_stage$gene_entrez[p_vals_gene_NAT_stage$p_adj < P_THRESH]
  sig_genes <- sig_genes[!is.na(sig_genes)]

  if(length(sig_genes) == 0){return(NULL)}
  #over representation analysis:
  enrich_res_stage <- enrichKEGG(sig_genes, organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                 universe = p_vals_gene_NAT_stage$gene_entrez,
                                 minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  barplot(enrich_res_stage, showCategory=40, title = high_skewness_cancer)
  return(enrich_res_stage@result[enrich_res_stage@result$p.adjust < 0.1, c("ID", "Description", "p.adjust")])
})
names(enrich_pwys_per_cancer) <- high_skewness_cancers

#number of NAT samples per stage:
stage_distribution_NAT_per_cancer <- do.call('rbind', lapply(rel_cancers_NAT, function(rel_cancer_NAT){
  print(rel_cancer_NAT)
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",rel_cancer_NAT), 
                                       paste0(rel_cancer_NAT, "_TCGA.rds")))
  
  stage_col_name <- colnames(colData(exp_data))[str_detect(colnames(colData(exp_data)), "ajcc_pathologic_stage|figo_stage|primary_gleason_grade")]
  exp_data$stage_short <- str_extract(colData(exp_data)[,stage_col_name], "[I]+V*|[0-9]+")
  
  relevant_stages <- c("I", "II", "III", "IV")
  if(rel_cancer_NAT == "PRAD"){relevant_stages <- c("2", "3", "4", "5")}
  
  rel_samples_healthy <- colnames(exp_data)[exp_data$stage_short %in% relevant_stages &
                                              (exp_data$definition == "Solid Tissue Normal")]
  cancer_type_RnaseqSE_stage <- sapply(rel_samples_healthy, function(sample){
    return(which(relevant_stages == colData(exp_data)[sample,'stage_short']))})
  return(data.frame(cancer_type = rel_cancer_NAT, stage = table(cancer_type_RnaseqSE_stage)))
}))
colnames(stage_distribution_NAT_per_cancer) <- c("cancer_type", "stage", "freq")
stage_distribution_NAT_per_cancer_wide <- dcast(stage_distribution_NAT_per_cancer, formula = cancer_type ~ stage, value.var = "freq")
rownames(stage_distribution_NAT_per_cancer_wide) <- stage_distribution_NAT_per_cancer_wide$cancer_type
stage_distribution_NAT_per_cancer_wide <- stage_distribution_NAT_per_cancer_wide[,-1]
pheatmap(stage_distribution_NAT_per_cancer_wide, cluster_cols = F, display_numbers = T)

#get a list of the hallmark gene set:
hallmark_genes_path <- "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt"
hallmark_genes <- read.gmt(hallmark_genes_path)

p_val_hallmark_pathway <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_TCGA.rds")))
  
  stage_col_name <- colnames(colData(exp_data))[str_detect(colnames(colData(exp_data)), "ajcc_pathologic_stage|figo_stage|primary_gleason_grade")]
  exp_data$stage_short <- str_extract(colData(exp_data)[,stage_col_name], "[I]+V*|[0-9]+")
  
  relevant_stages <- c("I", "II", "III", "IV")
  if(cancer_type == "PRAD"){relevant_stages <- c("2", "3", "4", "5")}
  
  rel_samples_healthy <- colnames(exp_data)[exp_data$stage_short %in% relevant_stages &
                                              (exp_data$definition == "Solid Tissue Normal")]
  cancer_type_RnaseqSE_stage <- sapply(rel_samples_healthy, function(sample){
    return(which(relevant_stages == colData(exp_data)[sample,'stage_short']))})

  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data_scaled <- t(apply(assays(exp_data)$tpm_unstrand, 1, scale))
  colnames(exp_data_scaled) <- colnames(exp_data)
  
  do.call('rbind', lapply(as.character(unique(hallmark_genes$term)), function(hallmark_gene_set){
    print(hallmark_gene_set)
    genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
    
    exp_data_scaled_pathway <- exp_data_scaled[intersect(genes_in_pathway, rownames(exp_data_scaled)),]
    exp_data_scaled_pathway <- exp_data_scaled_pathway[complete.cases(exp_data_scaled_pathway), ]
    
    exp_data_scaled_pathway_mean <- apply(exp_data_scaled_pathway,2,mean)
    df_stage_mean_exp <- data.frame(mean_exp_gene_set = exp_data_scaled_pathway_mean[rel_samples_healthy], 
                                    stage = cancer_type_RnaseqSE_stage,
                                    age = colData(exp_data)[rel_samples_healthy,'days_to_birth'])
    df_stage_mean_exp$stage_met <- df_stage_mean_exp$stage == 4
    ggplot(df_stage_mean_exp, aes(x = factor(stage_met), y = mean_exp_gene_set)) + geom_boxplot()
    ggplot(df_stage_mean_exp, aes(x = factor(stage), y = mean_exp_gene_set)) + geom_boxplot() + geom_jitter()
    ggplot(df_stage_mean_exp, aes(x = age, y = mean_exp_gene_set)) + geom_point() 
    summary(lm(mean_exp_gene_set ~ age + stage, data = df_stage_mean_exp))
    if(sum(df_stage_mean_exp$stage_met) < 3){return(NULL)}
    tt_res <- t.test(df_stage_mean_exp$mean_exp_gene_set[df_stage_mean_exp$stage_met],
                     df_stage_mean_exp$mean_exp_gene_set[!df_stage_mean_exp$stage_met])
    
    # lm_res <- summary(lm(exp_data_scaled_pathway_mean[rel_samples_healthy] ~ cancer_type_RnaseqSE_stage))
    # plot(cancer_type_RnaseqSE_stage, exp_data_scaled_pathway_mean[rel_samples_healthy])
    
    return(data.frame(cancer_type = cancer_type, pathway = hallmark_gene_set, p_val = tt_res[[3]]))
  }))
}))
unique(p_val_hallmark_pathway$cancer_type)
hist(p_val_hallmark_pathway$p_val, 20)
p_val_hallmark_pathway[p_val_hallmark_pathway$p_val < 0.05,]
ggplot(p_val_hallmark_pathway, aes(x = cancer_type, y = p_val)) + geom_boxplot()

exp_data_across_pathways <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_TCGA.rds")))
  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data_tpm <- assays(exp_data)$unstranded
  return(data.frame(t(exp_data_tpm[intersect(inflammatory_genes, rownames(exp_data_tpm)),]), cancer_type = cancer_type))
}))
exp_data_across_pathways <- exp_data_across_pathways[exp_data_across_pathways$cancer_type %in% c("BRCA", "COAD", "LUAD", "PRAD"),]
exp_data_across_pathways <- data.frame(log(exp_data_across_pathways[,1:(ncol(exp_data_across_pathways) - 1)]), cancer_type = exp_data_across_pathways$cancer_type)
exp_data_scaled <- t(data.frame(apply(exp_data_across_pathways[,1:(ncol(exp_data_across_pathways) - 1)], 2, scale)))
exp_data_scaled <- exp_data_scaled[complete.cases(exp_data_scaled),]
colnames(exp_data_scaled) <- rownames(exp_data_across_pathways)

gene_exp_pathways <- do.call('rbind', lapply(c("BRCA", "COAD", "LUAD", "PRAD"), function(cancer_type){
  do.call('rbind', lapply(as.character(inflammatory_pathways), function(hallmark_gene_set){
    print(hallmark_gene_set)
    genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
    
    res <- apply(exp_data_scaled[intersect(rownames(exp_data_scaled), genes_in_pathway), exp_data_across_pathways$cancer_type == cancer_type], 1, mean)
    
  return(data.frame(cancer_type = cancer_type, hallmark_gene_set = hallmark_gene_set, res = res))
}))}))
ggplot(gene_exp_pathways, aes(x = cancer_type, y = res)) + geom_boxplot() + facet_wrap(~hallmark_gene_set)
mean_exp <- apply(exp_data_scaled, 1, mean)

do.call('rbind', lapply(as.character(unique(hallmark_genes$term)), function(hallmark_gene_set){
  print(hallmark_gene_set)
  genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
  
  exp_data_scaled_pathway <- exp_data_scaled[intersect(genes_in_pathway, rownames(exp_data_scaled)),]
  exp_data_scaled_pathway <- exp_data_scaled_pathway[complete.cases(exp_data_scaled_pathway), ]
  
  exp_data_scaled_pathway_mean <- apply(exp_data_scaled_pathway,2,mean)
  
  return(data.frame(cancer_type = cancer_type, pathway = hallmark_gene_set, p_val = tt_res[[3]]))
}))


# Generate a trajectory per inflammatory pathway --------------------------
#get a list of the hallmark gene set:
hallmark_genes_path <- "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt"
hallmark_genes <- read.gmt(hallmark_genes_path)

#calcuate the intra-pathway correlation between genes in the NAT samples across cancers:
cor_hallmark_pathway_vs_rand <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_NAT_TCGA.rds")))
  
  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data <- exp_data[!is.na(rownames(exp_data)),]
  
  exp_data_tpm <- assays(exp_data)$tpm_unstrand
  
  hallmark_median_cor = do.call('rbind', lapply(as.character(unique(hallmark_genes$term)), function(hallmark_gene_set){
    print(hallmark_gene_set)
    genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
    
    exp_data_pathway <- exp_data_tpm[intersect(genes_in_pathway, rownames(exp_data_tpm)),]
    cor_matrix_pathway <- abs(cor(t(exp_data_pathway)))
    
    rand_genes <- sample(rownames(exp_data_tpm), length(genes_in_pathway))
    
    exp_data_rand <- exp_data_tpm[rand_genes,]
    cor_matrix_rand <- abs(cor(t(exp_data_rand)))
    
    return(data.frame(cancer_type = cancer_type, hallmark_pathway = hallmark_gene_set, 
                      median_cor_hallmark = median(cor_matrix_pathway[upper.tri(cor_matrix_pathway, diag = F)], na.rm = T),
                      median_cor_rand = median(cor_matrix_rand[upper.tri(cor_matrix_rand, diag = F)], na.rm = T)))
  }))
}))
cor_hallmark_pathway_vs_rand_long <- melt(cor_hallmark_pathway_vs_rand, id.vars = c("cancer_type", "hallmark_pathway"))
ggplot(cor_hallmark_pathway_vs_rand_long, aes(x = variable, y = value)) + geom_boxplot() + facet_wrap(~cancer_type)
ggplot(cor_hallmark_pathway_vs_rand_long, aes(x = variable, y = value)) + geom_boxplot() + facet_wrap(~hallmark_pathway)

#calculate a trajectory per signaling pathway using the NAT:
inflammatory_pathways <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                           "HALLMARK_IL2_STAT5_SIGNALING")

# Step 1: Create a Binary Matrix
binary_matrix <- table(hallmark_genes$term, hallmark_genes$gene)

# Step 2: Compute Pairwise Overlap (Jaccard Similarity)
terms <- rownames(binary_matrix)
overlap_matrix <- matrix(NA, nrow = length(terms), ncol = length(terms), dimnames = list(terms, terms))

for (i in 1:length(terms)) {
  for (j in 1:length(terms)) {
    if (i != j) {
      genes_term1 <- binary_matrix[terms[i], ]
      genes_term2 <- binary_matrix[terms[j], ]
      overlap_score <- sum(genes_term1 & genes_term2) / sum(genes_term1 | genes_term2)
      overlap_matrix[terms[i], terms[j]] <- overlap_score
    }
  }
}
overlap_matrix <- overlap_matrix[inflammatory_pathways,inflammatory_pathways]
pheatmap(overlap_matrix)

data_dir_out_NAT_analysis <- "D:\\Ayelet\\saved_data\\NAT_GE_by_stage"
cancer_types_out_NAT_analysis <- list.files(data_dir_out_NAT_analysis)
hallmark_gene_set <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
hallmark_gene_set <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
hallmark_gene_set <- "HALLMARK_IL2_STAT5_SIGNALING"
pathway_across_cancers <- do.call('rbind', lapply(inflammatory_pathways, function(hallmark_gene_set){
  print(hallmark_gene_set)
  genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
  
  pathway_across_cancers <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
    cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
    print(cancer_type)
    
    exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                         paste0(cancer_type, "_NAT_TCGA.rds")))
    
    rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
    exp_data <- exp_data[!is.na(rownames(exp_data)),]
    
    exp_data_tpm <- assays(exp_data)$tpm_unstrand
    
    exp_data_pathway <- exp_data_tpm[intersect(genes_in_pathway, rownames(exp_data_tpm)),]
    
    # pca_res <- data.frame(prcomp(t(exp_data_pathway), scale. = T, center = T)$x[,c("PC1", "PC2")])
    # biplot(prcomp(t(exp_data_pathway), scale. = T, center = T))
    # 
    # exp_data_pathway_scaled <- t(apply(exp_data_pathway,1,scale))
    # colnames(exp_data_pathway_scaled) <- colnames(exp_data_pathway)
    # exp_data_pathway_scaled[exp_data_pathway_scaled > 2] <- 2
    # exp_data_pathway_scaled[exp_data_pathway_scaled < -2] <- -2
    # 
    # pca_res <- cbind(pca_res, data.frame(age_diagnosis = exp_data$age_at_diagnosis))
    # pca_res <- cbind(pca_res, t(exp_data_pathway_scaled))
    # pheatmap(pca_res[order(pca_res$PC1),rownames(exp_data_pathway_scaled)], cluster_rows = F)
    # 
    # clinical_data <- readRDS(file = "D:/Ayelet/bulk_TCGA/TCGA-BRCA/BRCA_clinical_data.rds")
    # 
    # # Extract relevant columns for hormone receptor status
    # hormone_receptor_status <- clinical_data[, c("bcr_patient_barcode", 
    #                                              "breast_carcinoma_estrogen_receptor_status", 
    #                                              "breast_carcinoma_progesterone_receptor_status",
    #                                              "her2_immunohistochemistry_level_result",
    #                                              "her2_erbb_pos_finding_cell_percent_category",
    #                                              "her2_neu_chromosone_17_signal_ratio_value", 
    #                                              "lab_procedure_her2_neu_in_situ_hybrid_outcome_type")]
    # hormone_receptor_status_NAT <- hormone_receptor_status[hormone_receptor_status$bcr_patient_barcode %in% sub("^([^\\-]+\\-[^\\-]+\\-[^\\-]+).*", "\\1", colnames(exp_data)),]
    # hormone_receptor_status_NAT <- unique(hormone_receptor_status_NAT)
    # rownames(hormone_receptor_status_NAT) <- hormone_receptor_status_NAT$bcr_patient_barcode
    # 
    # rownames(pca_res) <- sub("^([^\\-]+\\-[^\\-]+\\-[^\\-]+).*", "\\1", rownames(pca_res))
    # pca_res <- cbind(pca_res, hormone_receptor_status_NAT[rownames(pca_res),])
    # 
    # ggplot(pca_res, aes(x = PC1, y = PC2, color = breast_carcinoma_estrogen_receptor_status)) + geom_point()
    # ggplot(pca_res, aes(x = breast_carcinoma_estrogen_receptor_status, y = PC1)) + geom_boxplot()
    # t.test(pca_res$PC1[pca_res$breast_carcinoma_estrogen_receptor_status == "Positive"],
    #        pca_res$PC1[pca_res$breast_carcinoma_estrogen_receptor_status == "Negative"])
    # 
    # ggplot(pca_res, aes(x = PC1, y = PC2, color = her2_immunohistochemistry_level_result)) + geom_point()
    # ggplot(pca_res, aes(x = her2_immunohistochemistry_level_result, y = PC1)) + geom_boxplot()
    # t.test(pca_res$PC1[pca_res$her2_immunohistochemistry_level_result == "Positive"],
    #        pca_res$PC1[pca_res$her2_immunohistochemistry_level_result == "Negative"])
    # 
    # 
    # ggplot(pca_res, aes(x = PC1, y = PC2, color = breast_carcinoma_progesterone_receptor_status)) + geom_point()
    # ggplot(pca_res, aes(x = breast_carcinoma_progesterone_receptor_status, y = PC1)) + geom_boxplot()
    # t.test(pca_res$PC1[pca_res$breast_carcinoma_progesterone_receptor_status == "Positive"],
    #        pca_res$PC1[pca_res$breast_carcinoma_progesterone_receptor_status == "Negative"])
    # 
    # sum(sub("^([^\\-]+\\-[^\\-]+\\-[^\\-]+).*", "\\1", colnames(exp_data)) %in% hormone_receptor_status$bcr_patient_barcode)
    # hormone_receptor_status[,c("breast_carcinoma_estrogen_receptor_status", "breast_carcinoma_progesterone_receptor_status")]
    # ggplot(pca_res, aes(x = PC1, y = PC2, color = log(TRIM5))) + geom_point()
    # ggplot(pca_res, aes(x = PC1, y = PC2, color = age_diagnosis)) + geom_point()
    # ggplot(pca_res, aes(x = age_diagnosis, y = PC1)) + geom_point()
    # ggplot(pca_res, aes(x = age_diagnosis, y = PC2)) + geom_point()
    # 
    # summary(lm(PC1 ~ age_diagnosis, data = pca_res))
    # summary(lm(PC2 ~ age_diagnosis, data = pca_res))
    #

    # pca_res <- prcomp(pathway_across_cancers[,1:(ncol(pathway_across_cancers) - 1)], scale. = T)$x[,c("PC1", "PC2")]
    # pca_res <- data.frame(pca_res, cancer_type = pathway_across_cancers$cancer_type)
    # biplot(prcomp(t(exp_data_pathway), scale. = T, center = T))
    # ggplot(pca_res, aes(x = cancer_type, y = PC1)) + geom_boxplot()
    return(data.frame(t(exp_data_pathway), cancer_type = cancer_type))
  }))
  
  pca_res <- data.frame(prcomp(pathway_across_cancers[,intersect(genes_in_pathway, colnames(pathway_across_cancers))],
                               scale. = T, center = T)$x[,c("PC1", "PC2", "PC3")],
                        cancer_type = pathway_across_cancers$cancer_type,
                        hallmark_gene_set = hallmark_gene_set)
  # biplot(prcomp(pathway_across_cancers[,intersect(genes_in_pathway, colnames(pathway_across_cancers))], scale. = T, center = T))
  # ggplot(pca_res, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point() + ggtitle(hallmark_gene_set) + theme_classic()
  # ggplot(pca_res, aes(x = cancer_type, y = PC1)) + geom_boxplot() + ggtitle(hallmark_gene_set) + theme_classic()
  return(data.frame(pca_res, patient_id = rownames(pca_res)))
}))
ggplot(traj_per_cancer, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point() + ggtitle(hallmark_gene_set) + 
  theme_classic() + facet_wrap(~hallmark_gene_set, scales = "free")
dcast(pathway_across_cancers, patient_id ~ hallmark_gene_set + cancer_type, value.var = 'PC1')
plot(pathway_across_cancers$PC1[pathway_across_cancers$hallmark_gene_set == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"],
     pathway_across_cancers$PC1[pathway_across_cancers$hallmark_gene_set == "HALLMARK_IL2_STAT5_SIGNALING"])


pathway_across_cancers <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type),
                                       paste0(cancer_type, "_NAT_TCGA.rds")))
  
  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data <- exp_data[!is.na(rownames(exp_data)),]
  
  # exp_data_tpm <- t(assays(exp_data)$tpm_unstrand)
  exp_data_tpm <- t(assays(exp_data)$unstrand)
  
  return(data.frame(exp_data_tpm[,c("ACTB", "GAPDH")], cancer_type = cancer_type))
}))
ggplot(pathway_across_cancers, aes(x = cancer_type, y = ACTB)) + geom_boxplot()
ggplot(pathway_across_cancers, aes(x = log(GAPDH), y = log(ACTB), color = cancer_type)) + geom_point() + geom_smooth(method = 'lm', se = F)

cancer_type_out_NAT_analysis <- "BRCA__per_gene_NAT_stage.rds"
clinical_data <- readRDS(file = "D:/Ayelet/bulk_TCGA/TCGA-BRCA/BRCA_clinical_data.rds")

# Extract relevant columns for hormone receptor status
hormone_receptor_status <- clinical_data[, c("bcr_patient_barcode", 
                                             "breast_carcinoma_estrogen_receptor_status", 
                                             "breast_carcinoma_progesterone_receptor_status",
                                             "her2_immunohistochemistry_level_result",
                                             "her2_erbb_pos_finding_cell_percent_category",
                                             "her2_neu_chromosone_17_signal_ratio_value", 
                                             "lab_procedure_her2_neu_in_situ_hybrid_outcome_type", 
                                             "menopause_status")]


pathways_BRCA <- do.call('rbind', lapply(inflammatory_pathways, function(hallmark_gene_set){
  print(hallmark_gene_set)
  genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == hallmark_gene_set]
  
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_NAT_TCGA.rds")))
  
  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data <- exp_data[!is.na(rownames(exp_data)),]
  
  exp_data_tpm <- assays(exp_data)$tpm_unstrand
  
  exp_data_pathway <- exp_data_tpm[intersect(genes_in_pathway, rownames(exp_data_tpm)),]
  
  pca_res <- data.frame(prcomp(t(exp_data_pathway), scale. = T, center = T)$x[,c("PC1", "PC2")])
  
  exp_data_pathway_scaled <- t(apply(exp_data_pathway,1,scale))
  colnames(exp_data_pathway_scaled) <- colnames(exp_data_pathway)
  # exp_data_pathway_scaled[exp_data_pathway_scaled > 2] <- 2
  # exp_data_pathway_scaled[exp_data_pathway_scaled < -2] <- -2
  
  pca_res <- cbind(pca_res, data.frame(age_diagnosis = exp_data$age_at_diagnosis))
  pca_res <- cbind(pca_res, t(exp_data_pathway_scaled))
  rownames(pca_res) <- sub("^([^\\-]+\\-[^\\-]+\\-[^\\-]+).*", "\\1", rownames(pca_res))
  # pheatmap(pca_res[order(pca_res$PC1),rownames(exp_data_pathway_scaled)], cluster_rows = F)
  
  hormone_receptor_status_NAT <- hormone_receptor_status[hormone_receptor_status$bcr_patient_barcode %in% sub("^([^\\-]+\\-[^\\-]+\\-[^\\-]+).*", "\\1", colnames(exp_data)),]
  hormone_receptor_status_NAT <- unique(hormone_receptor_status_NAT)
  rownames(hormone_receptor_status_NAT) <- hormone_receptor_status_NAT$bcr_patient_barcode
    
  pca_res <- cbind(pca_res, hormone_receptor_status_NAT[rownames(pca_res),])

  # ggplot(pca_res, aes(x = PC1, y = PC2, color = breast_carcinoma_estrogen_receptor_status)) + geom_point()
  # ggplot(pca_res, aes(x = breast_carcinoma_estrogen_receptor_status, y = PC1)) + geom_boxplot()
  tt_res_PC1 =  t.test(pca_res$PC1[pca_res$breast_carcinoma_estrogen_receptor_status == "Positive"],
                       pca_res$PC1[pca_res$breast_carcinoma_estrogen_receptor_status == "Negative"])
  tt_res_PC2 =  t.test(pca_res$PC2[pca_res$breast_carcinoma_estrogen_receptor_status == "Positive"],
                       pca_res$PC2[pca_res$breast_carcinoma_estrogen_receptor_status == "Negative"])
  
  lmres_PC1_age = summary(lm(PC1 ~ age_diagnosis, data = pca_res))
  lmres_PC2_age = summary(lm(PC2 ~ age_diagnosis, data = pca_res))
  
  return(data.frame(hallmark_gene_set = hallmark_gene_set,
                    p_PC1_estrogen = tt_res_PC1[[3]], 
                    p_PC2_estrogen = tt_res_PC2[[3]], 
                    coef_PC1_age = coefficients(lmres_PC1_age)[2,4],
                    coef_PC2_age = coefficients(lmres_PC2_age)[2,4]))
    
  
  # ggplot(pca_res, aes(x = PC1, y = PC2, color = her2_immunohistochemistry_level_result)) + geom_point()
  # ggplot(pca_res, aes(x = her2_immunohistochemistry_level_result, y = PC1)) + geom_boxplot()
  # t.test(pca_res$PC1[pca_res$her2_immunohistochemistry_level_result == "Positive"],
  #        pca_res$PC1[pca_res$her2_immunohistochemistry_level_result == "Negative"])
  # 
  # 
  # ggplot(pca_res, aes(x = PC1, y = PC2, color = breast_carcinoma_progesterone_receptor_status)) + geom_point()
  # ggplot(pca_res, aes(x = breast_carcinoma_progesterone_receptor_status, y = PC1)) + geom_boxplot()
  # t.test(pca_res$PC1[pca_res$breast_carcinoma_progesterone_receptor_status == "Positive"],
  #        pca_res$PC1[pca_res$breast_carcinoma_progesterone_receptor_status == "Negative"])
  # 
  #   hormone_receptor_status[,c("breast_carcinoma_estrogen_receptor_status", "breast_carcinoma_progesterone_receptor_status")]
  #   ggplot(pca_res, aes(x = PC1, y = PC2, color = log(TRIM5))) + geom_point()
  #   ggplot(pca_res, aes(x = PC1, y = PC2, color = age_diagnosis)) + geom_point()
  #   ggplot(pca_res, aes(x = age_diagnosis, y = PC1)) + geom_point()
  #   ggplot(pca_res, aes(x = age_diagnosis, y = PC2)) + geom_point()
  #   
  #   summary(lm(PC1 ~ age_diagnosis, data = pca_res))
  #   summary(lm(PC2 ~ age_diagnosis, data = pca_res))
  #   
  #   return(data.frame(t(exp_data_pathway), cancer_type = cancer_type))
  # }))
  
  # pca_res <- data.frame(prcomp(pathway_across_cancers[,intersect(genes_in_pathway, colnames(pathway_across_cancers))], 
  #                              scale. = T, center = T)$x[,c("PC1", "PC2")], 
  #                       cancer_type = pathway_across_cancers$cancer_type,
  #                       hallmark_gene_set = hallmark_gene_set)
  # ggplot(pca_res, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point() + ggtitle(hallmark_gene_set) + theme_classic()
  # return(data.frame(pathway_across_cancers, hallmark_gene_set = hallmark_gene_set))
}))


# Kidney cancer - KIRC, KIRP ----------------------------------------------
hallmark_genes_path <- "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt"
hallmark_genes <- read.gmt(hallmark_genes_path)
data_dir_out_NAT_analysis <- "D:\\Ayelet\\saved_data\\NAT_GE_by_stage"
inflammatory_pathways <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                           "HALLMARK_IL2_STAT5_SIGNALING")
inflammatory_genes <- unique(hallmark_genes$gene)

cancer_types <- c("KIRC", "KIRP")
exp_data_kidney <- do.call('rbind', lapply(cancer_types, function(cancer_type){
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_NAT_TCGA.rds")))
  rownames(exp_data) <- ensembl2sym(str_replace(rownames(exp_data), "\\.[0-9]+",""))
  exp_data <- exp_data[!is.na(rownames(exp_data)),]
  exp_data <- exp_data[intersect(inflammatory_genes, rownames(exp_data)),]
  exp_data_tpm <- data.frame(t(assays(exp_data)$tpm_unstrand), check.names = F)
  
  exp_data_tpm$cancer_type <- cancer_type
  exp_data_tpm
}))
remove_col <- which(apply(exp_data_kidney[,1:(ncol(exp_data_kidney) - 1)],2,var) == 0)
exp_data_kidney <- exp_data_kidney[,-remove_col]

inflammatory_pathway_PCA <- do.call('rbind', lapply(inflammatory_pathways, function(inflammatory_pathway){
  print(inflammatory_pathway)
  genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == inflammatory_pathway]
  exp_data_kidney_pathway <- exp_data_kidney[,intersect(genes_in_pathway, colnames(exp_data_kidney))]
  pca_res_kidney_pathway <- data.frame(prcomp(exp_data_kidney_pathway, center = T, scale. = T)$x[,c("PC1", "PC2")])
  pca_res_kidney_pathway$cancer_type <- exp_data_kidney$cancer_type
  
  ggplot(pca_res_kidney_pathway, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point() + ggtitle(inflammatory_pathway)
  ggplot(pca_res_kidney, aes(x = cancer_type, y = PC1)) + geom_boxplot() + ggtitle(inflammatory_pathway)
  
  ttest_resPC1 <- t.test(pca_res_kidney_pathway$PC1[pca_res_kidney_pathway$cancer_type == "KIRC"],
                      pca_res_kidney_pathway$PC1[pca_res_kidney_pathway$cancer_type == "KIRP"])[[3]]
  ttest_resPC2 <- t.test(pca_res_kidney_pathway$PC2[pca_res_kidney_pathway$cancer_type == "KIRC"],
                         pca_res_kidney_pathway$PC2[pca_res_kidney_pathway$cancer_type == "KIRP"])[[3]]
  return(data.frame(inflammatory_pathway = inflammatory_pathway, p_val_tt_PC1 = ttest_resPC1, p_val_tt_PC2 = ttest_resPC2))
}))

pca_res_kidney <- data.frame(prcomp(exp_data_kidney[,1:(ncol(exp_data_kidney) - 1)], center = T, scale. = T)$x[,c("PC1", "PC2")])
pca_res_kidney$cancer_type <- exp_data_kidney$cancer_type
ggplot(pca_res_kidney, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point()
ggplot(pca_res_kidney, aes(x = cancer_type, y = PC1)) + geom_boxplot()
t.test(pca_res_kidney$PC1[pca_res_kidney$cancer_type == "KIRC"], 
       pca_res_kidney$PC1[pca_res_kidney$cancer_type == "KIRP"])
ggplot(pca_res_kidney, aes(x = cancer_type, y = PC2)) + geom_boxplot()


# check age distribution across cancers -----------------------------------
data_dir_out_NAT_analysis <- "D:\\Ayelet\\saved_data\\NAT_GE_by_stage"
cancer_types_out_NAT_analysis <- list.files(data_dir_out_NAT_analysis)

age_dist_across_cancers <- do.call('rbind', lapply(cancer_types_out_NAT_analysis, function(cancer_type_out_NAT_analysis){
  cancer_type <- sub("__.*", "", cancer_type_out_NAT_analysis)
  print(cancer_type)
  
  exp_data <- readRDS(file = file.path(data_dir_bulk_tissue, paste0("TCGA-",cancer_type), 
                                       paste0(cancer_type, "_NAT_TCGA.rds")))
  return(data.frame(cancer_type = cancer_type, age = exp_data$age_at_diagnosis, gender = exp_data$gender))
}))
age_dist_across_cancers$age <- age_dist_across_cancers$age/365
ggplot(age_dist_across_cancers, aes(x = cancer_type, y = age, color = cancer_type)) + geom_violin(alpha = 0.2) + geom_jitter()

gender_dist_across_Cancers <- do.call('rbind', lapply(unique(age_dist_across_cancers$cancer_type), function(cancer_type){
  return(data.frame(cancer_type = cancer_type, 
                    percent_male = sum(age_dist_across_cancers$gender == "male" & age_dist_across_cancers$cancer_type == cancer_type)/sum(age_dist_across_cancers$cancer_type == cancer_type)))
}))
ggplot(gender_dist_across_Cancers, aes(x = cancer_type, y = percent_male)) + geom_bar(stat = "identity")

detach("package:openssl", unload = TRUE)

# plasma proteomics -------------------------------------------------------
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
protein_name <- c("IFNG", "IL6", "IL2", "TNF")
ggplot(data_proteomics[data_proteomics$protein_name %in% protein_name,], aes(x = Cancer, y = NPX)) + 
  geom_boxplot() + geom_jitter(alpha = 0.2) + facet_wrap(~protein_name, scales = "free")
mean_exp_cytokine_cancer <- do.call('rbind', lapply(protein_name, function(protein){
  do.call('rbind', lapply(solid_cancers, function(cancer_type){
    mean_exp = mean(data_proteomics$NPX[data_proteomics$protein_name == protein & data_proteomics$Cancer == cancer_type])
    return(data.frame(cancer_type = cancer_type, protein = protein, mean_exp = mean_exp))
  }))
}))
t.test(data_proteomics$NPX[data_proteomics$protein_name == protein_name & data_proteomics$Cancer == "GLIOM"],
       data_proteomics$NPX[data_proteomics$protein_name == protein_name & data_proteomics$Cancer == "LUNGC"])

hallmark_genes_path <- "~/MDPhD/oncology/BayesPrismLUAD/Data/hallmark_genesets_gene_symbols.gmt"
hallmark_genes <- read.gmt(hallmark_genes_path)
data_dir_out_NAT_analysis <- "D:\\Ayelet\\saved_data\\NAT_GE_by_stage"
inflammatory_pathways <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                           "HALLMARK_IL2_STAT5_SIGNALING")
inflammatory_genes <- unique(hallmark_genes$gene[hallmark_genes$term %in% inflammatory_pathways])
sum(inflammatory_genes %in% data_proteomics$protein_name)

fraction_proteins_inflammatory_pathway <- do.call('rbind', lapply(inflammatory_pathways, function(inflammatory_pathway){
  fraction_in_prot_data = sum(hallmark_genes$gene[hallmark_genes$term == inflammatory_pathway] %in% data_proteomics$protein_name)/sum(hallmark_genes$term == inflammatory_pathway)
  return(data.frame(inflammatory_pathway = inflammatory_pathway,  fraction_in_prot_data = fraction_in_prot_data))
}))

ggplot(fraction_proteins_inflammatory_pathway, aes(x = inflammatory_pathway, y = fraction_in_prot_data)) +
  geom_point() +
  geom_vline(aes(xintercept = fraction_in_prot_data), linetype = "dashed", color = "red") + coord_flip()


# Classification of cancer classes based on inflammatory genes only:
pathway <- inflammatory_pathways[1]
genes_in_pathway <- hallmark_genes$gene[hallmark_genes$term == inflammatory_pathway]
data_proteomics_inflammation <- data_proteomics[data_proteomics$protein_name %in% inflammatory_genes, c("Sample_ID", "Cancer", "protein_name", "NPX")]
data_proteomics_inflammation_wide <- dcast(data_proteomics_inflammation, formula = Sample_ID ~ protein_name, value.var = "NPX")
data_proteomics_inflammation_wide <- data_proteomics_inflammation_wide[complete.cases(data_proteomics_inflammation_wide),]
class <- sub("_(.*)", "", data_proteomics_inflammation_wide$Sample_ID)
data_proteomics_inflammation_wide <- data_proteomics_inflammation_wide[,-1]

fit <- cv.glmnet(as.matrix(data_proteomics_inflammation_wide), class, family = "multinomial")
print(fit)

#test the model o the training data:
AUC_inflammation <- sapply(1:8, function(class_i){
  pred_prob <- predict(fit, newx = as.matrix(data_proteomics_inflammation_wide), type = "response")
  pred_prob_pos_class <- pred_prob[, class_i,]
  roc_curve <- multiclass.roc(class, pred_prob_pos_class)
  auc_value <- auc(roc_curve)
})

rand_genes <- sample(data_proteomics$protein_name, sum(inflammatory_genes %in% data_proteomics$protein_name))
data_proteomics_rand <- data_proteomics[data_proteomics$protein_name %in% rand_genes, c("Sample_ID", "Cancer", "protein_name", "NPX")]
data_proteomics_rand_wide <- dcast(data_proteomics_rand, formula = Sample_ID ~ protein_name, value.var = "NPX")
data_proteomics_rand_wide <- data_proteomics_rand_wide[complete.cases(data_proteomics_rand_wide),]
class <- sub("_(.*)", "", data_proteomics_rand_wide$Sample_ID)
data_proteomics_rand_wide <- data_proteomics_rand_wide[,-1]
fit_rand <- cv.glmnet(as.matrix(data_proteomics_rand_wide), class, family = "multinomial")
# print(fit)

#test the model o the training data:
AUC_rand <- sapply(1:8, function(class_i){
  pred_prob <- predict(fit_rand, newx = as.matrix(data_proteomics_rand_wide), type = "response")
  pred_prob_pos_class <- pred_prob[, class_i,]
  roc_curve <- multiclass.roc(class, pred_prob_pos_class)
  auc_value <- auc(roc_curve)
})
boxplot(AUC_inflammation, AUC_rand)

#association of the infammatory proteins with disease type:
names(inflammatory_pathways) <- c("IL6", "IFNG", "TNF", "TGFB1", "IFNG", "IL2")
pathway_activation <- lapply(inflammatory_pathways, function(inflammatory_pathway){
  proteins_in_pathway <- hallmark_genes$gene[hallmark_genes$term == inflammatory_pathway]
  data_proteomics_pathway <- data_proteomics[data_proteomics$protein_name %in% proteins_in_pathway, c("Sample_ID", "Cancer", "protein_name", "NPX")]
  data_proteomics_pathway_wide <- dcast(data_proteomics_pathway, formula = Sample_ID ~ protein_name, value.var = "NPX")
  data_proteomics_pathway_wide <- data_proteomics_pathway_wide[complete.cases(data_proteomics_pathway_wide),]
  pca_res_prot <- data.frame(prcomp(data_proteomics_pathway_wide[,2:ncol(data_proteomics_pathway_wide)], scale. = T)$x[,paste0("PC", 1:10)])
  rownames(pca_res_prot) <- data_proteomics_pathway_wide$Sample_ID
  pca_res_prot$cancer_type <- sub("_(.*)", "", data_proteomics_pathway_wide$Sample_ID)
  cor_PCs = sapply(1:10, function(i){
    cor(pca_res_prot[,paste0("PC", i)], data_proteomics_wide[rownames(pca_res_prot),names(inflammatory_pathways)[inflammatory_pathways == inflammatory_pathway]],
        method = "spearman")
  })
  res <- pca_res_prot[,paste0("PC", which.max(abs(cor_PCs)))]*sign(cor_PCs[which.max(abs(cor_PCs))])
  names(res) <- rownames(pca_res_prot)
  return(res)
  
  # plot(pca_res_prot$PC3, data_proteomics_wide[rownames(pca_res_prot),names(inflammatory_pathways)[inflammatory_pathways == inflammatory_pathway]])
  # 
  # ggplot(pca_res_prot, aes(x = PC1, y = PC2, color = cancer_type)) + geom_point() + ggtitle(inflammatory_pathway)
  # ggplot(pca_res_prot, aes(x = cancer_type, y = PC1, color = cancer_type)) + geom_boxplot() + ggtitle(inflammatory_pathway) + geom_jitter(alpha = 0.2)
  # ggplot(pca_res_prot, aes(x = cancer_type, y = PC2, color = cancer_type)) + geom_boxplot() + ggtitle(inflammatory_pathway) + geom_jitter(alpha = 0.2)
  # pca_res_prot$lead_cytokine <- sapply(rownames(pca_res_prot), function(sample_ID){return(data_proteomics)})
  # ggplot(pca_res_prot, aes(x = IL6, y = PC1, color = cancer_type)) + geom_point()
  # ggplot(pca_res_prot, aes(x = IL6, y = ..density.., fill = cancer_type)) + geom_density(alpha = 0.2)
  # summary(lm(IL6 ~ PC1, data = pca_res_prot))
})

# Extract the names of vectors
common_patients <- Reduce(intersect, sapply(pathway_activation, names))
pathways_activation <- do.call('rbind', lapply(pathway_activation, function(x){return(x[common_patients])}))
pheatmap(cor(t(pathways_activation)), display_numbers = T, fontsize = 14, main = "correlations between inflammatory pathways\nestimated using PCA")


# Unsupervised approach to identifying modules of cytokines ---------------
inflammatory_proteins <- unique(data_proteomics$protein_name[data_proteomics$Panel == "Inflammation"])
inflammatory_proteins <- intersect(inflammatory_proteins, colnames(data_proteomics_wide))
umap_result <- umap(t(data_proteomics_wide[,inflammatory_proteins]), n_neighbors = 30, n_epochs = 5, learning_rate = 0.5)
plot(umap_result$layout[, 1], umap_result$layout[, 2], col = rainbow(30)[1], pch = 16)

library(Rtsne)
tsne_result <- Rtsne(t(data_proteomics_wide[,3:ncol(data_proteomics_wide)]), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
tsne_res <- data.frame(tSNE1 = tsne_result$Y[,1], tSNE2 = tsne_result$Y[,2])
rownames(tsne_res) <- colnames(data_proteomics_wide[,3:ncol(data_proteomics_wide)])
tsne_res$class <- sapply(rownames(tsne_res), function(x){return(unique(data_proteomics$Panel[data_proteomics$protein_name == x]))})
ggplot(tsne_res, aes(x = tSNE1, y = tSNE2, color = class)) + geom_point()

cor_proteins <- cor(data_proteomics_wide[data_proteomics_wide$Cancer == "BRC",inflammatory_proteins], method = "spearman")
cor_proteins <- cor(data_proteomics_wide[,inflammatory_proteins], method = "spearman")

pheatmap(cor_proteins)

hc_result <- hclust(dist(cor_proteins))

# Cut the dendrogram to get clusters (k = 2)
cutree_result <- cutree(hc_result, k = 3)

# Visualize the clusters
plot(hc_result, main = "Hierarchical Clustering Dendrogram")
rect.hclust(hc_result, k = 3, border = 2:3)

desired_cluster <- 2

# Get names of members in the desired cluster
members_in_cluster <- names(cutree_result[cutree_result == desired_cluster])
pheatmap(cor_proteins[members_in_cluster,members_in_cluster])

ggplot(data_proteomics_wide, aes(x = TNFSF11, y = TNFRSF11B, color = Cancer)) + geom_point() + geom_smooth(method = "lm")
summary(lm(TNFSF11 ~ TNFRSF11B + factor(Cancer), data = data_proteomics_wide))
main_cytokines <- c("IFNG", "IL2", "TNF", "TGFB1", "IL6")
lapply(main_cytokines, function(cytokine){
  print(cytokine)
  return(rownames(cor_proteins)[cor_proteins[cytokine, ] > 0.3])
})

