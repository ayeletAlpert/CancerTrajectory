#  FUNCTIONS FOR MAPPING GENE SYMBOLS TO EMSEBLE IDS ------------------------

# Connect to Ensembl database
ensembl <- useMart("ensembl")

# Specify the dataset and species
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Get all Ensembl gene IDs
genes <- getBM(attributes = c("ensembl_gene_id"), mart = ensembl)

res_symbol <- mapIds(org.Hs.eg.db, keys <- gsub("\\.\\d+", "", genes$ensembl_gene_id), 
                     column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
ensembl2sym <- function(ensembl){return(res_symbol[ensembl])}
res_op_sym <- names(res_symbol)
names(res_op_sym) <- res_symbol
sym2ensembl <- function(sym){return(res_op_sym[sym])}
