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

# Define a function to rename columns and clean their content
clean_and_rename <- function(data_frame) {
  # Define a function to clean column values
  clean_column <- function(column_values) {
    # Use regular expressions to extract the relevant information
    cleaned_values <- gsub(".*: ", "", column_values)
    return(cleaned_values)
  }
  
  # Rename columns based on their content
  column_names <- c("title", "patient_identifier", "tumor_number", "primary_progression", 
                    "tumor_stage", "grade", "molecular_subtype")
  colnames(data_frame) <- column_names
  
  # Apply the clean_column function to each column
  for (i in 2:ncol(data_frame)) {
    data_frame[, i] <- clean_column(data_frame[, i])
  }
  
  return(data_frame)
}
