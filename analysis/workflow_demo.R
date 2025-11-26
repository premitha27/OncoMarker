# ==============================================================================
# OncoMarker Workflow Demo
# ==============================================================================

# 1. Setup
# Load your newly installed package
library(OncoMarker)

# 2. Load Raw Data
print("Loading Data...")

# FIX: We use 'read.delim' instead of 'read.table' to strictly enforce Tab separation
normal_df <- read.delim("raw_data/BC-TCGA-Normal.txt", header = TRUE, row.names = 1, check.names = FALSE)
tumor_df  <- read.delim("raw_data/BC-TCGA-Tumor.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Verify it loaded correctly
print(paste("Normal Samples:", ncol(normal_df)))
print(paste("Tumor Samples:", ncol(tumor_df)))

# 3. Prepare the Matrix (ETL)
# Bind columns: Normal first, then Tumor
expression_matrix <- as.matrix(cbind(normal_df, tumor_df))

# CRITICAL FIX: The data loaded as "Text" (Characters). We must force it to be "Numeric".
# This removes the quotes so R can do math on it.
storage.mode(expression_matrix) <- "numeric"

# 4. Prepare Metadata
# We manually create the labels because we know the order of binding
n_normal <- ncol(normal_df)
n_tumor  <- ncol(tumor_df)

meta_df <- data.frame(
  SampleID = colnames(expression_matrix),
  Diagnosis = factor(c(rep("Normal", n_normal), rep("Tumor", n_tumor))),
  row.names = colnames(expression_matrix)
)

# 5. Subset to Targeted Panel (PAM50 + Biomarkers)
# We filter the 17,000 genes down to just the ones we care about
target_genes <- c("ESR1", "PGR", "ERBB2", "MKI67", "TP53", "BRCA1", "BRCA2", "PTEN",
                  "AKT1", "PIK3CA", "MYC", "CCND1", "EGFR", "CDH1", "ATM", "CHEK2")

# Find which genes exist in the file
genes_found <- intersect(target_genes, rownames(expression_matrix))
final_matrix <- expression_matrix[genes_found, ]

# 6. Create the OncoMarker Object
# This uses the S4 class we built
brca_panel <- new("GenePanel",
                  expression_data = final_matrix,
                  patient_metadata = meta_df,
                  cancer_type = "TCGA-BRCA")

print("Object created successfully!")
# Show a summary of the object
show(brca_panel)

# 7. Run Analysis (Differential Expression)
print("Calculating Fold Changes...")
results <- calc_fold_change(brca_panel)
head(results)

# 8. Visualization (Volcano Plot)
print("Generating Volcano Plot...")
p <- plot_volcano(brca_panel)
print(p)

# 9. Risk Prediction (Clinical Stratification)
print("Predicting Risk for TP53...")
risk_scores <- predict_risk(brca_panel, gene = "TP53", direction = "low_risk_high_expr")
table(risk_scores)
