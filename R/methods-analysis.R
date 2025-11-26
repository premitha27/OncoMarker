# ==============================================================================
# STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================
#' Calculate Fold Change
#'
#' Performs differential expression analysis using Welch's t-test.
#' @param object A GenePanel object.
#' @return A dataframe of results with Log2FC and P-values.
#' @export
#' @docType methods
#' @rdname calc_fold_change-methods
# 1. Define the Generic Method
# This tells R that "calc_fold_change" is a function that can behave differently
# depending on the object it is given.
setGeneric("calc_fold_change", function(object) standardGeneric("calc_fold_change"))

# 2. Implement the Method for GenePanel
#' @param object A GenePanel S4 object
#' @return A dataframe with Gene, Log2FC, PValue, and FDR
setMethod("calc_fold_change", "GenePanel", function(object) {

  # A. Extract Data using our Accessors
  mat <- exprs(object)
  meta <- meta(object)

  print(paste("Analyzing", nrow(mat), "genes across", ncol(mat), "samples..."))

  # B. Identify Groups
  # We use regex to be robust against "Primary Tumor" vs "Tumor"
  tumor_idx <- which(grepl("Tumor|Malignant", meta$Diagnosis, ignore.case = TRUE))
  normal_idx <- which(grepl("Normal|Benign", meta$Diagnosis, ignore.case = TRUE))

  # Safety Check: Need at least 3 samples per group for a valid t-test
  if (length(tumor_idx) < 3 || length(normal_idx) < 3) {
    stop("Error: Insufficient samples (n < 3) in Tumor or Normal groups.")
  }

  # C. Pre-allocate Results Container
  genes <- rownames(mat)
  results <- data.frame(
    Gene = genes,
    Log2FC = numeric(length(genes)),
    PValue = numeric(length(genes)),
    stringsAsFactors = FALSE
  )

  # D. The Calculation Loop
  # We iterate through every gene in the panel
  for (i in seq_along(genes)) {

    # Extract the vector of expression values for the current gene
    tumor_vals <- mat[i, tumor_idx]
    normal_vals <- mat[i, normal_idx]

    # 1. Calculate Log2 Fold Change
    # Logic: Log(A/B) = Log(A) - Log(B)
    # We assume input data is already log-transformed (e.g., Log2(TPM+1))
    results$Log2FC[i] <- mean(tumor_vals, na.rm=TRUE) - mean(normal_vals, na.rm=TRUE)

    # 2. Perform Welch's T-test
    # We use tryCatch because if a gene has 0 variance (all zeros), t.test throws an error.
    test_res <- tryCatch({
      t.test(tumor_vals, normal_vals, var.equal = FALSE) # var.equal=FALSE = Welch's
    }, error = function(e) {
      return(list(p.value = NA)) # Handle edge case gracefully
    })

    results$PValue[i] <- test_res$p.value
  }

  # E. Multiple Testing Correction (Benjamini-Hochberg)
  # This adjusts for the fact that we asked 50 questions (50 genes), reducing false positives.
  results$FDR <- p.adjust(results$PValue, method = "BH")

  # Sort by statistical significance (most significant at top)
  results <- results[order(results$PValue), ]

  return(results)
})

print("Method 'calc_fold_change' registered.")
