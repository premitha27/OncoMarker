# ==============================================================================
# STEP 3: S4 CLASS DEFINITION
# ==============================================================================

# 1. Define the Generic Accessors (Section 3.3)
# We define these empty functions first so R knows they exist.
setGeneric("exprs", function(object) standardGeneric("exprs"))
setGeneric("meta", function(object) standardGeneric("meta"))
setGeneric("cohort", function(object) standardGeneric("cohort"))

# 2. Define the Class Schema (Section 3.2.1)
# This strictly enforces that expression_data MUST be a matrix, etc.
setClass("GenePanel",
         slots = c(
           expression_data = "matrix",
           patient_metadata = "data.frame",
           cancer_type = "character"
         ),
         # Define a default empty state (Prototype)
         prototype = list(
           expression_data = matrix(0, 0, 0),
           patient_metadata = data.frame(),
           cancer_type = "Unknown"
         )
)

# 3. Define the Validity Checker (Section 3.2.2)
# This function runs AUTOMATICALLY every time you say new("GenePanel", ...)
setValidity("GenePanel", function(object) {
  errors <- character()

  # Check A: Dimensional Consistency (Crucial!)
  # The number of columns in the matrix (patients) must match rows in metadata
  if (ncol(object@expression_data) != nrow(object@patient_metadata)) {
    msg <- paste("Dimension Mismatch: expression_data has",
                 ncol(object@expression_data), "columns but patient_metadata has",
                 nrow(object@patient_metadata), "rows.")
    errors <- c(errors, msg)
  }

  # Check B: Label Consistency
  # Ensure the patient IDs actually match
  if (!all(colnames(object@expression_data) == rownames(object@patient_metadata))) {
    errors <- c(errors, "Sample IDs do not align! Check colnames(matrix) vs rownames(metadata).")
  }

  # Check C: Data Sufficiency (The Diagnosis Dichotomy)
  # We require the metadata to have the column "Diagnosis" we created in Step 2
  if (!"Diagnosis" %in% colnames(object@patient_metadata)) {
    errors <- c(errors, "Metadata is missing the required 'Diagnosis' column.")
  } else {
    # Check that we have both Tumor and Normal samples to do a t-test
    levels_present <- unique(object@patient_metadata$Diagnosis)
    if (length(levels_present) < 2) {
      errors <- c(errors, "Diagnosis column must have at least two levels (e.g., Tumor and Normal).")
    }
  }

  # Return errors if any found, otherwise TRUE
  if (length(errors) == 0) TRUE else errors
})

# 4. Implement the Accessor Methods (Section 3.3)
# These allow users to get data without using the '@' symbol.

#' @describeIn GenePanel Get the expression matrix
setMethod("exprs", "GenePanel", function(object) {
  return(object@expression_data)
})

#' @describeIn GenePanel Get the clinical metadata
setMethod("meta", "GenePanel", function(object) {
  return(object@patient_metadata)
})

#' @describeIn GenePanel Get the cancer type label
setMethod("cohort", "GenePanel", function(object) {
  return(object@cancer_type)
})

print("Class 'GenePanel' defined successfully with strict validation.")
