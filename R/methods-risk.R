# ==============================================================================
# STEP 6: RISK PREDICTION (The "Prognostic Model")
# ==============================================================================
#' Predict Risk
#'
#' Stratifies patients based on biomarker expression.
#'
#' @param object A GenePanel object.
#' @param gene Character string. The gene symbol to use (must be in the matrix).
#' @param direction String. Either 'low_risk_high_expr' (Tumor Suppressors) or 'high_risk_high_expr' (Oncogenes).
#' @return A factor vector of risk labels ("High Risk", "Low Risk").
#' @export
#' @docType methods
#' @rdname predict_risk-methods
#' @name predict_risk
# 1. Define the Generic
setGeneric("predict_risk", function(object, gene, direction = "low_risk_high_expr")
  standardGeneric("predict_risk"))

# 2. Implement the Method
#' @rdname predict_risk-methods
#' @aliases predict_risk,GenePanel-method
setMethod("predict_risk", "GenePanel", function(object, gene, direction) {

  # A. Validation
  # Ensure the requested gene is actually in our 50-gene panel
  if (!gene %in% rownames(exprs(object))) {
    stop(paste("Error: Gene", gene, "not found in the GenePanel."))
  }

  # B. Extract Expression Data
  # We grab the row of data for just this one gene across all patients
  expr_vals <- exprs(object)[gene, ]

  # C. Determine the Cutoff
  # We use the median of the cohort as the dividing line.
  cutoff <- median(expr_vals, na.rm = TRUE)
  print(paste("Stratifying cohort based on", gene, "median expression:", round(cutoff, 2)))

  # D. Apply Logic based on Biological Function
  risk_vector <- rep(NA, length(expr_vals))

  if (direction == "low_risk_high_expr") {
    # --- TUMOR SUPPRESSOR LOGIC (e.g., TP53, BRCA1) ---
    # Biology: We need these genes. Loss of them is dangerous.
    # Logic: If expression is LOW (< median), risk is HIGH.
    risk_vector <- ifelse(expr_vals < cutoff, "High Risk", "Low Risk")

  } else {
    # --- ONCOGENE LOGIC (e.g., HER2/ERBB2, MYC) ---
    # Biology: These genes drive cancer. Too much is dangerous.
    # Logic: If expression is HIGH (> median), risk is HIGH.
    risk_vector <- ifelse(expr_vals > cutoff, "High Risk", "Low Risk")
  }

  # E. Return as a Factor
  # Setting levels ensures "Low Risk" is the reference level (0) for stats
  return(factor(risk_vector, levels = c("Low Risk", "High Risk")))
})

print("Method 'predict_risk' registered.")
