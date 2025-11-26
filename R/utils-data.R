#' Standardize Clinical Metadata
#'
#' Implements the logic from Section 2.3 of the OncoMarker architecture.
#' It scans the metadata for diagnosis labels and maps them to a strict
#' Tumor/Normal factor for downstream differential analysis.
#'
#' @param metadata A data.frame containing patient clinical data.
#' @return A data.frame with a standardized 'Diagnosis' column.
#' @export
standardize_metadata <- function(metadata) {

  # 1. Identify the Diagnosis Column
  # We look for common column names used in Kaggle/TCGA datasets
  potential_cols <- c("Diagnosis", "diagnosis", "Sample_Type", "sample_type",
                      "Label", "label", "Condition", "condition", "Target")

  # Find the first match in the dataframe
  match_col <- intersect(potential_cols, colnames(metadata))[1]

  if (is.na(match_col)) {
    stop("Error: Could not find a recognizable diagnosis column (e.g., 'Diagnosis', 'Label', 'Sample_Type').")
  }

  print(paste("Found diagnosis information in column:", match_col))

  # Extract the raw vector
  raw_labels <- as.character(metadata[[match_col]])
  normalized_labels <- rep(NA, length(raw_labels))

  # 2. Apply the Mapping Logic (The "Diagnosis Dichotomy")
  # Mappings defined in Table 2.3 of the design doc

  # Group A: Tumor Mappings (Numerator in FC)
  # Matches: "Primary Tumor", "Malignant", "Metastatic", "Glioblastoma", "01" (TCGA code)
  tumor_regex <- "Tumor|Malignant|Metastatic|Glioblastoma|Cancer|Positive|^01$"
  normalized_labels[grepl(tumor_regex, raw_labels, ignore.case = TRUE)] <- "Tumor"

  # Group B: Normal Mappings (Denominator in FC)
  # Matches: "Solid Tissue Normal", "Benign", "Control", "Healthy", "11" (TCGA code)
  normal_regex <- "Normal|Benign|Control|Healthy|Negative|^11$"
  normalized_labels[grepl(normal_regex, raw_labels, ignore.case = TRUE)] <- "Normal"

  # 3. Validation
  if (any(is.na(normalized_labels))) {
    warning("Some samples could not be classified and were set to NA. Check input labels.")
  }

  # 4. Standardize Structure
  # Create the strict factor required for t-tests
  metadata$Diagnosis <- factor(normalized_labels, levels = c("Normal", "Tumor"))

  return(metadata)
}
