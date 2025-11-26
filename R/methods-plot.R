# ==============================================================================
# STEP 5: VISUALIZATION (The "Volcano Plot")
# ==============================================================================
#' Volcano Plot
#'
#' Generates a volcano plot for differential expression results.
#' @param object A GenePanel object.
#' @return A ggplot object.
#' @export
#' @docType methods
#' @rdname plot_volcano-methods
library(ggplot2)

# 1. Define the Generic
setGeneric("plot_volcano", function(object) standardGeneric("plot_volcano"))

# 2. Implement the Method
#' @param object A GenePanel S4 object
#' @return A ggplot object
setMethod("plot_volcano", "GenePanel", function(object) {

  # A. Get the Stats
  # We call the method we wrote in Step 4 internally
  stats <- calc_fold_change(object)

  # B. Define Significance Categories (Section 5.2)
  # We classify every gene into one of three buckets for coloring
  stats$Category <- "Not Significant"

  # Upregulated: Log2FC > 1 (2x higher) AND Significant (p < 0.05)
  stats$Category[stats$Log2FC > 1 & stats$PValue < 0.05] <- "Upregulated"

  # Downregulated: Log2FC < -1 (2x lower) AND Significant (p < 0.05)
  stats$Category[stats$Log2FC < -1 & stats$PValue < 0.05] <- "Downregulated"

  # C. Construct the Plot using ggplot2
  p <- ggplot(stats, aes(x = Log2FC, y = -log10(PValue), color = Category, label = Gene)) +

    # Geometry: Scatter plot points
    geom_point(alpha = 0.7, size = 2) +

    # Colors: Strict Red/Blue mapping (Section 5.2)
    scale_color_manual(values = c(
      "Upregulated" = "#E41A1C",   # Red
      "Downregulated" = "#377EB8", # Blue
      "Not Significant" = "#999999" # Grey
    )) +

    # Threshold Lines: Visual cues for the cutoffs
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Fold Change Cutoff
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # P-value Cutoff

    # Labels and Theme
    labs(
      title = paste("Volcano Plot:", cohort(object)),
      subtitle = "Thresholds: p < 0.05 and |Log2FC| > 1",
      x = "Log2 Fold Change (Tumor/Normal)",
      y = "-Log10 P-Value"
    ) +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)
})

print("Method 'plot_volcano' registered.")
