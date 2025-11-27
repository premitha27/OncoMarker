# OncoMarker
**Targeted Genomic Biomarker Discovery Framework**<br>
**Author:** Premitha Pagadala<br>
**Email:** premitha@usf.edu<br>
**Blog Link:** [Final Project Package Blog Post](https://premithapagadala.blogspot.com/2025/11/assignment-11.html)

## Overview
OncoMarker is an R package designed for targeted biomarker discovery, differential expression analysis, volcano-plot visualization, and simple risk stratification across cancer cohorts.<br>
Unlike traditional RNA-seq pipelines requiring FASTQ → alignment → quantification, OncoMarker operates directly on pre-processed Level-3 expression matrices (CSV / RData).
This dramatically reduces compute requirements and makes biomarker discovery accessible to researchers using standard laptops.

## Key Features
- Accepts normalized gene-expression matrices (TPM, FPKM, RSEM, etc.)
- Strict S4 GenePanel class ensures integrity of genomic + clinical data
- Differential expression (Welch’s t-test + FDR correction)
- Publication-quality Volcano plots
- Gene-based risk prediction using tumor suppressor / oncogene logic
- Fully reproducible & Bioconductor-style architecture

## Installation
- Install devtools in R Studio
- Enter the following commands in the console:
```r
devtools::install_github("premitha27/OncoMarker")
library(OncoMarker)
```

## Before Using the Package – Important Setup
Set your working directory to the data folder. You can do this by going into the top toolbar in R:<br>
**Session --> Set Working Directory --> Choose Directory**<br>
Choose the data directory which you downloaded from this repository.<br>
Make sure all the data sets or files you want to analyze are in .txt format

## The GenePanel Object
All analysis starts with:
```r
gp <- new("GenePanel",
          expression_data = expr_matrix,
          patient_metadata = metadata,
          cancer_type = "TCGA-BRCA")
```
**Requirements:**
Rows = genes, columns = patients
Metadata rownames must match patient IDs
Metadata must include Diagnosis (e.g., Tumor/Normal)

## Functions
The OncoMarker package currently includes 3 core functions:
calc_fold_change() – differential expression using Welch’s t-test + FDR
**plot_volcano()** – publication-ready volcano plot
**predict_risk()** – simple risk stratification using gene expression

## Differential Expression
```r
results <- calc_fold_change(gp)
head(results)
```

## Volcano Plot
```r
plot_volcano(gp)
```
Automatically highlights:
Significant upregulated genes
Significant downregulated genes
Uses |Log2FC| > 1 and p < 0.05 thresholds

## Risk Prediction
Works for tumor suppressors or oncogenes.
```r
risk <- predict_risk(gp, gene = "TP53", direction = "low_risk_high_expr")
gp@patient_metadata$Risk <- risk
```
Labels each patient as “Low Risk” or “High Risk” based on median expression.

## Data Source
Kaggle Dataset: Gene Expression Profiles of Breast Cancer
https://www.kaggle.com/datasets/orvile/gene-expression-profiles-of-breast-cancer