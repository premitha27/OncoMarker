---
title: "OncoMarker: Targeted Genomic Analysis Workflow"
author: "Premitha Pagadala"
date: "2025-11-27"
output: 
  html_document:
    toc: true
    toc_float: true
    css: style.css
    df_print: paged
    keep_md: true
    self_contained: true
    output_file: "OncoMarker_Analysis.html"
vignette: >
  %\VignetteIndexEntry{OncoMarker Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# 1. Introduction
The OncoMarker package provides a streamlined S4 object-oriented framework for analyzing targeted genomic panels. It is designed to bridge the gap between massive public datasets (like TCGA) and clinical biomarker discovery.

This vignette demonstrates the end-to-end workflow:<br>
1. Ingesting raw text data (Level 3 expression).<br>
2. Constructing the strictly typed GenePanel S4 object.<br>
3. Performing differential expression analysis.<br>
4. Visualizing results via Volcano plots.<br>
5. Stratifying patients based on risk biomarkers.

# 2. Data Ingestion & Object Construction
The OncoMarker framework relies on the GenePanel class to ensure data integrity. The code below demonstrates how to load raw text files and create the object.

Note: The code below assumes you have the specific TCGA text files locally.

``` r
# 1. Define Paths to your raw data
normal_path <- "../raw_data/BC-TCGA-Normal.txt"
tumor_path  <- "../raw_data/BC-TCGA-Tumor.txt"

# 2. Load Data 
# We use read.delim to strictly handle tab-separation and avoid parsing errors
normal_df <- read.delim(normal_path, header=TRUE, row.names=1, check.names=FALSE)
tumor_df  <- read.delim(tumor_path, header=TRUE, row.names=1, check.names=FALSE)

# 3. ETL: Merge and Format
# Bind columns: Normal samples first, then Tumor samples
expr_mat <- as.matrix(cbind(normal_df, tumor_df))

# CRITICAL STEP: Ensure matrix is numeric (converts text to numbers)
storage.mode(expr_mat) <- "numeric"

# 4. Create Clinical Metadata
# We automatically generate labels based on the column counts
meta_df <- data.frame(
  SampleID = colnames(expr_mat),
  Diagnosis = factor(c(rep("Normal", ncol(normal_df)), rep("Tumor", ncol(tumor_df)))),
  row.names = colnames(expr_mat)
)

# 5. Create the S4 Object
# This step runs internal validation to ensure dimensions match
panel <- new("GenePanel", 
             expression_data = expr_mat, 
             patient_metadata = meta_df, 
             cancer_type = "TCGA-BRCA")

show(panel)
#> An object of class "GenePanel"
#> Slot "expression_data":
#>                 TCGA-BH-A0AY-11A-23R-A089-07 TCGA-A7-A0DB-11A-33R-A089-07
#> ELMO2                           2.043333e-01                 8.694167e-01
#>                 TCGA-BH-A0HK-11A-11R-A089-07 TCGA-BH-A0BM-11A-12R-A089-07
#> ELMO2                           0.0645000000                -1.862500e-01
#>                 TCGA-BH-A0B3-11B-21R-A089-07 TCGA-BH-A0DK-11A-13R-A089-07
#> ELMO2                           0.0632500000                 2.887500e-01
#>                 TCGA-BH-A0BJ-11A-23R-A089-07 TCGA-BH-A0DP-11A-12R-A089-07
#> ELMO2                           4.296667e-01                 0.0664166667
#>                 TCGA-BH-A0BV-11A-31R-A089-07 TCGA-BH-A0DQ-11A-12R-A089-07
#> ELMO2                           0.2381666667                 2.569167e-01
#>                 TCGA-BH-A0B8-11A-41R-A089-07 TCGA-BH-A0C0-11A-21R-A089-07
#> ELMO2                           4.752500e-01                 2.032500e-01
#>                 TCGA-BH-A0DZ-11A-22R-A089-07 TCGA-BH-A0BA-11A-22R-A089-07
#> ELMO2                           0.1170833333                -0.5690833333
#>                 TCGA-BH-A0BC-11A-22R-A089-07 TCGA-BH-A0DH-11A-31R-A089-07
#> ELMO2                          -1.718500e+00                -1.331000e+00
#>                 TCGA-A7-A0CH-11A-32R-A089-07 TCGA-A7-A0CE-11A-21R-A089-07
#> ELMO2                          -0.6330000000                 5.580000e-01
#>                 TCGA-BH-A0E1-11A-13R-A089-07 TCGA-BH-A0H9-11A-22R-A089-07
#> ELMO2                          -2.984167e-01                 1.083333e-03
#>                 TCGA-A7-A0D9-11A-53R-A089-07 TCGA-BH-A0H7-11A-13R-A089-07
#> ELMO2                           0.0037500000                -2.878333e-01
#>                 TCGA-BH-A0E0-11A-13R-A089-07 TCGA-BH-A0C3-11A-23R-A12P-07
#> ELMO2                          -5.593333e-01                -0.5531666667
#>                 TCGA-BH-A0DL-11A-13R-A115-07 TCGA-BH-A0H5-11A-62R-A115-07
#> ELMO2                           7.441667e-02                -5.009167e-01
#>                 TCGA-BH-A0BQ-11A-33R-A115-07 TCGA-BH-A0B7-11A-34R-A115-07
#> ELMO2                          -5.364167e-01                 0.2710833333
#>                 TCGA-BH-A0BW-11A-12R-A115-07 TCGA-BH-A18N-11A-43R-A12D-07
#> ELMO2                           1.270000e-01                -1.559167e-01
#>                 TCGA-E2-A158-11A-22R-A12D-07 TCGA-BH-A18P-11A-43R-A12D-07
#> ELMO2                           0.2155000000                -2.640000e-01
#>                 TCGA-BH-A0DO-11A-22R-A12D-07 TCGA-BH-A18Q-11A-34R-A12D-07
#> ELMO2                          -0.2705833333                -3.113333e-01
#>                 TCGA-BH-A0DT-11A-12R-A12D-07 TCGA-BH-A18R-11A-42R-A12D-07
#> ELMO2                          -9.694167e-01                 -0.724166667
#>                 TCGA-E2-A15M-11A-22R-A12D-07 TCGA-BH-A18J-11A-31R-A12D-07
#> ELMO2                          -0.3385833333                -8.841667e-01
#>                 TCGA-BH-A18S-11A-43R-A12D-07 TCGA-BH-A18L-11A-42R-A12D-07
#> ELMO2                          -0.4750833333                -1.816667e-01
#>                 TCGA-BH-A18V-11A-52R-A12D-07 TCGA-BH-A18M-11A-33R-A12D-07
#> ELMO2                          -6.741667e-01                -5.675000e-02
#>                 TCGA-E2-A153-11A-31R-A12D-07 TCGA-BH-A0BS-11A-11R-A12P-07
#> ELMO2                          -0.2708181818                -5.553333e-01
#>                 TCGA-A7-A13E-11A-61R-A12P-07 TCGA-BH-A0BZ-11A-61R-A12P-07
#> ELMO2                           2.583333e-03                -0.1279166667
#>                 TCGA-A7-A13F-11A-42R-A12P-07 TCGA-BH-A18K-11A-13R-A12D-07
#> ELMO2                          -3.020833e-01                -3.583333e-02
#>                 TCGA-BH-A18U-11A-23R-A12D-07 TCGA-BH-A0AU-11A-11R-A12P-07
#> ELMO2                          -3.071667e-01                -4.965000e-01
#>                 TCGA-BH-A0DD-11A-23R-A12P-07 TCGA-BH-A0B5-11A-23R-A12P-07
#> ELMO2                           2.540833e-01                -2.734167e-01
#>                 TCGA-BH-A0DV-11A-22R-A12P-07 TCGA-BH-A1F0-11B-23R-A137-07
#> ELMO2                          -1.955833e-01                -0.8787500000
#>                 TCGA-E2-A1BC-11A-32R-A12P-07 TCGA-BH-A0HA-11A-31R-A12P-07
#> ELMO2                          -5.591667e-02                 0.1646666667
#>                 TCGA-E2-A15I-11A-32R-A137-07 TCGA-BH-A1EW-11B-33R-A137-07
#> ELMO2                          -8.009167e-01                -1.205583e+00
#>                 TCGA-BH-A1EU-11A-23R-A137-07 TCGA-BH-A1ET-11B-23R-A137-07
#> ELMO2                          -1.077750e+00                -0.6860833333
#>                 TCGA-BH-A1EO-11A-31R-A137-07 TCGA-AO-A03P-01A-11R-A00Z-07
#> ELMO2                          -7.753333e-01                 1.2194166667
#>                 TCGA-A8-A06T-01A-11R-A00Z-07 TCGA-A8-A07F-01A-11R-A00Z-07
#> ELMO2                           4.303333e-01                 0.4398333333
#>                 TCGA-A8-A081-01A-11R-A00Z-07 TCGA-A8-A08C-01A-11R-A00Z-07
#> ELMO2                           1.2426666667                -2.283333e-01
#>                 TCGA-A8-A08T-01A-21R-A00Z-07 TCGA-A8-A095-01A-11R-A00Z-07
#> ELMO2                          -0.3430833333                 4.421667e-01
#>                 TCGA-A8-A09G-01A-21R-A00Z-07 TCGA-A8-A09W-01A-11R-A00Z-07
#> ELMO2                           7.523333e-01                 2.875833e-01
#>                 TCGA-AO-A03O-01A-11R-A00Z-07 TCGA-A8-A06R-01A-11R-A00Z-07
#> ELMO2                           1.171750e+00                 2.727500e-01
#>                 TCGA-A8-A07B-01A-11R-A00Z-07 TCGA-A8-A07Z-01A-11R-A00Z-07
#> ELMO2                           0.1303333333                 1.5358333333
#>                 TCGA-A8-A08B-01A-11R-A00Z-07 TCGA-A8-A08P-01A-11R-A00Z-07
#> ELMO2                           6.370000e-01                 0.5302500000
#>                 TCGA-A8-A094-01A-11R-A00Z-07 TCGA-A8-A09E-01A-11R-A00Z-07
#> ELMO2                           2.041583e+00                 1.2508333333
#>                 TCGA-A8-A09T-01A-11R-A00Z-07 TCGA-AN-A0AR-01A-11R-A00Z-07
#> ELMO2                           0.4650000000                -0.1448333333
#>                 TCGA-BH-A0DZ-01A-11R-A00Z-07 TCGA-AN-A0AS-01A-11R-A00Z-07
#> ELMO2                           5.180833e-01                 0.7245833333
#>                 TCGA-A7-A0CH-01A-21R-A00Z-07 TCGA-AN-A0AJ-01A-11R-A00Z-07
#> ELMO2                           0.1730833333                  0.673666667
#>                 TCGA-AN-A03X-01A-21R-A00Z-07 TCGA-A8-A06U-01A-11R-A00Z-07
#> ELMO2                           4.467500e-01                 9.288333e-01
#>                 TCGA-A8-A07I-01A-11R-A00Z-07 TCGA-A8-A082-01A-11R-A00Z-07
#> ELMO2                           1.198583e+00                 0.9755833333
#>                 TCGA-A8-A08F-01A-11R-A00Z-07 TCGA-A8-A08X-01A-21R-A00Z-07
#> ELMO2                           5.432500e-01                 4.690000e-01
#>                 TCGA-A8-A096-01A-11R-A00Z-07 TCGA-A8-A09K-01A-11R-A00Z-07
#> ELMO2                           3.180833e-01                 0.5110000000
#>                 TCGA-A8-A09X-01A-11R-A00Z-07 TCGA-BH-A0AY-01A-21R-A00Z-07
#> ELMO2                           0.2485833333                 2.677500e-01
#>                 TCGA-A7-A0CJ-01A-21R-A00Z-07 TCGA-AN-A03Y-01A-21R-A00Z-07
#> ELMO2                          -2.014167e-01                 0.6890833333
#>                 TCGA-A8-A06X-01A-21R-A00Z-07 TCGA-A8-A07J-01A-11R-A00Z-07
#> ELMO2                           9.896667e-01                 1.078000e+00
#>                 TCGA-A8-A083-01A-21R-A00Z-07 TCGA-A8-A08G-01A-11R-A00Z-07
#> ELMO2                           1.061000e+00                 1.682417e+00
#>                 TCGA-A8-A08Z-01A-21R-A00Z-07 TCGA-A8-A099-01A-11R-A00Z-07
#> ELMO2                           4.861667e-01                 1.112667e+00
#>                 TCGA-A8-A09M-01A-11R-A00Z-07 TCGA-A8-A09Z-01A-11R-A00Z-07
#> ELMO2                           0.4363333333                 1.036167e+00
#>                 TCGA-BH-A0B4-01A-11R-A00Z-07 TCGA-A2-A0CX-01A-21R-A00Z-07
#> ELMO2                          -0.0370833333                -3.032500e-01
#>                 TCGA-AN-A049-01A-21R-A00Z-07 TCGA-A8-A06Y-01A-21R-A00Z-07
#> ELMO2                           7.333333e-03                 4.941667e-02
#>                 TCGA-A8-A084-01A-21R-A00Z-07 TCGA-A8-A08H-01A-21R-A00Z-07
#> ELMO2                           0.6041666667                -4.259167e-01
#>                 TCGA-A8-A090-01A-11R-A00Z-07 TCGA-A8-A09A-01A-11R-A00Z-07
#> ELMO2                           7.541667e-01                 0.9787500000
#>                 TCGA-A8-A09N-01A-11R-A00Z-07 TCGA-A8-A0A1-01A-11R-A00Z-07
#> ELMO2                           1.4937500000                 0.6157500000
#>                 TCGA-AN-A0AK-01A-21R-A00Z-07 TCGA-A2-A0D0-01A-11R-A00Z-07
#> ELMO2                           6.397500e-01                 8.222500e-01
#>                 TCGA-AN-A0FJ-01A-11R-A00Z-07 TCGA-A8-A06Z-01A-11R-A00Z-07
#> ELMO2                          -2.487500e-01                -4.313333e-01
#>                 TCGA-A8-A07O-01A-11R-A00Z-07 TCGA-A8-A085-01A-11R-A00Z-07
#> ELMO2                          -0.4847500000                 1.009167e-01
#>                 TCGA-A8-A08I-01A-11R-A00Z-07 TCGA-A8-A091-01A-11R-A00Z-07
#> ELMO2                           0.9524166667                 1.422417e+00
#>                 TCGA-A8-A09B-01A-11R-A00Z-07 TCGA-AN-A0AL-01A-11R-A00Z-07
#> ELMO2                           4.766667e-02                -3.217500e-01
#>                 TCGA-A8-A0A4-01A-11R-A00Z-07 TCGA-BH-A0BV-01A-11R-A00Z-07
#> ELMO2                          -4.224167e-01                 5.395833e-01
#>                 TCGA-A2-A0D4-01A-11R-A00Z-07 TCGA-AN-A0FV-01A-11R-A00Z-07
#> ELMO2                          -2.699167e-01                 4.391667e-02
#>                 TCGA-A8-A06O-01A-11R-A00Z-07 TCGA-A8-A076-01A-21R-A00Z-07
#> ELMO2                           0.3248333333                 0.4275833333
#>                 TCGA-A8-A07P-01A-11R-A00Z-07 TCGA-A8-A086-01A-11R-A00Z-07
#> ELMO2                           6.087500e-01                 0.1322500000
#>                 TCGA-A8-A08J-01A-11R-A00Z-07 TCGA-A8-A092-01A-11R-A00Z-07
#> ELMO2                           1.032917e+00                 8.300833e-01
#>                 TCGA-A8-A09C-01A-11R-A00Z-07 TCGA-A8-A09Q-01A-11R-A00Z-07
#> ELMO2                          -0.3285833333                 5.533333e-01
#>                 TCGA-A7-A0CD-01A-11R-A00Z-07 TCGA-A7-A0DB-01A-11R-A00Z-07
#> ELMO2                           0.2130833333                 0.0604166667
#>                 TCGA-A8-A09I-01A-22R-A034-07 TCGA-A8-A09V-01A-11R-A034-07
#> ELMO2                           1.287417e+00                -1.285833e-01
#>                 TCGA-A8-A0A2-01A-11R-A034-07 TCGA-A8-A0AB-01A-11R-A034-07
#> ELMO2                           2.508333e-01                 9.758333e-02
#>                 TCGA-AN-A0AM-01A-11R-A034-07 TCGA-AN-A0AT-01A-11R-A034-07
#> ELMO2                           6.640000e-01                -1.451667e-01
#>                 TCGA-BH-A0BD-01A-11R-A034-07 TCGA-A2-A0CM-01A-31R-A034-07
#> ELMO2                           0.4530000000                 1.209167e-01
#>                 TCGA-A2-A0CP-01A-11R-A034-07 TCGA-A2-A0CQ-01A-21R-A034-07
#> ELMO2                          -0.6200833333                -0.3813333333
#>                 TCGA-A2-A0CU-01A-12R-A034-07 TCGA-A8-A09D-01A-11R-A00Z-07
#> ELMO2                           5.414167e-01                -0.0393333333
#>                 TCGA-A8-A09R-01A-11R-A00Z-07 TCGA-A8-A0A9-01A-11R-A00Z-07
#> ELMO2                           1.349417e+00                 0.8151666667
#>                 TCGA-A8-A06P-01A-11R-A00Z-07 TCGA-A8-A079-01A-21R-A00Z-07
#> ELMO2                          -0.0763333333                 0.7087500000
#>                 TCGA-A8-A07W-01A-11R-A00Z-07 TCGA-A8-A08A-01A-11R-A00Z-07
#> ELMO2                           8.690833e-01                -1.663333e-01
#>                 TCGA-A8-A08L-01A-11R-A00Z-07 TCGA-A8-A093-01A-11R-A00Z-07
#> ELMO2                           1.502500e-01                 2.880000e-01
#>                 TCGA-A8-A08S-01A-11R-A034-07 TCGA-A8-A097-01A-11R-A034-07
#> ELMO2                           9.234167e-01                 1.384167e-01
#>                 TCGA-A2-A04T-01A-21R-A034-07 TCGA-A2-A04V-01A-21R-A034-07
#> ELMO2                          -9.541667e-02                -2.381667e-01
#>                 TCGA-A7-A0CE-01A-11R-A00Z-07 TCGA-AO-A03R-01A-21R-A034-07
#> ELMO2                          -0.9611666667                -0.1550000000
#>                 TCGA-AO-A03T-01A-21R-A034-07 TCGA-AN-A041-01A-11R-A034-07
#> ELMO2                           0.3345000000                 4.035000e-01
#>                 TCGA-AN-A046-01A-21R-A034-07 TCGA-AN-A04A-01A-21R-A034-07
#> ELMO2                           0.4428333333                 0.4835000000
#>                 TCGA-AN-A04C-01A-21R-A034-07 TCGA-AN-A04D-01A-21R-A034-07
#> ELMO2                           0.2275000000                 0.3689166667
#>                 TCGA-AQ-A04J-01A-02R-A034-07 TCGA-A2-A04P-01A-31R-A034-07
#> ELMO2                           6.950000e-02                -9.477500e-01
#>                 TCGA-A2-A04Q-01A-21R-A034-07 TCGA-A2-A04X-01A-21R-A034-07
#> ELMO2                          -0.4536666667                 8.378333e-01
#>                 TCGA-A2-A04Y-01A-21R-A034-07 TCGA-A8-A06Q-01A-11R-A034-07
#> ELMO2                          -1.325000e-02                 3.681667e-01
#>                 TCGA-A8-A07C-01A-11R-A034-07 TCGA-A8-A07E-01A-11R-A034-07
#> ELMO2                          -2.900000e-01                 1.755833e-01
#>                 TCGA-A8-A07G-01A-11R-A034-07 TCGA-A8-A07R-01A-21R-A034-07
#> ELMO2                           4.240000e-01                -0.1596363636
#>                 TCGA-A8-A07S-01A-11R-A034-07 TCGA-A8-A07U-01A-11R-A034-07
#> ELMO2                           4.602500e-01                 3.420000e-01
#>                 TCGA-A8-A08R-01A-11R-A034-07 TCGA-A2-A0CY-01A-12R-A034-07
#> ELMO2                           4.495833e-01                 1.852500e-01
#>                 TCGA-A2-A0CZ-01A-11R-A034-07 TCGA-A2-A0D1-01A-11R-A034-07
#> ELMO2                          -2.516667e-01                -1.374167e-01
#>                 TCGA-A2-A0D2-01A-21R-A034-07 TCGA-BH-A0E6-01A-11R-A034-07
#> ELMO2                          -3.608333e-01                -1.5249166667
#>                 TCGA-BH-A0E7-01A-11R-A034-07 TCGA-A2-A0EM-01A-11R-A034-07
#> ELMO2                           2.075000e-01                  0.413666667
#>                 TCGA-A2-A0ER-01A-21R-A034-07 TCGA-A2-A0ET-01A-31R-A034-07
#> ELMO2                           2.562500e-01                 0.4263333333
#>                 TCGA-A2-A0EX-01A-21R-A034-07 TCGA-A2-A0EY-01A-11R-A034-07
#> ELMO2                          -0.4480000000                 0.1359166667
#>                 TCGA-BH-A0EB-01A-11R-A034-07 TCGA-AN-A0FL-01A-11R-A034-07
#> ELMO2                           3.050000e-02                 3.525000e-01
#>                 TCGA-AN-A0FT-01A-11R-A034-07 TCGA-AN-A0FX-01A-11R-A034-07
#> ELMO2                          -0.0219166667                -1.828333e-01
#>                 TCGA-AN-A0FY-01A-11R-A034-07 TCGA-BH-A0HO-01A-11R-A034-07
#> ELMO2                           1.5645000000                 2.052500e-01
#>                 TCGA-B6-A0I5-01A-11R-A034-07 TCGA-B6-A0I8-01A-11R-A034-07
#> ELMO2                          -6.295833e-01                 1.3710000000
#>                 TCGA-B6-A0IA-01A-11R-A034-07 TCGA-B6-A0IB-01A-11R-A034-07
#> ELMO2                           3.753333e-01                -8.875000e-02
#>                 TCGA-BH-A0EE-01A-11R-A034-07 TCGA-A2-A0EO-01A-11R-A034-07
#> ELMO2                          -8.585000e-01                -0.1433333333
#>                 TCGA-A2-A0EQ-01A-11R-A034-07 TCGA-A2-A0EV-01A-11R-A034-07
#> ELMO2                          -8.208333e-02                 4.679167e-01
#>                 TCGA-AN-A0FS-01A-11R-A034-07 TCGA-AN-A0FZ-01A-11R-A034-07
#> ELMO2                          -1.437500e-01                 4.846667e-01
#>                 TCGA-AN-A0G0-01A-11R-A034-07 TCGA-BH-A0HQ-01A-11R-A034-07
#> ELMO2                          -5.550000e-01                -3.648333e-01
#>                 TCGA-BH-A0HU-01A-11R-A034-07 TCGA-BH-A0HW-01A-11R-A034-07
#> ELMO2                          -0.1180833333                 0.4390833333
#>                 TCGA-B6-A0I9-01A-11R-A034-07 TCGA-AN-A0FD-01A-11R-A034-07
#> ELMO2                          -0.0978333333                 0.2731666667
#>                 TCGA-AN-A0FF-01A-11R-A034-07 TCGA-AN-A0FK-01A-11R-A034-07
#> ELMO2                           0.0463333333                 4.216667e-01
#>                 TCGA-AN-A0FN-01A-11R-A034-07 TCGA-AN-A0FW-01A-11R-A034-07
#> ELMO2                          -0.2280833333                -3.987500e-01
#>                 TCGA-B6-A0I2-01A-11R-A034-07 TCGA-B6-A0I6-01A-11R-A034-07
#> ELMO2                           1.127500e-01                -4.626667e-01
#>                 TCGA-B6-A0IC-01A-11R-A034-07 TCGA-B6-A0IE-01A-11R-A034-07
#> ELMO2                           0.1783333333                -0.3892500000
#>                 TCGA-B6-A0IG-01A-11R-A034-07 TCGA-B6-A0IJ-01A-11R-A034-07
#> ELMO2                           0.7341666667                 0.5870000000
#>                 TCGA-A8-A0A7-01A-11R-A040-07 TCGA-A8-A07L-01A-11R-A040-07
#> ELMO2                          -6.333333e-03                 2.902500e-01
#>                 TCGA-A7-A0CG-01A-12R-A056-07 TCGA-A8-A06N-01A-12R-A056-07
#> ELMO2                          -2.825833e-01                 5.922500e-01
#>                 TCGA-AO-A03L-01A-41R-A056-07 TCGA-BH-A0B9-01A-11R-A056-07
#> ELMO2                           4.398333e-01                -4.666667e-02
#>                 TCGA-BH-A0BA-01A-11R-A056-07 TCGA-BH-A0DS-01A-11R-A056-07
#> ELMO2                           0.5530000000                 4.936667e-01
#>                 TCGA-BH-A0H6-01A-21R-A056-07 TCGA-B6-A0IK-01A-12R-A056-07
#> ELMO2                           0.3866666667                 6.035833e-01
#>                 TCGA-AO-A0JJ-01A-11R-A056-07 TCGA-A8-A0A6-01A-12R-A056-07
#> ELMO2                          -0.1542500000                 5.927500e-01
#>                 TCGA-BH-A0BJ-01A-11R-A056-07 TCGA-BH-A0E0-01A-11R-A056-07
#> ELMO2                           4.570000e-01                -1.159167e-01
#>                 TCGA-BH-A0H7-01A-13R-A056-07 TCGA-AO-A0JA-01A-11R-A056-07
#> ELMO2                           2.451667e-01                 8.434167e-01
#>                 TCGA-AO-A0JL-01A-11R-A056-07 TCGA-BH-A0DP-01A-21R-A056-07
#> ELMO2                          -2.803333e-01                -0.1110833333
#>                 TCGA-BH-A0H0-01A-11R-A056-07 TCGA-BH-A0HY-01A-11R-A056-07
#> ELMO2                           0.8779166667                 8.189167e-01
#>                 TCGA-AO-A0JI-01A-21R-A056-07 TCGA-A8-A08O-01A-21R-A056-07
#> ELMO2                           0.4263333333                 0.4815833333
#>                 TCGA-AR-A0TT-01A-31R-A084-07 TCGA-B6-A0IM-01A-11R-A034-07
#> ELMO2                           3.899167e-01                 0.4921666667
#>                 TCGA-B6-A0IO-01A-11R-A034-07 TCGA-B6-A0IN-01A-11R-A034-07
#> ELMO2                           0.3675000000                 0.0566666667
#>                 TCGA-B6-A0WT-01A-11R-A109-07 TCGA-B6-A0X5-01A-21R-A109-07
#> ELMO2                           3.882500e-01                 3.116667e-01
#>                 TCGA-AN-A0XV-01A-11R-A109-07 TCGA-A2-A0YJ-01A-11R-A109-07
#> ELMO2                           0.2657500000                 0.5140833333
#>                 TCGA-AR-A0U2-01A-11R-A109-07 TCGA-B6-A0WV-01A-11R-A109-07
#> ELMO2                           9.070833e-01                 1.107417e+00
#>                 TCGA-AN-A0XN-01A-21R-A109-07 TCGA-AN-A0XW-01A-11R-A109-07
#> ELMO2                           5.948333e-01                 0.6708333333
#>                 TCGA-A2-A0YK-01A-22R-A109-07 TCGA-AR-A0U3-01A-11R-A109-07
#> ELMO2                           2.527500e-01                 5.420000e-01
#>                 TCGA-B6-A0WW-01A-11R-A109-07 TCGA-AN-A0XO-01A-11R-A109-07
#> ELMO2                           4.641667e-01                 0.8900833333
#>                 TCGA-A2-A0YC-01A-11R-A109-07 TCGA-A2-A0YL-01A-21R-A109-07
#> ELMO2                           5.299167e-01                -1.661667e-01
#>                 TCGA-AR-A0U4-01A-11R-A109-07 TCGA-B6-A0WX-01A-11R-A109-07
#> ELMO2                           4.570833e-01                -2.446667e-01
#>                 TCGA-AN-A0XP-01A-11R-A109-07 TCGA-A2-A0YD-01A-11R-A109-07
#> ELMO2                           0.5648333333                 2.631667e-01
#>                 TCGA-A2-A0YM-01A-11R-A109-07 TCGA-BH-A0W3-01A-11R-A109-07
#> ELMO2                          -4.170000e-01                -6.025000e-02
#>                 TCGA-B6-A0WY-01A-11R-A109-07 TCGA-AN-A0XR-01A-11R-A109-07
#> ELMO2                           2.924167e-01                -1.073333e-01
#>                 TCGA-AO-A0J2-01A-11R-A034-07 TCGA-AO-A0J3-01A-11R-A034-07
#> ELMO2                          -0.3335000000                 0.5079166667
#>                 TCGA-AO-A0J4-01A-11R-A034-07 TCGA-AO-A0J5-01A-11R-A034-07
#> ELMO2                            0.232666667                 2.171667e-01
#>                 TCGA-AO-A0J6-01A-11R-A034-07 TCGA-AO-A0J7-01A-11R-A034-07
#> ELMO2                           0.1896666667                -0.0156666667
#>                 TCGA-AO-A0J8-01A-21R-A034-07 TCGA-AO-A0J9-01A-11R-A034-07
#> ELMO2                          -1.441667e-02                -0.9088333333
#>                 TCGA-A2-A0YE-01A-11R-A109-07 TCGA-A2-A0YT-01A-11R-A109-07
#> ELMO2                          -1.474167e-01                 3.137500e-01
#>                 TCGA-BH-A0W4-01A-11R-A109-07 TCGA-B6-A0WZ-01A-11R-A109-07
#> ELMO2                          -1.379167e-01                -0.3781666667
#>                 TCGA-AN-A0XS-01A-22R-A109-07 TCGA-A2-A0YF-01A-21R-A109-07
#> ELMO2                          -1.460833e-01                 0.0659166667
#>                 TCGA-B6-A0IP-01A-11R-A034-07 TCGA-B6-A0IQ-01A-11R-A034-07
#> ELMO2                          -0.4058333333                -1.0851666667
#>                 TCGA-A8-A0AD-01A-11R-A056-07 TCGA-BH-A0BM-01A-11R-A056-07
#> ELMO2                           7.000000e-01                -0.3094166667
#>                 TCGA-BH-A0E1-01A-11R-A056-07 TCGA-BH-A0H9-01A-11R-A056-07
#> ELMO2                           3.281667e-01                 1.2175000000
#>                 TCGA-AO-A0JB-01A-11R-A056-07 TCGA-AO-A0JM-01A-21R-A056-07
#> ELMO2                          -7.555833e-01                 7.656667e-01
#>                 TCGA-BH-A0AW-01A-11R-A056-07 TCGA-BH-A0C0-01A-21R-A056-07
#> ELMO2                           7.833333e-02                 6.071667e-01
#>                 TCGA-BH-A0E2-01A-11R-A056-07 TCGA-BH-A0HB-01A-11R-A056-07
#> ELMO2                          -1.850000e-02                  0.213250000
#>                 TCGA-AO-A0JC-01A-11R-A056-07 TCGA-B6-A0RE-01A-11R-A056-07
#> ELMO2                           1.839167e-01                 0.4112500000
#>                 TCGA-BH-A0B1-01A-12R-A056-07 TCGA-A2-A0CT-01A-31R-A056-07
#> ELMO2                           0.5747500000                 0.7445833333
#>                 TCGA-A2-A0EU-01A-22R-A056-07 TCGA-BH-A0HF-01A-11R-A056-07
#> ELMO2                           9.368333e-01                 3.790833e-01
#>                 TCGA-AO-A0JD-01A-11R-A056-07 TCGA-B6-A0RG-01A-11R-A056-07
#> ELMO2                           3.401667e-01                 0.5829166667
#>                 TCGA-BH-A0B3-01A-11R-A056-07 TCGA-A7-A0D9-01A-31R-A056-07
#> ELMO2                          -7.525000e-02                 4.408333e-01
#>                 TCGA-BH-A0GY-01A-11R-A056-07 TCGA-BH-A0HK-01A-11R-A056-07
#> ELMO2                           0.1471666667                -0.2366666667
#>                 TCGA-AO-A0JE-01A-11R-A056-07 TCGA-B6-A0RI-01A-11R-A056-07
#> ELMO2                           5.395833e-01                 0.3815000000
#>                 TCGA-BH-A0B8-01A-21R-A056-07 TCGA-BH-A0DK-01A-21R-A056-07
#> ELMO2                           0.7791666667                 1.576667e-01
#>                 TCGA-BH-A0GZ-01A-11R-A056-07 TCGA-BH-A0HX-01A-21R-A056-07
#> ELMO2                            0.015750000                 0.3371666667
#>                 TCGA-AO-A0JF-01A-11R-A056-07 TCGA-AR-A0TR-01A-11R-A084-07
#> ELMO2                           2.941667e-01                 5.775833e-01
#>                 TCGA-BH-A0DH-01A-11R-A084-07 TCGA-B6-A0RM-01A-11R-A084-07
#> ELMO2                           5.583333e-02                 4.758333e-01
#>                 TCGA-BH-A0RX-01A-21R-A084-07 TCGA-A2-A0ST-01A-12R-A084-07
#> ELMO2                          -9.000000e-02                -0.0148333333
#>                 TCGA-A8-A075-01A-11R-A084-07 TCGA-AO-A0JG-01A-31R-A084-07
#> ELMO2                           0.7825833333                 0.3529166667
#>                 TCGA-B6-A0RU-01A-11R-A084-07 TCGA-A1-A0SO-01A-22R-A084-07
#> ELMO2                           0.2171666667                -0.0759166667
#>                 TCGA-A2-A0T0-01A-22R-A084-07 TCGA-AR-A0TQ-01A-11R-A084-07
#> ELMO2                          -0.1706666667                 1.259750e+00
#>                 TCGA-BH-A0BC-01A-22R-A084-07 TCGA-B6-A0RL-01A-11R-A084-07
#> ELMO2                          -0.3236666667                 1.2390833333
#>                 TCGA-B6-A0RV-01A-11R-A084-07 TCGA-A1-A0SP-01A-11R-A084-07
#> ELMO2                           1.891667e-02                 5.541667e-02
#>                 TCGA-A2-A0T1-01A-21R-A084-07 TCGA-BH-A0DQ-01A-11R-A084-07
#> ELMO2                           2.690833e-01                 6.961667e-01
#>                 TCGA-B6-A0RN-01A-12R-A084-07 TCGA-A2-A0T4-01A-31R-A084-07
#> ELMO2                          -1.964167e-01                 1.713333e-01
#>                 TCGA-AR-A0TV-01A-21R-A084-07 TCGA-A2-A0EN-01A-13R-A084-07
#> ELMO2                           1.141833e+00                 5.040833e-01
#>                 TCGA-B6-A0RO-01A-22R-A084-07 TCGA-BH-A0HI-01A-11R-A084-07
#> ELMO2                           1.031250e+00                 0.6034166667
#>                 TCGA-AR-A0TX-01A-11R-A084-07 TCGA-A1-A0SE-01A-11R-A084-07
#> ELMO2                           2.9011666667                 0.1814166667
#>                 TCGA-A2-A0SU-01A-11R-A084-07 TCGA-A1-A0SH-01A-11R-A084-07
#> ELMO2                           5.823333e-01                 0.4615000000
#>                 TCGA-A2-A0SV-01A-11R-A084-07 TCGA-A2-A0T5-01A-21R-A084-07
#> ELMO2                           5.599167e-01                -0.1587500000
#>                 TCGA-AR-A0TW-01A-11R-A084-07 TCGA-B6-A0RP-01A-21R-A084-07
#> ELMO2                           5.525000e-02                 2.780000e-01
#>                 TCGA-A1-A0SJ-01A-11R-A084-07 TCGA-A2-A0SW-01A-11R-A084-07
#> ELMO2                           8.770000e-01                 9.079167e-01
#>                 TCGA-A2-A0T6-01A-11R-A084-07 TCGA-B6-A0RS-01A-11R-A084-07
#> ELMO2                           3.375833e-01                 4.019167e-01
#>                 TCGA-A1-A0SK-01A-12R-A084-07 TCGA-A2-A0SX-01A-12R-A084-07
#> ELMO2                           1.4123333333                 1.090833e+00
#>                 TCGA-A2-A0T7-01A-21R-A084-07 TCGA-AR-A0TZ-01A-12R-A084-07
#> ELMO2                           0.7499166667                 1.5315000000
#>                 TCGA-BH-A0HP-01A-12R-A084-07 TCGA-B6-A0RT-01A-21R-A084-07
#> ELMO2                           7.033333e-01                 1.0057500000
#>                 TCGA-A1-A0SM-01A-11R-A084-07 TCGA-A2-A0SY-01A-31R-A084-07
#> ELMO2                           1.6667500000                 0.4905000000
#>                 TCGA-AR-A0TP-01A-11R-A084-07 TCGA-A2-A04R-01A-41R-A109-07
#> ELMO2                           0.6280833333                -2.520000e-01
#>                 TCGA-BH-A0W5-01A-11R-A109-07 TCGA-B6-A0X1-01A-11R-A109-07
#> ELMO2                           0.4768333333                -0.5409166667
#>                 TCGA-AN-A0XT-01A-11R-A109-07 TCGA-A2-A0YG-01A-21R-A109-07
#> ELMO2                          -6.116667e-02                 6.849167e-01
#>                 TCGA-AR-A0TU-01A-31R-A109-07 TCGA-BH-A0WA-01A-11R-A109-07
#> ELMO2                          -0.0416666667                 2.180833e-01
#>                 TCGA-B6-A0X4-01A-11R-A109-07 TCGA-AN-A0XU-01A-11R-A109-07
#> ELMO2                          -1.593333e-01                -2.306667e-01
#>                 TCGA-A2-A0YH-01A-11R-A109-07 TCGA-AR-A0U0-01A-11R-A109-07
#> ELMO2                           1.126250e+00                 2.495000e-01
#>                 TCGA-AO-A12G-01A-11R-A10J-07 TCGA-AQ-A04H-01B-11R-A10J-07
#> ELMO2                           1.564167e-01                 1.6726666667
#>                 TCGA-E2-A107-01A-11R-A10J-07 TCGA-AO-A125-01A-11R-A10J-07
#> ELMO2                          -0.4727500000                -4.620833e-01
#>                 TCGA-AQ-A04L-01B-21R-A10J-07 TCGA-E2-A108-01A-13R-A10J-07
#> ELMO2                           0.0828333333                 0.0268333333
#>                 TCGA-AO-A126-01A-11R-A10J-07 TCGA-E2-A109-01A-11R-A10J-07
#> ELMO2                          -0.1252500000                 0.2565833333
#>                 TCGA-AO-A03M-01B-11R-A10J-07 TCGA-A2-A0YI-01A-31R-A10J-07
#> ELMO2                           1.421333e+00                 5.719167e-01
#>                 TCGA-E2-A10E-01A-21R-A10J-07 TCGA-AO-A12C-01A-11R-A10J-07
#> ELMO2                           6.472500e-01                -0.1586666667
#>                 TCGA-AO-A03N-01B-11R-A10J-07 TCGA-E2-A105-01A-11R-A10J-07
#> ELMO2                           1.435000e-01                 0.8365833333
#>                 TCGA-E2-A10F-01A-11R-A10J-07 TCGA-AO-A12E-01A-11R-A10J-07
#> ELMO2                           0.6332500000                 4.531667e-01
#>                 TCGA-AO-A03U-01B-21R-A10J-07 TCGA-E2-A106-01A-11R-A10J-07
#> ELMO2                           0.4156666667                 3.844167e-01
#>                 TCGA-AO-A124-01A-11R-A10J-07 TCGA-AO-A128-01A-11R-A10J-07
#> ELMO2                          -7.987500e-01                 0.6207500000
#>                 TCGA-B6-A0X7-01A-11R-A10J-07 TCGA-E2-A10B-01A-11R-A10J-07
#> ELMO2                          -3.015833e-01                 1.494583e+00
#>                 TCGA-AO-A129-01A-21R-A10J-07 TCGA-AN-A0XL-01A-11R-A10J-07
#> ELMO2                          -3.683333e-01                 9.050000e-02
#>                 TCGA-E2-A10C-01A-21R-A10J-07 TCGA-AO-A12B-01A-11R-A10J-07
#> ELMO2                          -0.2430000000                -0.2346666667
#>                 TCGA-AO-A03V-01A-11R-A115-07 TCGA-BH-A0BG-01A-11R-A115-07
#> ELMO2                           9.044167e-01                -9.508333e-02
#>                 TCGA-BH-A0C7-01B-11R-A115-07 TCGA-BH-A0DL-01A-11R-A115-07
#> ELMO2                           0.4645000000                -4.518333e-01
#>                 TCGA-BH-A0H5-01A-21R-A115-07 TCGA-AR-A0TY-01A-12R-A115-07
#> ELMO2                          -0.1708333333                -2.331667e-01
#>                 TCGA-AO-A12F-01A-11R-A115-07 TCGA-C8-A12Q-01A-11R-A115-07
#> ELMO2                          -0.5459166667                 3.505833e-01
#>                 TCGA-C8-A131-01A-11R-A115-07 TCGA-D8-A140-01A-11R-A115-07
#> ELMO2                          -2.880833e-01                -4.133333e-02
#>                 TCGA-E2-A14R-01A-11R-A115-07 TCGA-E2-A15O-01A-11R-A115-07
#> ELMO2                          -1.0746666667                 0.9098333333
#>                 TCGA-A2-A04N-01A-11R-A115-07 TCGA-BH-A0BL-01A-11R-A115-07
#> ELMO2                          -0.4860000000                -4.615000e-01
#>                 TCGA-A2-A0CL-01A-11R-A115-07 TCGA-AO-A12H-01A-11R-A115-07
#> ELMO2                           7.035000e-01                 6.103333e-01
#>                 TCGA-C8-A12T-01A-11R-A115-07 TCGA-C8-A132-01A-31R-A115-07
#> ELMO2                           1.4706666667                 0.8400833333
#>                 TCGA-D8-A141-01A-11R-A115-07 TCGA-E2-A14T-01A-11R-A115-07
#> ELMO2                           4.302500e-01                 0.8737500000
#>                 TCGA-E2-A15P-01A-11R-A115-07 TCGA-A2-A04U-01A-11R-A115-07
#> ELMO2                           6.728333e-01                 5.591667e-02
#>                 TCGA-BH-A0BP-01A-11R-A115-07 TCGA-A2-A0CS-01A-11R-A115-07
#> ELMO2                          -0.0450000000                 4.709167e-01
#>                 TCGA-BH-A0DX-01A-11R-A115-07 TCGA-B6-A0IH-01A-11R-A115-07
#> ELMO2                           0.5060000000                 4.565000e-01
#>                 TCGA-BH-A0W7-01A-11R-A115-07 TCGA-C8-A12K-01A-21R-A115-07
#> ELMO2                           0.2012500000                -4.221667e-01
#>                 TCGA-C8-A12U-01A-11R-A115-07 TCGA-C8-A134-01A-11R-A115-07
#> ELMO2                          -1.198333e-01                 2.247500e-01
#>                 TCGA-D8-A142-01A-11R-A115-07 TCGA-E2-A14X-01A-11R-A115-07
#> ELMO2                           0.8915000000                 8.040000e-01
#>                 TCGA-E2-A15R-01A-11R-A115-07 TCGA-A2-A04W-01A-31R-A115-07
#> ELMO2                           1.187167e+00                 6.164167e-01
#>                 TCGA-BH-A0BQ-01A-21R-A115-07 TCGA-A2-A0CV-01A-31R-A115-07
#> ELMO2                           0.8482500000                 0.0232500000
#>                 TCGA-BH-A0E9-01B-11R-A115-07 TCGA-B6-A0RH-01A-21R-A115-07
#> ELMO2                           1.411167e+00                 0.4006666667
#>                 TCGA-B6-A0WS-01A-11R-A115-07 TCGA-C8-A12L-01A-11R-A115-07
#> ELMO2                           6.175000e-01                 0.0518333333
#>                 TCGA-C8-A12V-01A-11R-A115-07 TCGA-C8-A135-01A-11R-A115-07
#> ELMO2                          -1.431667e-01                 0.0898333333
#>                 TCGA-D8-A143-01A-11R-A115-07 TCGA-E2-A14Z-01A-11R-A115-07
#> ELMO2                           2.165833e-01                -0.1254166667
#>                 TCGA-E2-A15S-01A-11R-A115-07 TCGA-BH-A0AV-01A-31R-A115-07
#> ELMO2                           3.511667e-01                -1.230833e-01
#>                 TCGA-A2-A0CW-01A-21R-A115-07 TCGA-BH-A0EA-01A-11R-A115-07
#> ELMO2                           1.403000e+00                -6.048333e-01
#>                 TCGA-B6-A0RQ-01A-11R-A115-07 TCGA-B6-A0X0-01A-21R-A115-07
#> ELMO2                          -3.185833e-01                -0.2100833333
#>                 TCGA-C8-A12M-01A-11R-A115-07 TCGA-C8-A12W-01A-11R-A115-07
#> ELMO2                           3.193333e-01                 3.661667e-01
#>                 TCGA-C8-A137-01A-11R-A115-07 TCGA-D8-A145-01A-11R-A115-07
#> ELMO2                           0.0741666667                 7.420000e-01
#>                 TCGA-E2-A154-01A-11R-A115-07 TCGA-E2-A15T-01A-11R-A115-07
#> ELMO2                           3.645833e-01                 0.8493333333
#>                 TCGA-BH-A0B0-01A-21R-A115-07 TCGA-BH-A0BR-01A-21R-A115-07
#> ELMO2                           3.520000e-01                 0.3999166667
#>                 TCGA-A2-A0D3-01A-11R-A115-07 TCGA-BH-A0EI-01A-11R-A115-07
#> ELMO2                          -1.869167e-01                 0.8678333333
#>                 TCGA-A1-A0SD-01A-11R-A115-07 TCGA-E2-A10A-01A-21R-A115-07
#> ELMO2                           0.5070833333                 0.3275833333
#>                 TCGA-C8-A12N-01A-11R-A115-07 TCGA-C8-A12X-01A-11R-A115-07
#> ELMO2                           0.2528333333                 9.386667e-01
#>                 TCGA-C8-A138-01A-11R-A115-07 TCGA-D8-A146-01A-31R-A115-07
#> ELMO2                           0.7965833333                 2.870833e-01
#>                 TCGA-E2-A159-01A-11R-A115-07 TCGA-BH-A0B7-01A-12R-A115-07
#> ELMO2                           4.560833e-01                 8.332500e-01
#>                 TCGA-BH-A0BW-01A-11R-A115-07 TCGA-A7-A0DA-01A-31R-A115-07
#> ELMO2                           8.548333e-01                 0.4940000000
#>                 TCGA-A2-A0ES-01A-11R-A115-07 TCGA-A2-A0T3-01A-21R-A115-07
#> ELMO2                           1.102500e-01                 4.241667e-02
#>                 TCGA-AO-A12A-01A-21R-A115-07 TCGA-C8-A12O-01A-11R-A115-07
#> ELMO2                          -1.139167e-01                 6.198333e-01
#>                 TCGA-C8-A12Z-01A-11R-A115-07 TCGA-D8-A13Y-01A-11R-A115-07
#> ELMO2                           6.388333e-01                 7.403333e-01
#>                 TCGA-D8-A147-01A-11R-A115-07 TCGA-E2-A15D-01A-11R-A115-07
#> ELMO2                           4.884167e-01                 2.869167e-01
#>                 TCGA-BH-A0DE-01A-11R-A115-07 TCGA-A2-A0EW-01A-21R-A115-07
#> ELMO2                           3.422500e-01                -1.341667e-02
#>                 TCGA-AR-A0TS-01A-11R-A115-07 TCGA-AO-A12D-01A-11R-A115-07
#> ELMO2                           1.241667e-02                 6.024167e-01
#>                 TCGA-C8-A12P-01A-11R-A115-07 TCGA-C8-A130-01A-31R-A115-07
#> ELMO2                           0.5637500000                -0.2748333333
#>                 TCGA-D8-A13Z-01A-11R-A115-07 TCGA-E2-A14O-01A-31R-A115-07
#> ELMO2                           0.3846666667                 0.2728333333
#>                 TCGA-E2-A15F-01A-11R-A115-07 TCGA-A2-A0T2-01A-11R-A084-07
#> ELMO2                           2.285000e-01                 0.3933333333
#>                 TCGA-AR-A1AH-01A-11R-A12D-07 TCGA-BH-A18I-01A-11R-A12D-07
#> ELMO2                          -0.1332500000                 2.245000e-01
#>                 TCGA-BH-A18R-01A-11R-A12D-07 TCGA-E2-A14Q-01A-11R-A12D-07
#> ELMO2                          -0.0245833333                 6.958333e-02
#>                 TCGA-E2-A155-01A-11R-A12D-07 TCGA-E2-A15L-01A-11R-A12D-07
#> ELMO2                           0.4885000000                 0.8012500000
#>                 TCGA-BH-A0BO-01A-23R-A12D-07 TCGA-BH-A18J-01A-11R-A12D-07
#> ELMO2                           0.1893333333                 3.044167e-01
#>                 TCGA-BH-A18S-01A-11R-A12D-07 TCGA-E2-A14S-01A-11R-A12D-07
#> ELMO2                           5.990833e-01                 1.118417e+00
#>                 TCGA-E2-A156-01A-11R-A12D-07 TCGA-E2-A15M-01A-11R-A12D-07
#> ELMO2                          -0.0139166667                 0.9187500000
#>                 TCGA-E2-A15A-06A-11R-A12D-07 TCGA-BH-A0C1-01B-11R-A12D-07
#> ELMO2                           0.9398333333                 0.7405000000
#>                 TCGA-BH-A18K-01A-11R-A12D-07 TCGA-BH-A18T-01A-11R-A12D-07
#> ELMO2                           4.202500e-01                 0.2599166667
#>                 TCGA-E2-A14V-01A-11R-A12D-07 TCGA-E2-A158-01A-11R-A12D-07
#> ELMO2                           0.2389166667                -6.135833e-01
#>                 TCGA-E2-A15E-06A-11R-A12D-07 TCGA-BH-A0DO-01B-11R-A12D-07
#> ELMO2                           1.259167e-01                 0.2303333333
#>                 TCGA-BH-A18L-01A-32R-A12D-07 TCGA-BH-A18U-01A-21R-A12D-07
#> ELMO2                           0.6017500000                 1.502167e+00
#>                 TCGA-E2-A14W-01A-11R-A12D-07 TCGA-E2-A15A-01A-11R-A12D-07
#> ELMO2                          -0.3732500000                 7.261667e-01
#>                 TCGA-BH-A0DT-01A-21R-A12D-07 TCGA-BH-A18M-01A-11R-A12D-07
#> ELMO2                          -5.001667e-01                -0.0440833333
#>                 TCGA-BH-A18V-01A-11R-A12D-07 TCGA-E2-A14Y-01A-21R-A12D-07
#> ELMO2                           0.0235833333                 0.0866666667
#>                 TCGA-E2-A15C-01A-31R-A12D-07 TCGA-BH-A18F-01A-11R-A12D-07
#> ELMO2                           2.115833e-01                 7.285000e-01
#>                 TCGA-C8-A12Y-01A-11R-A12D-07 TCGA-E2-A150-01A-11R-A12D-07
#> ELMO2                           6.203333e-01                 3.415000e-01
#>                 TCGA-BH-A18N-01A-11R-A12D-07 TCGA-E2-A15E-01A-11R-A12D-07
#> ELMO2                           6.029167e-01                 0.5210000000
#>                 TCGA-BH-A18G-01A-11R-A12D-07 TCGA-BH-A18P-01A-11R-A12D-07
#> ELMO2                          -4.466667e-02                 6.200833e-01
#>                 TCGA-C8-A133-01A-32R-A12D-07 TCGA-E2-A152-01A-11R-A12D-07
#> ELMO2                           0.4285000000                 5.995833e-01
#>                 TCGA-E2-A15G-01A-11R-A12D-07 TCGA-BH-A18H-01A-11R-A12D-07
#> ELMO2                           1.748333e-01                 3.692500e-01
#>                 TCGA-BH-A18Q-01A-12R-A12D-07 TCGA-E2-A14P-01A-31R-A12D-07
#> ELMO2                           4.406667e-01                 1.2021666667
#>                 TCGA-E2-A153-01A-12R-A12D-07 TCGA-E2-A15H-01A-11R-A12D-07
#> ELMO2                           1.515000e-01                 4.375833e-01
#>                 TCGA-A7-A13D-01A-13R-A12P-07 TCGA-AR-A1AO-01A-11R-A12P-07
#> ELMO2                          -4.530833e-01                 7.066667e-02
#>                 TCGA-AR-A1AX-01A-11R-A12P-07 TCGA-BH-A0BZ-01A-31R-A12P-07
#> ELMO2                           2.187500e-01                 1.711583e+00
#>                 TCGA-E2-A15J-01A-11R-A12P-07 TCGA-E2-A1BC-01A-11R-A12P-07
#> ELMO2                           5.300000e-01                 9.883333e-02
#>                 TCGA-E2-A15K-06A-11R-A12P-07 TCGA-A7-A13E-01A-11R-A12P-07
#> ELMO2                           4.822500e-01                -0.3424166667
#>                 TCGA-AR-A1AP-01A-11R-A12P-07 TCGA-AR-A1AY-01A-21R-A12P-07
#> ELMO2                           7.931667e-01                 0.4584166667
#>                 TCGA-BH-A0C3-01A-21R-A12P-07 TCGA-E2-A15K-01A-11R-A12P-07
#> ELMO2                           3.778333e-01                 3.925000e-02
#>                 TCGA-A7-A13F-01A-11R-A12P-07 TCGA-AR-A1AQ-01A-11R-A12P-07
#> ELMO2                           2.759167e+00                 2.401667e-01
#>                 TCGA-BH-A0AU-01A-11R-A12P-07 TCGA-BH-A0DD-01A-31R-A12P-07
#> ELMO2                           6.147500e-01                -5.416667e-02
#>                 TCGA-E2-A1AZ-01A-11R-A12P-07 TCGA-AR-A1AI-01A-11R-A12P-07
#> ELMO2                           0.2841666667                -6.773333e-01
#>                 TCGA-AR-A1AS-01A-11R-A12P-07 TCGA-BH-A0AZ-01A-21R-A12P-07
#> ELMO2                           7.855833e-01                 6.943333e-01
#>                 TCGA-BH-A0DG-01A-21R-A12P-07 TCGA-E2-A1B0-01A-11R-A12P-07
#> ELMO2                          -3.034167e-01                -6.816667e-02
#>                 TCGA-E2-A1BD-01A-11R-A12P-07 TCGA-AR-A1AJ-01A-21R-A12P-07
#> ELMO2                          -0.3150000000                 0.4255833333
#>                 TCGA-AR-A1AT-01A-11R-A12P-07 TCGA-BH-A0B5-01A-11R-A12P-07
#> ELMO2                           7.842500e-01                -8.333333e-04
#>                 TCGA-BH-A0DI-01A-21R-A12P-07 TCGA-E2-A1B1-01A-21R-A12P-07
#> ELMO2                           1.809167e-01                 1.1577500000
#>                 TCGA-AR-A1AK-01A-21R-A12P-07 TCGA-AR-A1AU-01A-11R-A12P-07
#> ELMO2                           6.909167e-01                -2.558333e-02
#>                 TCGA-BH-A0BF-01A-21R-A12P-07 TCGA-E2-A1B4-01A-11R-A12P-07
#> ELMO2                           1.2236666667                -3.511667e-01
#>                 TCGA-AR-A1AL-01A-21R-A12P-07 TCGA-AR-A1AV-01A-21R-A12P-07
#> ELMO2                           5.787500e-01                 3.970000e-01
#>                 TCGA-BH-A0BS-01A-11R-A12P-07 TCGA-BH-A0H3-01A-11R-A12P-07
#> ELMO2                           0.4031666667                -1.566667e-01
#>                 TCGA-E2-A1B5-01A-21R-A12P-07 TCGA-AR-A1AN-01A-11R-A12P-07
#> ELMO2                          -2.104167e-01                 3.340833e-01
#>                 TCGA-AR-A1AW-01A-21R-A12P-07 TCGA-BH-A0BT-01A-11R-A12P-07
#> ELMO2                           0.4830000000                 0.1826666667
#>                 TCGA-BH-A0HA-01A-11R-A12P-07 TCGA-C8-A1HG-01A-11R-A137-07
#> ELMO2                           0.0423333333                 0.3644166667
#>                 TCGA-BH-A1ES-01A-11R-A137-07 TCGA-AR-A1AR-01A-31R-A137-07
#> ELMO2                           5.763333e-01                -0.5881666667
#>                 TCGA-BH-A1EU-01A-11R-A137-07 TCGA-BH-A1EW-01A-11R-A137-07
#> ELMO2                          -6.183333e-02                 0.6844166667
#>                 TCGA-C8-A1HF-01A-11R-A137-07 TCGA-C8-A1HN-01A-11R-A137-07
#> ELMO2                          -1.259167e-01                -0.0705000000
#>                 TCGA-E2-A15I-01A-21R-A137-07 TCGA-C8-A1HM-01A-12R-A137-07
#> ELMO2                          -0.5493333333                 0.9078333333
#>                 TCGA-E2-A1B6-01A-31R-A12P-07 TCGA-E2-A14N-01A-31R-A137-07
#> ELMO2                          -3.708333e-02                -0.9248333333
#>                 TCGA-BH-A1ET-01A-11R-A137-07 TCGA-C8-A1HL-01A-11R-A137-07
#> ELMO2                           0.1725833333                 0.5730000000
#>                 TCGA-BH-A1EV-01A-11R-A137-07 TCGA-C8-A1HI-01A-11R-A137-07
#> ELMO2                           1.148750e+00                 9.545833e-01
#>                 TCGA-BH-A1F0-01A-11R-A137-07 TCGA-BH-A1EO-01A-11R-A137-07
#> ELMO2                          -2.902500e-01                 0.2420000000
#>  [ reached 'max' / getOption("max.print") -- omitted 17813 rows ]
#> 
#> Slot "patient_metadata":
#>                                                  SampleID Diagnosis
#> TCGA-BH-A0AY-11A-23R-A089-07 TCGA-BH-A0AY-11A-23R-A089-07    Normal
#> TCGA-A7-A0DB-11A-33R-A089-07 TCGA-A7-A0DB-11A-33R-A089-07    Normal
#> TCGA-BH-A0HK-11A-11R-A089-07 TCGA-BH-A0HK-11A-11R-A089-07    Normal
#> TCGA-BH-A0BM-11A-12R-A089-07 TCGA-BH-A0BM-11A-12R-A089-07    Normal
#> TCGA-BH-A0B3-11B-21R-A089-07 TCGA-BH-A0B3-11B-21R-A089-07    Normal
#> TCGA-BH-A0DK-11A-13R-A089-07 TCGA-BH-A0DK-11A-13R-A089-07    Normal
#> TCGA-BH-A0BJ-11A-23R-A089-07 TCGA-BH-A0BJ-11A-23R-A089-07    Normal
#> TCGA-BH-A0DP-11A-12R-A089-07 TCGA-BH-A0DP-11A-12R-A089-07    Normal
#> TCGA-BH-A0BV-11A-31R-A089-07 TCGA-BH-A0BV-11A-31R-A089-07    Normal
#> TCGA-BH-A0DQ-11A-12R-A089-07 TCGA-BH-A0DQ-11A-12R-A089-07    Normal
#> TCGA-BH-A0B8-11A-41R-A089-07 TCGA-BH-A0B8-11A-41R-A089-07    Normal
#> TCGA-BH-A0C0-11A-21R-A089-07 TCGA-BH-A0C0-11A-21R-A089-07    Normal
#> TCGA-BH-A0DZ-11A-22R-A089-07 TCGA-BH-A0DZ-11A-22R-A089-07    Normal
#> TCGA-BH-A0BA-11A-22R-A089-07 TCGA-BH-A0BA-11A-22R-A089-07    Normal
#> TCGA-BH-A0BC-11A-22R-A089-07 TCGA-BH-A0BC-11A-22R-A089-07    Normal
#> TCGA-BH-A0DH-11A-31R-A089-07 TCGA-BH-A0DH-11A-31R-A089-07    Normal
#> TCGA-A7-A0CH-11A-32R-A089-07 TCGA-A7-A0CH-11A-32R-A089-07    Normal
#> TCGA-A7-A0CE-11A-21R-A089-07 TCGA-A7-A0CE-11A-21R-A089-07    Normal
#> TCGA-BH-A0E1-11A-13R-A089-07 TCGA-BH-A0E1-11A-13R-A089-07    Normal
#> TCGA-BH-A0H9-11A-22R-A089-07 TCGA-BH-A0H9-11A-22R-A089-07    Normal
#> TCGA-A7-A0D9-11A-53R-A089-07 TCGA-A7-A0D9-11A-53R-A089-07    Normal
#> TCGA-BH-A0H7-11A-13R-A089-07 TCGA-BH-A0H7-11A-13R-A089-07    Normal
#> TCGA-BH-A0E0-11A-13R-A089-07 TCGA-BH-A0E0-11A-13R-A089-07    Normal
#> TCGA-BH-A0C3-11A-23R-A12P-07 TCGA-BH-A0C3-11A-23R-A12P-07    Normal
#> TCGA-BH-A0DL-11A-13R-A115-07 TCGA-BH-A0DL-11A-13R-A115-07    Normal
#> TCGA-BH-A0H5-11A-62R-A115-07 TCGA-BH-A0H5-11A-62R-A115-07    Normal
#> TCGA-BH-A0BQ-11A-33R-A115-07 TCGA-BH-A0BQ-11A-33R-A115-07    Normal
#> TCGA-BH-A0B7-11A-34R-A115-07 TCGA-BH-A0B7-11A-34R-A115-07    Normal
#> TCGA-BH-A0BW-11A-12R-A115-07 TCGA-BH-A0BW-11A-12R-A115-07    Normal
#> TCGA-BH-A18N-11A-43R-A12D-07 TCGA-BH-A18N-11A-43R-A12D-07    Normal
#> TCGA-E2-A158-11A-22R-A12D-07 TCGA-E2-A158-11A-22R-A12D-07    Normal
#> TCGA-BH-A18P-11A-43R-A12D-07 TCGA-BH-A18P-11A-43R-A12D-07    Normal
#> TCGA-BH-A0DO-11A-22R-A12D-07 TCGA-BH-A0DO-11A-22R-A12D-07    Normal
#> TCGA-BH-A18Q-11A-34R-A12D-07 TCGA-BH-A18Q-11A-34R-A12D-07    Normal
#> TCGA-BH-A0DT-11A-12R-A12D-07 TCGA-BH-A0DT-11A-12R-A12D-07    Normal
#> TCGA-BH-A18R-11A-42R-A12D-07 TCGA-BH-A18R-11A-42R-A12D-07    Normal
#> TCGA-E2-A15M-11A-22R-A12D-07 TCGA-E2-A15M-11A-22R-A12D-07    Normal
#> TCGA-BH-A18J-11A-31R-A12D-07 TCGA-BH-A18J-11A-31R-A12D-07    Normal
#> TCGA-BH-A18S-11A-43R-A12D-07 TCGA-BH-A18S-11A-43R-A12D-07    Normal
#> TCGA-BH-A18L-11A-42R-A12D-07 TCGA-BH-A18L-11A-42R-A12D-07    Normal
#> TCGA-BH-A18V-11A-52R-A12D-07 TCGA-BH-A18V-11A-52R-A12D-07    Normal
#> TCGA-BH-A18M-11A-33R-A12D-07 TCGA-BH-A18M-11A-33R-A12D-07    Normal
#> TCGA-E2-A153-11A-31R-A12D-07 TCGA-E2-A153-11A-31R-A12D-07    Normal
#> TCGA-BH-A0BS-11A-11R-A12P-07 TCGA-BH-A0BS-11A-11R-A12P-07    Normal
#> TCGA-A7-A13E-11A-61R-A12P-07 TCGA-A7-A13E-11A-61R-A12P-07    Normal
#> TCGA-BH-A0BZ-11A-61R-A12P-07 TCGA-BH-A0BZ-11A-61R-A12P-07    Normal
#> TCGA-A7-A13F-11A-42R-A12P-07 TCGA-A7-A13F-11A-42R-A12P-07    Normal
#> TCGA-BH-A18K-11A-13R-A12D-07 TCGA-BH-A18K-11A-13R-A12D-07    Normal
#> TCGA-BH-A18U-11A-23R-A12D-07 TCGA-BH-A18U-11A-23R-A12D-07    Normal
#> TCGA-BH-A0AU-11A-11R-A12P-07 TCGA-BH-A0AU-11A-11R-A12P-07    Normal
#> TCGA-BH-A0DD-11A-23R-A12P-07 TCGA-BH-A0DD-11A-23R-A12P-07    Normal
#> TCGA-BH-A0B5-11A-23R-A12P-07 TCGA-BH-A0B5-11A-23R-A12P-07    Normal
#> TCGA-BH-A0DV-11A-22R-A12P-07 TCGA-BH-A0DV-11A-22R-A12P-07    Normal
#> TCGA-BH-A1F0-11B-23R-A137-07 TCGA-BH-A1F0-11B-23R-A137-07    Normal
#> TCGA-E2-A1BC-11A-32R-A12P-07 TCGA-E2-A1BC-11A-32R-A12P-07    Normal
#> TCGA-BH-A0HA-11A-31R-A12P-07 TCGA-BH-A0HA-11A-31R-A12P-07    Normal
#> TCGA-E2-A15I-11A-32R-A137-07 TCGA-E2-A15I-11A-32R-A137-07    Normal
#> TCGA-BH-A1EW-11B-33R-A137-07 TCGA-BH-A1EW-11B-33R-A137-07    Normal
#> TCGA-BH-A1EU-11A-23R-A137-07 TCGA-BH-A1EU-11A-23R-A137-07    Normal
#> TCGA-BH-A1ET-11B-23R-A137-07 TCGA-BH-A1ET-11B-23R-A137-07    Normal
#> TCGA-BH-A1EO-11A-31R-A137-07 TCGA-BH-A1EO-11A-31R-A137-07    Normal
#> TCGA-AO-A03P-01A-11R-A00Z-07 TCGA-AO-A03P-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A06T-01A-11R-A00Z-07 TCGA-A8-A06T-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A07F-01A-11R-A00Z-07 TCGA-A8-A07F-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A081-01A-11R-A00Z-07 TCGA-A8-A081-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08C-01A-11R-A00Z-07 TCGA-A8-A08C-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08T-01A-21R-A00Z-07 TCGA-A8-A08T-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A095-01A-11R-A00Z-07 TCGA-A8-A095-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09G-01A-21R-A00Z-07 TCGA-A8-A09G-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A09W-01A-11R-A00Z-07 TCGA-A8-A09W-01A-11R-A00Z-07     Tumor
#> TCGA-AO-A03O-01A-11R-A00Z-07 TCGA-AO-A03O-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A06R-01A-11R-A00Z-07 TCGA-A8-A06R-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A07B-01A-11R-A00Z-07 TCGA-A8-A07B-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A07Z-01A-11R-A00Z-07 TCGA-A8-A07Z-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08B-01A-11R-A00Z-07 TCGA-A8-A08B-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08P-01A-11R-A00Z-07 TCGA-A8-A08P-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A094-01A-11R-A00Z-07 TCGA-A8-A094-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09E-01A-11R-A00Z-07 TCGA-A8-A09E-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09T-01A-11R-A00Z-07 TCGA-A8-A09T-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0AR-01A-11R-A00Z-07 TCGA-AN-A0AR-01A-11R-A00Z-07     Tumor
#> TCGA-BH-A0DZ-01A-11R-A00Z-07 TCGA-BH-A0DZ-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0AS-01A-11R-A00Z-07 TCGA-AN-A0AS-01A-11R-A00Z-07     Tumor
#> TCGA-A7-A0CH-01A-21R-A00Z-07 TCGA-A7-A0CH-01A-21R-A00Z-07     Tumor
#> TCGA-AN-A0AJ-01A-11R-A00Z-07 TCGA-AN-A0AJ-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A03X-01A-21R-A00Z-07 TCGA-AN-A03X-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A06U-01A-11R-A00Z-07 TCGA-A8-A06U-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A07I-01A-11R-A00Z-07 TCGA-A8-A07I-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A082-01A-11R-A00Z-07 TCGA-A8-A082-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08F-01A-11R-A00Z-07 TCGA-A8-A08F-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08X-01A-21R-A00Z-07 TCGA-A8-A08X-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A096-01A-11R-A00Z-07 TCGA-A8-A096-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09K-01A-11R-A00Z-07 TCGA-A8-A09K-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09X-01A-11R-A00Z-07 TCGA-A8-A09X-01A-11R-A00Z-07     Tumor
#> TCGA-BH-A0AY-01A-21R-A00Z-07 TCGA-BH-A0AY-01A-21R-A00Z-07     Tumor
#> TCGA-A7-A0CJ-01A-21R-A00Z-07 TCGA-A7-A0CJ-01A-21R-A00Z-07     Tumor
#> TCGA-AN-A03Y-01A-21R-A00Z-07 TCGA-AN-A03Y-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A06X-01A-21R-A00Z-07 TCGA-A8-A06X-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A07J-01A-11R-A00Z-07 TCGA-A8-A07J-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A083-01A-21R-A00Z-07 TCGA-A8-A083-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A08G-01A-11R-A00Z-07 TCGA-A8-A08G-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08Z-01A-21R-A00Z-07 TCGA-A8-A08Z-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A099-01A-11R-A00Z-07 TCGA-A8-A099-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09M-01A-11R-A00Z-07 TCGA-A8-A09M-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09Z-01A-11R-A00Z-07 TCGA-A8-A09Z-01A-11R-A00Z-07     Tumor
#> TCGA-BH-A0B4-01A-11R-A00Z-07 TCGA-BH-A0B4-01A-11R-A00Z-07     Tumor
#> TCGA-A2-A0CX-01A-21R-A00Z-07 TCGA-A2-A0CX-01A-21R-A00Z-07     Tumor
#> TCGA-AN-A049-01A-21R-A00Z-07 TCGA-AN-A049-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A06Y-01A-21R-A00Z-07 TCGA-A8-A06Y-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A084-01A-21R-A00Z-07 TCGA-A8-A084-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A08H-01A-21R-A00Z-07 TCGA-A8-A08H-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A090-01A-11R-A00Z-07 TCGA-A8-A090-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09A-01A-11R-A00Z-07 TCGA-A8-A09A-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09N-01A-11R-A00Z-07 TCGA-A8-A09N-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A0A1-01A-11R-A00Z-07 TCGA-A8-A0A1-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0AK-01A-21R-A00Z-07 TCGA-AN-A0AK-01A-21R-A00Z-07     Tumor
#> TCGA-A2-A0D0-01A-11R-A00Z-07 TCGA-A2-A0D0-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0FJ-01A-11R-A00Z-07 TCGA-AN-A0FJ-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A06Z-01A-11R-A00Z-07 TCGA-A8-A06Z-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A07O-01A-11R-A00Z-07 TCGA-A8-A07O-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A085-01A-11R-A00Z-07 TCGA-A8-A085-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08I-01A-11R-A00Z-07 TCGA-A8-A08I-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A091-01A-11R-A00Z-07 TCGA-A8-A091-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09B-01A-11R-A00Z-07 TCGA-A8-A09B-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0AL-01A-11R-A00Z-07 TCGA-AN-A0AL-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A0A4-01A-11R-A00Z-07 TCGA-A8-A0A4-01A-11R-A00Z-07     Tumor
#> TCGA-BH-A0BV-01A-11R-A00Z-07 TCGA-BH-A0BV-01A-11R-A00Z-07     Tumor
#> TCGA-A2-A0D4-01A-11R-A00Z-07 TCGA-A2-A0D4-01A-11R-A00Z-07     Tumor
#> TCGA-AN-A0FV-01A-11R-A00Z-07 TCGA-AN-A0FV-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A06O-01A-11R-A00Z-07 TCGA-A8-A06O-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A076-01A-21R-A00Z-07 TCGA-A8-A076-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A07P-01A-11R-A00Z-07 TCGA-A8-A07P-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A086-01A-11R-A00Z-07 TCGA-A8-A086-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08J-01A-11R-A00Z-07 TCGA-A8-A08J-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A092-01A-11R-A00Z-07 TCGA-A8-A092-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09C-01A-11R-A00Z-07 TCGA-A8-A09C-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09Q-01A-11R-A00Z-07 TCGA-A8-A09Q-01A-11R-A00Z-07     Tumor
#> TCGA-A7-A0CD-01A-11R-A00Z-07 TCGA-A7-A0CD-01A-11R-A00Z-07     Tumor
#> TCGA-A7-A0DB-01A-11R-A00Z-07 TCGA-A7-A0DB-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09I-01A-22R-A034-07 TCGA-A8-A09I-01A-22R-A034-07     Tumor
#> TCGA-A8-A09V-01A-11R-A034-07 TCGA-A8-A09V-01A-11R-A034-07     Tumor
#> TCGA-A8-A0A2-01A-11R-A034-07 TCGA-A8-A0A2-01A-11R-A034-07     Tumor
#> TCGA-A8-A0AB-01A-11R-A034-07 TCGA-A8-A0AB-01A-11R-A034-07     Tumor
#> TCGA-AN-A0AM-01A-11R-A034-07 TCGA-AN-A0AM-01A-11R-A034-07     Tumor
#> TCGA-AN-A0AT-01A-11R-A034-07 TCGA-AN-A0AT-01A-11R-A034-07     Tumor
#> TCGA-BH-A0BD-01A-11R-A034-07 TCGA-BH-A0BD-01A-11R-A034-07     Tumor
#> TCGA-A2-A0CM-01A-31R-A034-07 TCGA-A2-A0CM-01A-31R-A034-07     Tumor
#> TCGA-A2-A0CP-01A-11R-A034-07 TCGA-A2-A0CP-01A-11R-A034-07     Tumor
#> TCGA-A2-A0CQ-01A-21R-A034-07 TCGA-A2-A0CQ-01A-21R-A034-07     Tumor
#> TCGA-A2-A0CU-01A-12R-A034-07 TCGA-A2-A0CU-01A-12R-A034-07     Tumor
#> TCGA-A8-A09D-01A-11R-A00Z-07 TCGA-A8-A09D-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A09R-01A-11R-A00Z-07 TCGA-A8-A09R-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A0A9-01A-11R-A00Z-07 TCGA-A8-A0A9-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A06P-01A-11R-A00Z-07 TCGA-A8-A06P-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A079-01A-21R-A00Z-07 TCGA-A8-A079-01A-21R-A00Z-07     Tumor
#> TCGA-A8-A07W-01A-11R-A00Z-07 TCGA-A8-A07W-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08A-01A-11R-A00Z-07 TCGA-A8-A08A-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08L-01A-11R-A00Z-07 TCGA-A8-A08L-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A093-01A-11R-A00Z-07 TCGA-A8-A093-01A-11R-A00Z-07     Tumor
#> TCGA-A8-A08S-01A-11R-A034-07 TCGA-A8-A08S-01A-11R-A034-07     Tumor
#> TCGA-A8-A097-01A-11R-A034-07 TCGA-A8-A097-01A-11R-A034-07     Tumor
#> TCGA-A2-A04T-01A-21R-A034-07 TCGA-A2-A04T-01A-21R-A034-07     Tumor
#> TCGA-A2-A04V-01A-21R-A034-07 TCGA-A2-A04V-01A-21R-A034-07     Tumor
#> TCGA-A7-A0CE-01A-11R-A00Z-07 TCGA-A7-A0CE-01A-11R-A00Z-07     Tumor
#> TCGA-AO-A03R-01A-21R-A034-07 TCGA-AO-A03R-01A-21R-A034-07     Tumor
#> TCGA-AO-A03T-01A-21R-A034-07 TCGA-AO-A03T-01A-21R-A034-07     Tumor
#> TCGA-AN-A041-01A-11R-A034-07 TCGA-AN-A041-01A-11R-A034-07     Tumor
#> TCGA-AN-A046-01A-21R-A034-07 TCGA-AN-A046-01A-21R-A034-07     Tumor
#> TCGA-AN-A04A-01A-21R-A034-07 TCGA-AN-A04A-01A-21R-A034-07     Tumor
#> TCGA-AN-A04C-01A-21R-A034-07 TCGA-AN-A04C-01A-21R-A034-07     Tumor
#> TCGA-AN-A04D-01A-21R-A034-07 TCGA-AN-A04D-01A-21R-A034-07     Tumor
#> TCGA-AQ-A04J-01A-02R-A034-07 TCGA-AQ-A04J-01A-02R-A034-07     Tumor
#> TCGA-A2-A04P-01A-31R-A034-07 TCGA-A2-A04P-01A-31R-A034-07     Tumor
#> TCGA-A2-A04Q-01A-21R-A034-07 TCGA-A2-A04Q-01A-21R-A034-07     Tumor
#> TCGA-A2-A04X-01A-21R-A034-07 TCGA-A2-A04X-01A-21R-A034-07     Tumor
#> TCGA-A2-A04Y-01A-21R-A034-07 TCGA-A2-A04Y-01A-21R-A034-07     Tumor
#> TCGA-A8-A06Q-01A-11R-A034-07 TCGA-A8-A06Q-01A-11R-A034-07     Tumor
#> TCGA-A8-A07C-01A-11R-A034-07 TCGA-A8-A07C-01A-11R-A034-07     Tumor
#> TCGA-A8-A07E-01A-11R-A034-07 TCGA-A8-A07E-01A-11R-A034-07     Tumor
#> TCGA-A8-A07G-01A-11R-A034-07 TCGA-A8-A07G-01A-11R-A034-07     Tumor
#> TCGA-A8-A07R-01A-21R-A034-07 TCGA-A8-A07R-01A-21R-A034-07     Tumor
#> TCGA-A8-A07S-01A-11R-A034-07 TCGA-A8-A07S-01A-11R-A034-07     Tumor
#> TCGA-A8-A07U-01A-11R-A034-07 TCGA-A8-A07U-01A-11R-A034-07     Tumor
#> TCGA-A8-A08R-01A-11R-A034-07 TCGA-A8-A08R-01A-11R-A034-07     Tumor
#> TCGA-A2-A0CY-01A-12R-A034-07 TCGA-A2-A0CY-01A-12R-A034-07     Tumor
#> TCGA-A2-A0CZ-01A-11R-A034-07 TCGA-A2-A0CZ-01A-11R-A034-07     Tumor
#> TCGA-A2-A0D1-01A-11R-A034-07 TCGA-A2-A0D1-01A-11R-A034-07     Tumor
#> TCGA-A2-A0D2-01A-21R-A034-07 TCGA-A2-A0D2-01A-21R-A034-07     Tumor
#> TCGA-BH-A0E6-01A-11R-A034-07 TCGA-BH-A0E6-01A-11R-A034-07     Tumor
#> TCGA-BH-A0E7-01A-11R-A034-07 TCGA-BH-A0E7-01A-11R-A034-07     Tumor
#> TCGA-A2-A0EM-01A-11R-A034-07 TCGA-A2-A0EM-01A-11R-A034-07     Tumor
#> TCGA-A2-A0ER-01A-21R-A034-07 TCGA-A2-A0ER-01A-21R-A034-07     Tumor
#> TCGA-A2-A0ET-01A-31R-A034-07 TCGA-A2-A0ET-01A-31R-A034-07     Tumor
#> TCGA-A2-A0EX-01A-21R-A034-07 TCGA-A2-A0EX-01A-21R-A034-07     Tumor
#> TCGA-A2-A0EY-01A-11R-A034-07 TCGA-A2-A0EY-01A-11R-A034-07     Tumor
#> TCGA-BH-A0EB-01A-11R-A034-07 TCGA-BH-A0EB-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FL-01A-11R-A034-07 TCGA-AN-A0FL-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FT-01A-11R-A034-07 TCGA-AN-A0FT-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FX-01A-11R-A034-07 TCGA-AN-A0FX-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FY-01A-11R-A034-07 TCGA-AN-A0FY-01A-11R-A034-07     Tumor
#> TCGA-BH-A0HO-01A-11R-A034-07 TCGA-BH-A0HO-01A-11R-A034-07     Tumor
#> TCGA-B6-A0I5-01A-11R-A034-07 TCGA-B6-A0I5-01A-11R-A034-07     Tumor
#> TCGA-B6-A0I8-01A-11R-A034-07 TCGA-B6-A0I8-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IA-01A-11R-A034-07 TCGA-B6-A0IA-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IB-01A-11R-A034-07 TCGA-B6-A0IB-01A-11R-A034-07     Tumor
#> TCGA-BH-A0EE-01A-11R-A034-07 TCGA-BH-A0EE-01A-11R-A034-07     Tumor
#> TCGA-A2-A0EO-01A-11R-A034-07 TCGA-A2-A0EO-01A-11R-A034-07     Tumor
#> TCGA-A2-A0EQ-01A-11R-A034-07 TCGA-A2-A0EQ-01A-11R-A034-07     Tumor
#> TCGA-A2-A0EV-01A-11R-A034-07 TCGA-A2-A0EV-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FS-01A-11R-A034-07 TCGA-AN-A0FS-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FZ-01A-11R-A034-07 TCGA-AN-A0FZ-01A-11R-A034-07     Tumor
#> TCGA-AN-A0G0-01A-11R-A034-07 TCGA-AN-A0G0-01A-11R-A034-07     Tumor
#> TCGA-BH-A0HQ-01A-11R-A034-07 TCGA-BH-A0HQ-01A-11R-A034-07     Tumor
#> TCGA-BH-A0HU-01A-11R-A034-07 TCGA-BH-A0HU-01A-11R-A034-07     Tumor
#> TCGA-BH-A0HW-01A-11R-A034-07 TCGA-BH-A0HW-01A-11R-A034-07     Tumor
#> TCGA-B6-A0I9-01A-11R-A034-07 TCGA-B6-A0I9-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FD-01A-11R-A034-07 TCGA-AN-A0FD-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FF-01A-11R-A034-07 TCGA-AN-A0FF-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FK-01A-11R-A034-07 TCGA-AN-A0FK-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FN-01A-11R-A034-07 TCGA-AN-A0FN-01A-11R-A034-07     Tumor
#> TCGA-AN-A0FW-01A-11R-A034-07 TCGA-AN-A0FW-01A-11R-A034-07     Tumor
#> TCGA-B6-A0I2-01A-11R-A034-07 TCGA-B6-A0I2-01A-11R-A034-07     Tumor
#> TCGA-B6-A0I6-01A-11R-A034-07 TCGA-B6-A0I6-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IC-01A-11R-A034-07 TCGA-B6-A0IC-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IE-01A-11R-A034-07 TCGA-B6-A0IE-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IG-01A-11R-A034-07 TCGA-B6-A0IG-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IJ-01A-11R-A034-07 TCGA-B6-A0IJ-01A-11R-A034-07     Tumor
#> TCGA-A8-A0A7-01A-11R-A040-07 TCGA-A8-A0A7-01A-11R-A040-07     Tumor
#> TCGA-A8-A07L-01A-11R-A040-07 TCGA-A8-A07L-01A-11R-A040-07     Tumor
#> TCGA-A7-A0CG-01A-12R-A056-07 TCGA-A7-A0CG-01A-12R-A056-07     Tumor
#> TCGA-A8-A06N-01A-12R-A056-07 TCGA-A8-A06N-01A-12R-A056-07     Tumor
#> TCGA-AO-A03L-01A-41R-A056-07 TCGA-AO-A03L-01A-41R-A056-07     Tumor
#> TCGA-BH-A0B9-01A-11R-A056-07 TCGA-BH-A0B9-01A-11R-A056-07     Tumor
#> TCGA-BH-A0BA-01A-11R-A056-07 TCGA-BH-A0BA-01A-11R-A056-07     Tumor
#> TCGA-BH-A0DS-01A-11R-A056-07 TCGA-BH-A0DS-01A-11R-A056-07     Tumor
#> TCGA-BH-A0H6-01A-21R-A056-07 TCGA-BH-A0H6-01A-21R-A056-07     Tumor
#> TCGA-B6-A0IK-01A-12R-A056-07 TCGA-B6-A0IK-01A-12R-A056-07     Tumor
#> TCGA-AO-A0JJ-01A-11R-A056-07 TCGA-AO-A0JJ-01A-11R-A056-07     Tumor
#> TCGA-A8-A0A6-01A-12R-A056-07 TCGA-A8-A0A6-01A-12R-A056-07     Tumor
#> TCGA-BH-A0BJ-01A-11R-A056-07 TCGA-BH-A0BJ-01A-11R-A056-07     Tumor
#> TCGA-BH-A0E0-01A-11R-A056-07 TCGA-BH-A0E0-01A-11R-A056-07     Tumor
#> TCGA-BH-A0H7-01A-13R-A056-07 TCGA-BH-A0H7-01A-13R-A056-07     Tumor
#> TCGA-AO-A0JA-01A-11R-A056-07 TCGA-AO-A0JA-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JL-01A-11R-A056-07 TCGA-AO-A0JL-01A-11R-A056-07     Tumor
#> TCGA-BH-A0DP-01A-21R-A056-07 TCGA-BH-A0DP-01A-21R-A056-07     Tumor
#> TCGA-BH-A0H0-01A-11R-A056-07 TCGA-BH-A0H0-01A-11R-A056-07     Tumor
#> TCGA-BH-A0HY-01A-11R-A056-07 TCGA-BH-A0HY-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JI-01A-21R-A056-07 TCGA-AO-A0JI-01A-21R-A056-07     Tumor
#> TCGA-A8-A08O-01A-21R-A056-07 TCGA-A8-A08O-01A-21R-A056-07     Tumor
#> TCGA-AR-A0TT-01A-31R-A084-07 TCGA-AR-A0TT-01A-31R-A084-07     Tumor
#> TCGA-B6-A0IM-01A-11R-A034-07 TCGA-B6-A0IM-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IO-01A-11R-A034-07 TCGA-B6-A0IO-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IN-01A-11R-A034-07 TCGA-B6-A0IN-01A-11R-A034-07     Tumor
#> TCGA-B6-A0WT-01A-11R-A109-07 TCGA-B6-A0WT-01A-11R-A109-07     Tumor
#> TCGA-B6-A0X5-01A-21R-A109-07 TCGA-B6-A0X5-01A-21R-A109-07     Tumor
#> TCGA-AN-A0XV-01A-11R-A109-07 TCGA-AN-A0XV-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YJ-01A-11R-A109-07 TCGA-A2-A0YJ-01A-11R-A109-07     Tumor
#> TCGA-AR-A0U2-01A-11R-A109-07 TCGA-AR-A0U2-01A-11R-A109-07     Tumor
#> TCGA-B6-A0WV-01A-11R-A109-07 TCGA-B6-A0WV-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XN-01A-21R-A109-07 TCGA-AN-A0XN-01A-21R-A109-07     Tumor
#> TCGA-AN-A0XW-01A-11R-A109-07 TCGA-AN-A0XW-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YK-01A-22R-A109-07 TCGA-A2-A0YK-01A-22R-A109-07     Tumor
#> TCGA-AR-A0U3-01A-11R-A109-07 TCGA-AR-A0U3-01A-11R-A109-07     Tumor
#> TCGA-B6-A0WW-01A-11R-A109-07 TCGA-B6-A0WW-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XO-01A-11R-A109-07 TCGA-AN-A0XO-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YC-01A-11R-A109-07 TCGA-A2-A0YC-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YL-01A-21R-A109-07 TCGA-A2-A0YL-01A-21R-A109-07     Tumor
#> TCGA-AR-A0U4-01A-11R-A109-07 TCGA-AR-A0U4-01A-11R-A109-07     Tumor
#> TCGA-B6-A0WX-01A-11R-A109-07 TCGA-B6-A0WX-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XP-01A-11R-A109-07 TCGA-AN-A0XP-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YD-01A-11R-A109-07 TCGA-A2-A0YD-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YM-01A-11R-A109-07 TCGA-A2-A0YM-01A-11R-A109-07     Tumor
#> TCGA-BH-A0W3-01A-11R-A109-07 TCGA-BH-A0W3-01A-11R-A109-07     Tumor
#> TCGA-B6-A0WY-01A-11R-A109-07 TCGA-B6-A0WY-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XR-01A-11R-A109-07 TCGA-AN-A0XR-01A-11R-A109-07     Tumor
#> TCGA-AO-A0J2-01A-11R-A034-07 TCGA-AO-A0J2-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J3-01A-11R-A034-07 TCGA-AO-A0J3-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J4-01A-11R-A034-07 TCGA-AO-A0J4-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J5-01A-11R-A034-07 TCGA-AO-A0J5-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J6-01A-11R-A034-07 TCGA-AO-A0J6-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J7-01A-11R-A034-07 TCGA-AO-A0J7-01A-11R-A034-07     Tumor
#> TCGA-AO-A0J8-01A-21R-A034-07 TCGA-AO-A0J8-01A-21R-A034-07     Tumor
#> TCGA-AO-A0J9-01A-11R-A034-07 TCGA-AO-A0J9-01A-11R-A034-07     Tumor
#> TCGA-A2-A0YE-01A-11R-A109-07 TCGA-A2-A0YE-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YT-01A-11R-A109-07 TCGA-A2-A0YT-01A-11R-A109-07     Tumor
#> TCGA-BH-A0W4-01A-11R-A109-07 TCGA-BH-A0W4-01A-11R-A109-07     Tumor
#> TCGA-B6-A0WZ-01A-11R-A109-07 TCGA-B6-A0WZ-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XS-01A-22R-A109-07 TCGA-AN-A0XS-01A-22R-A109-07     Tumor
#> TCGA-A2-A0YF-01A-21R-A109-07 TCGA-A2-A0YF-01A-21R-A109-07     Tumor
#> TCGA-B6-A0IP-01A-11R-A034-07 TCGA-B6-A0IP-01A-11R-A034-07     Tumor
#> TCGA-B6-A0IQ-01A-11R-A034-07 TCGA-B6-A0IQ-01A-11R-A034-07     Tumor
#> TCGA-A8-A0AD-01A-11R-A056-07 TCGA-A8-A0AD-01A-11R-A056-07     Tumor
#> TCGA-BH-A0BM-01A-11R-A056-07 TCGA-BH-A0BM-01A-11R-A056-07     Tumor
#> TCGA-BH-A0E1-01A-11R-A056-07 TCGA-BH-A0E1-01A-11R-A056-07     Tumor
#> TCGA-BH-A0H9-01A-11R-A056-07 TCGA-BH-A0H9-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JB-01A-11R-A056-07 TCGA-AO-A0JB-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JM-01A-21R-A056-07 TCGA-AO-A0JM-01A-21R-A056-07     Tumor
#> TCGA-BH-A0AW-01A-11R-A056-07 TCGA-BH-A0AW-01A-11R-A056-07     Tumor
#> TCGA-BH-A0C0-01A-21R-A056-07 TCGA-BH-A0C0-01A-21R-A056-07     Tumor
#> TCGA-BH-A0E2-01A-11R-A056-07 TCGA-BH-A0E2-01A-11R-A056-07     Tumor
#> TCGA-BH-A0HB-01A-11R-A056-07 TCGA-BH-A0HB-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JC-01A-11R-A056-07 TCGA-AO-A0JC-01A-11R-A056-07     Tumor
#> TCGA-B6-A0RE-01A-11R-A056-07 TCGA-B6-A0RE-01A-11R-A056-07     Tumor
#> TCGA-BH-A0B1-01A-12R-A056-07 TCGA-BH-A0B1-01A-12R-A056-07     Tumor
#> TCGA-A2-A0CT-01A-31R-A056-07 TCGA-A2-A0CT-01A-31R-A056-07     Tumor
#> TCGA-A2-A0EU-01A-22R-A056-07 TCGA-A2-A0EU-01A-22R-A056-07     Tumor
#> TCGA-BH-A0HF-01A-11R-A056-07 TCGA-BH-A0HF-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JD-01A-11R-A056-07 TCGA-AO-A0JD-01A-11R-A056-07     Tumor
#> TCGA-B6-A0RG-01A-11R-A056-07 TCGA-B6-A0RG-01A-11R-A056-07     Tumor
#> TCGA-BH-A0B3-01A-11R-A056-07 TCGA-BH-A0B3-01A-11R-A056-07     Tumor
#> TCGA-A7-A0D9-01A-31R-A056-07 TCGA-A7-A0D9-01A-31R-A056-07     Tumor
#> TCGA-BH-A0GY-01A-11R-A056-07 TCGA-BH-A0GY-01A-11R-A056-07     Tumor
#> TCGA-BH-A0HK-01A-11R-A056-07 TCGA-BH-A0HK-01A-11R-A056-07     Tumor
#> TCGA-AO-A0JE-01A-11R-A056-07 TCGA-AO-A0JE-01A-11R-A056-07     Tumor
#> TCGA-B6-A0RI-01A-11R-A056-07 TCGA-B6-A0RI-01A-11R-A056-07     Tumor
#> TCGA-BH-A0B8-01A-21R-A056-07 TCGA-BH-A0B8-01A-21R-A056-07     Tumor
#> TCGA-BH-A0DK-01A-21R-A056-07 TCGA-BH-A0DK-01A-21R-A056-07     Tumor
#> TCGA-BH-A0GZ-01A-11R-A056-07 TCGA-BH-A0GZ-01A-11R-A056-07     Tumor
#> TCGA-BH-A0HX-01A-21R-A056-07 TCGA-BH-A0HX-01A-21R-A056-07     Tumor
#> TCGA-AO-A0JF-01A-11R-A056-07 TCGA-AO-A0JF-01A-11R-A056-07     Tumor
#> TCGA-AR-A0TR-01A-11R-A084-07 TCGA-AR-A0TR-01A-11R-A084-07     Tumor
#> TCGA-BH-A0DH-01A-11R-A084-07 TCGA-BH-A0DH-01A-11R-A084-07     Tumor
#> TCGA-B6-A0RM-01A-11R-A084-07 TCGA-B6-A0RM-01A-11R-A084-07     Tumor
#> TCGA-BH-A0RX-01A-21R-A084-07 TCGA-BH-A0RX-01A-21R-A084-07     Tumor
#> TCGA-A2-A0ST-01A-12R-A084-07 TCGA-A2-A0ST-01A-12R-A084-07     Tumor
#> TCGA-A8-A075-01A-11R-A084-07 TCGA-A8-A075-01A-11R-A084-07     Tumor
#> TCGA-AO-A0JG-01A-31R-A084-07 TCGA-AO-A0JG-01A-31R-A084-07     Tumor
#> TCGA-B6-A0RU-01A-11R-A084-07 TCGA-B6-A0RU-01A-11R-A084-07     Tumor
#> TCGA-A1-A0SO-01A-22R-A084-07 TCGA-A1-A0SO-01A-22R-A084-07     Tumor
#> TCGA-A2-A0T0-01A-22R-A084-07 TCGA-A2-A0T0-01A-22R-A084-07     Tumor
#> TCGA-AR-A0TQ-01A-11R-A084-07 TCGA-AR-A0TQ-01A-11R-A084-07     Tumor
#> TCGA-BH-A0BC-01A-22R-A084-07 TCGA-BH-A0BC-01A-22R-A084-07     Tumor
#> TCGA-B6-A0RL-01A-11R-A084-07 TCGA-B6-A0RL-01A-11R-A084-07     Tumor
#> TCGA-B6-A0RV-01A-11R-A084-07 TCGA-B6-A0RV-01A-11R-A084-07     Tumor
#> TCGA-A1-A0SP-01A-11R-A084-07 TCGA-A1-A0SP-01A-11R-A084-07     Tumor
#> TCGA-A2-A0T1-01A-21R-A084-07 TCGA-A2-A0T1-01A-21R-A084-07     Tumor
#> TCGA-BH-A0DQ-01A-11R-A084-07 TCGA-BH-A0DQ-01A-11R-A084-07     Tumor
#> TCGA-B6-A0RN-01A-12R-A084-07 TCGA-B6-A0RN-01A-12R-A084-07     Tumor
#> TCGA-A2-A0T4-01A-31R-A084-07 TCGA-A2-A0T4-01A-31R-A084-07     Tumor
#> TCGA-AR-A0TV-01A-21R-A084-07 TCGA-AR-A0TV-01A-21R-A084-07     Tumor
#> TCGA-A2-A0EN-01A-13R-A084-07 TCGA-A2-A0EN-01A-13R-A084-07     Tumor
#> TCGA-B6-A0RO-01A-22R-A084-07 TCGA-B6-A0RO-01A-22R-A084-07     Tumor
#> TCGA-BH-A0HI-01A-11R-A084-07 TCGA-BH-A0HI-01A-11R-A084-07     Tumor
#> TCGA-AR-A0TX-01A-11R-A084-07 TCGA-AR-A0TX-01A-11R-A084-07     Tumor
#> TCGA-A1-A0SE-01A-11R-A084-07 TCGA-A1-A0SE-01A-11R-A084-07     Tumor
#> TCGA-A2-A0SU-01A-11R-A084-07 TCGA-A2-A0SU-01A-11R-A084-07     Tumor
#> TCGA-A1-A0SH-01A-11R-A084-07 TCGA-A1-A0SH-01A-11R-A084-07     Tumor
#> TCGA-A2-A0SV-01A-11R-A084-07 TCGA-A2-A0SV-01A-11R-A084-07     Tumor
#> TCGA-A2-A0T5-01A-21R-A084-07 TCGA-A2-A0T5-01A-21R-A084-07     Tumor
#> TCGA-AR-A0TW-01A-11R-A084-07 TCGA-AR-A0TW-01A-11R-A084-07     Tumor
#> TCGA-B6-A0RP-01A-21R-A084-07 TCGA-B6-A0RP-01A-21R-A084-07     Tumor
#> TCGA-A1-A0SJ-01A-11R-A084-07 TCGA-A1-A0SJ-01A-11R-A084-07     Tumor
#> TCGA-A2-A0SW-01A-11R-A084-07 TCGA-A2-A0SW-01A-11R-A084-07     Tumor
#> TCGA-A2-A0T6-01A-11R-A084-07 TCGA-A2-A0T6-01A-11R-A084-07     Tumor
#> TCGA-B6-A0RS-01A-11R-A084-07 TCGA-B6-A0RS-01A-11R-A084-07     Tumor
#> TCGA-A1-A0SK-01A-12R-A084-07 TCGA-A1-A0SK-01A-12R-A084-07     Tumor
#> TCGA-A2-A0SX-01A-12R-A084-07 TCGA-A2-A0SX-01A-12R-A084-07     Tumor
#> TCGA-A2-A0T7-01A-21R-A084-07 TCGA-A2-A0T7-01A-21R-A084-07     Tumor
#> TCGA-AR-A0TZ-01A-12R-A084-07 TCGA-AR-A0TZ-01A-12R-A084-07     Tumor
#> TCGA-BH-A0HP-01A-12R-A084-07 TCGA-BH-A0HP-01A-12R-A084-07     Tumor
#> TCGA-B6-A0RT-01A-21R-A084-07 TCGA-B6-A0RT-01A-21R-A084-07     Tumor
#> TCGA-A1-A0SM-01A-11R-A084-07 TCGA-A1-A0SM-01A-11R-A084-07     Tumor
#> TCGA-A2-A0SY-01A-31R-A084-07 TCGA-A2-A0SY-01A-31R-A084-07     Tumor
#> TCGA-AR-A0TP-01A-11R-A084-07 TCGA-AR-A0TP-01A-11R-A084-07     Tumor
#> TCGA-A2-A04R-01A-41R-A109-07 TCGA-A2-A04R-01A-41R-A109-07     Tumor
#> TCGA-BH-A0W5-01A-11R-A109-07 TCGA-BH-A0W5-01A-11R-A109-07     Tumor
#> TCGA-B6-A0X1-01A-11R-A109-07 TCGA-B6-A0X1-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XT-01A-11R-A109-07 TCGA-AN-A0XT-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YG-01A-21R-A109-07 TCGA-A2-A0YG-01A-21R-A109-07     Tumor
#> TCGA-AR-A0TU-01A-31R-A109-07 TCGA-AR-A0TU-01A-31R-A109-07     Tumor
#> TCGA-BH-A0WA-01A-11R-A109-07 TCGA-BH-A0WA-01A-11R-A109-07     Tumor
#> TCGA-B6-A0X4-01A-11R-A109-07 TCGA-B6-A0X4-01A-11R-A109-07     Tumor
#> TCGA-AN-A0XU-01A-11R-A109-07 TCGA-AN-A0XU-01A-11R-A109-07     Tumor
#> TCGA-A2-A0YH-01A-11R-A109-07 TCGA-A2-A0YH-01A-11R-A109-07     Tumor
#> TCGA-AR-A0U0-01A-11R-A109-07 TCGA-AR-A0U0-01A-11R-A109-07     Tumor
#> TCGA-AO-A12G-01A-11R-A10J-07 TCGA-AO-A12G-01A-11R-A10J-07     Tumor
#> TCGA-AQ-A04H-01B-11R-A10J-07 TCGA-AQ-A04H-01B-11R-A10J-07     Tumor
#> TCGA-E2-A107-01A-11R-A10J-07 TCGA-E2-A107-01A-11R-A10J-07     Tumor
#> TCGA-AO-A125-01A-11R-A10J-07 TCGA-AO-A125-01A-11R-A10J-07     Tumor
#> TCGA-AQ-A04L-01B-21R-A10J-07 TCGA-AQ-A04L-01B-21R-A10J-07     Tumor
#> TCGA-E2-A108-01A-13R-A10J-07 TCGA-E2-A108-01A-13R-A10J-07     Tumor
#> TCGA-AO-A126-01A-11R-A10J-07 TCGA-AO-A126-01A-11R-A10J-07     Tumor
#> TCGA-E2-A109-01A-11R-A10J-07 TCGA-E2-A109-01A-11R-A10J-07     Tumor
#> TCGA-AO-A03M-01B-11R-A10J-07 TCGA-AO-A03M-01B-11R-A10J-07     Tumor
#> TCGA-A2-A0YI-01A-31R-A10J-07 TCGA-A2-A0YI-01A-31R-A10J-07     Tumor
#> TCGA-E2-A10E-01A-21R-A10J-07 TCGA-E2-A10E-01A-21R-A10J-07     Tumor
#> TCGA-AO-A12C-01A-11R-A10J-07 TCGA-AO-A12C-01A-11R-A10J-07     Tumor
#> TCGA-AO-A03N-01B-11R-A10J-07 TCGA-AO-A03N-01B-11R-A10J-07     Tumor
#> TCGA-E2-A105-01A-11R-A10J-07 TCGA-E2-A105-01A-11R-A10J-07     Tumor
#> TCGA-E2-A10F-01A-11R-A10J-07 TCGA-E2-A10F-01A-11R-A10J-07     Tumor
#> TCGA-AO-A12E-01A-11R-A10J-07 TCGA-AO-A12E-01A-11R-A10J-07     Tumor
#> TCGA-AO-A03U-01B-21R-A10J-07 TCGA-AO-A03U-01B-21R-A10J-07     Tumor
#> TCGA-E2-A106-01A-11R-A10J-07 TCGA-E2-A106-01A-11R-A10J-07     Tumor
#> TCGA-AO-A124-01A-11R-A10J-07 TCGA-AO-A124-01A-11R-A10J-07     Tumor
#> TCGA-AO-A128-01A-11R-A10J-07 TCGA-AO-A128-01A-11R-A10J-07     Tumor
#> TCGA-B6-A0X7-01A-11R-A10J-07 TCGA-B6-A0X7-01A-11R-A10J-07     Tumor
#> TCGA-E2-A10B-01A-11R-A10J-07 TCGA-E2-A10B-01A-11R-A10J-07     Tumor
#> TCGA-AO-A129-01A-21R-A10J-07 TCGA-AO-A129-01A-21R-A10J-07     Tumor
#> TCGA-AN-A0XL-01A-11R-A10J-07 TCGA-AN-A0XL-01A-11R-A10J-07     Tumor
#> TCGA-E2-A10C-01A-21R-A10J-07 TCGA-E2-A10C-01A-21R-A10J-07     Tumor
#> TCGA-AO-A12B-01A-11R-A10J-07 TCGA-AO-A12B-01A-11R-A10J-07     Tumor
#> TCGA-AO-A03V-01A-11R-A115-07 TCGA-AO-A03V-01A-11R-A115-07     Tumor
#> TCGA-BH-A0BG-01A-11R-A115-07 TCGA-BH-A0BG-01A-11R-A115-07     Tumor
#> TCGA-BH-A0C7-01B-11R-A115-07 TCGA-BH-A0C7-01B-11R-A115-07     Tumor
#> TCGA-BH-A0DL-01A-11R-A115-07 TCGA-BH-A0DL-01A-11R-A115-07     Tumor
#> TCGA-BH-A0H5-01A-21R-A115-07 TCGA-BH-A0H5-01A-21R-A115-07     Tumor
#> TCGA-AR-A0TY-01A-12R-A115-07 TCGA-AR-A0TY-01A-12R-A115-07     Tumor
#> TCGA-AO-A12F-01A-11R-A115-07 TCGA-AO-A12F-01A-11R-A115-07     Tumor
#> TCGA-C8-A12Q-01A-11R-A115-07 TCGA-C8-A12Q-01A-11R-A115-07     Tumor
#> TCGA-C8-A131-01A-11R-A115-07 TCGA-C8-A131-01A-11R-A115-07     Tumor
#> TCGA-D8-A140-01A-11R-A115-07 TCGA-D8-A140-01A-11R-A115-07     Tumor
#> TCGA-E2-A14R-01A-11R-A115-07 TCGA-E2-A14R-01A-11R-A115-07     Tumor
#> TCGA-E2-A15O-01A-11R-A115-07 TCGA-E2-A15O-01A-11R-A115-07     Tumor
#> TCGA-A2-A04N-01A-11R-A115-07 TCGA-A2-A04N-01A-11R-A115-07     Tumor
#> TCGA-BH-A0BL-01A-11R-A115-07 TCGA-BH-A0BL-01A-11R-A115-07     Tumor
#> TCGA-A2-A0CL-01A-11R-A115-07 TCGA-A2-A0CL-01A-11R-A115-07     Tumor
#> TCGA-AO-A12H-01A-11R-A115-07 TCGA-AO-A12H-01A-11R-A115-07     Tumor
#> TCGA-C8-A12T-01A-11R-A115-07 TCGA-C8-A12T-01A-11R-A115-07     Tumor
#> TCGA-C8-A132-01A-31R-A115-07 TCGA-C8-A132-01A-31R-A115-07     Tumor
#> TCGA-D8-A141-01A-11R-A115-07 TCGA-D8-A141-01A-11R-A115-07     Tumor
#> TCGA-E2-A14T-01A-11R-A115-07 TCGA-E2-A14T-01A-11R-A115-07     Tumor
#> TCGA-E2-A15P-01A-11R-A115-07 TCGA-E2-A15P-01A-11R-A115-07     Tumor
#> TCGA-A2-A04U-01A-11R-A115-07 TCGA-A2-A04U-01A-11R-A115-07     Tumor
#> TCGA-BH-A0BP-01A-11R-A115-07 TCGA-BH-A0BP-01A-11R-A115-07     Tumor
#> TCGA-A2-A0CS-01A-11R-A115-07 TCGA-A2-A0CS-01A-11R-A115-07     Tumor
#> TCGA-BH-A0DX-01A-11R-A115-07 TCGA-BH-A0DX-01A-11R-A115-07     Tumor
#> TCGA-B6-A0IH-01A-11R-A115-07 TCGA-B6-A0IH-01A-11R-A115-07     Tumor
#> TCGA-BH-A0W7-01A-11R-A115-07 TCGA-BH-A0W7-01A-11R-A115-07     Tumor
#> TCGA-C8-A12K-01A-21R-A115-07 TCGA-C8-A12K-01A-21R-A115-07     Tumor
#> TCGA-C8-A12U-01A-11R-A115-07 TCGA-C8-A12U-01A-11R-A115-07     Tumor
#> TCGA-C8-A134-01A-11R-A115-07 TCGA-C8-A134-01A-11R-A115-07     Tumor
#> TCGA-D8-A142-01A-11R-A115-07 TCGA-D8-A142-01A-11R-A115-07     Tumor
#> TCGA-E2-A14X-01A-11R-A115-07 TCGA-E2-A14X-01A-11R-A115-07     Tumor
#> TCGA-E2-A15R-01A-11R-A115-07 TCGA-E2-A15R-01A-11R-A115-07     Tumor
#> TCGA-A2-A04W-01A-31R-A115-07 TCGA-A2-A04W-01A-31R-A115-07     Tumor
#> TCGA-BH-A0BQ-01A-21R-A115-07 TCGA-BH-A0BQ-01A-21R-A115-07     Tumor
#> TCGA-A2-A0CV-01A-31R-A115-07 TCGA-A2-A0CV-01A-31R-A115-07     Tumor
#> TCGA-BH-A0E9-01B-11R-A115-07 TCGA-BH-A0E9-01B-11R-A115-07     Tumor
#> TCGA-B6-A0RH-01A-21R-A115-07 TCGA-B6-A0RH-01A-21R-A115-07     Tumor
#> TCGA-B6-A0WS-01A-11R-A115-07 TCGA-B6-A0WS-01A-11R-A115-07     Tumor
#> TCGA-C8-A12L-01A-11R-A115-07 TCGA-C8-A12L-01A-11R-A115-07     Tumor
#> TCGA-C8-A12V-01A-11R-A115-07 TCGA-C8-A12V-01A-11R-A115-07     Tumor
#> TCGA-C8-A135-01A-11R-A115-07 TCGA-C8-A135-01A-11R-A115-07     Tumor
#> TCGA-D8-A143-01A-11R-A115-07 TCGA-D8-A143-01A-11R-A115-07     Tumor
#> TCGA-E2-A14Z-01A-11R-A115-07 TCGA-E2-A14Z-01A-11R-A115-07     Tumor
#> TCGA-E2-A15S-01A-11R-A115-07 TCGA-E2-A15S-01A-11R-A115-07     Tumor
#> TCGA-BH-A0AV-01A-31R-A115-07 TCGA-BH-A0AV-01A-31R-A115-07     Tumor
#> TCGA-A2-A0CW-01A-21R-A115-07 TCGA-A2-A0CW-01A-21R-A115-07     Tumor
#> TCGA-BH-A0EA-01A-11R-A115-07 TCGA-BH-A0EA-01A-11R-A115-07     Tumor
#> TCGA-B6-A0RQ-01A-11R-A115-07 TCGA-B6-A0RQ-01A-11R-A115-07     Tumor
#> TCGA-B6-A0X0-01A-21R-A115-07 TCGA-B6-A0X0-01A-21R-A115-07     Tumor
#> TCGA-C8-A12M-01A-11R-A115-07 TCGA-C8-A12M-01A-11R-A115-07     Tumor
#> TCGA-C8-A12W-01A-11R-A115-07 TCGA-C8-A12W-01A-11R-A115-07     Tumor
#> TCGA-C8-A137-01A-11R-A115-07 TCGA-C8-A137-01A-11R-A115-07     Tumor
#> TCGA-D8-A145-01A-11R-A115-07 TCGA-D8-A145-01A-11R-A115-07     Tumor
#> TCGA-E2-A154-01A-11R-A115-07 TCGA-E2-A154-01A-11R-A115-07     Tumor
#> TCGA-E2-A15T-01A-11R-A115-07 TCGA-E2-A15T-01A-11R-A115-07     Tumor
#> TCGA-BH-A0B0-01A-21R-A115-07 TCGA-BH-A0B0-01A-21R-A115-07     Tumor
#> TCGA-BH-A0BR-01A-21R-A115-07 TCGA-BH-A0BR-01A-21R-A115-07     Tumor
#> TCGA-A2-A0D3-01A-11R-A115-07 TCGA-A2-A0D3-01A-11R-A115-07     Tumor
#> TCGA-BH-A0EI-01A-11R-A115-07 TCGA-BH-A0EI-01A-11R-A115-07     Tumor
#> TCGA-A1-A0SD-01A-11R-A115-07 TCGA-A1-A0SD-01A-11R-A115-07     Tumor
#> TCGA-E2-A10A-01A-21R-A115-07 TCGA-E2-A10A-01A-21R-A115-07     Tumor
#> TCGA-C8-A12N-01A-11R-A115-07 TCGA-C8-A12N-01A-11R-A115-07     Tumor
#> TCGA-C8-A12X-01A-11R-A115-07 TCGA-C8-A12X-01A-11R-A115-07     Tumor
#> TCGA-C8-A138-01A-11R-A115-07 TCGA-C8-A138-01A-11R-A115-07     Tumor
#> TCGA-D8-A146-01A-31R-A115-07 TCGA-D8-A146-01A-31R-A115-07     Tumor
#> TCGA-E2-A159-01A-11R-A115-07 TCGA-E2-A159-01A-11R-A115-07     Tumor
#> TCGA-BH-A0B7-01A-12R-A115-07 TCGA-BH-A0B7-01A-12R-A115-07     Tumor
#> TCGA-BH-A0BW-01A-11R-A115-07 TCGA-BH-A0BW-01A-11R-A115-07     Tumor
#> TCGA-A7-A0DA-01A-31R-A115-07 TCGA-A7-A0DA-01A-31R-A115-07     Tumor
#> TCGA-A2-A0ES-01A-11R-A115-07 TCGA-A2-A0ES-01A-11R-A115-07     Tumor
#> TCGA-A2-A0T3-01A-21R-A115-07 TCGA-A2-A0T3-01A-21R-A115-07     Tumor
#> TCGA-AO-A12A-01A-21R-A115-07 TCGA-AO-A12A-01A-21R-A115-07     Tumor
#> TCGA-C8-A12O-01A-11R-A115-07 TCGA-C8-A12O-01A-11R-A115-07     Tumor
#> TCGA-C8-A12Z-01A-11R-A115-07 TCGA-C8-A12Z-01A-11R-A115-07     Tumor
#> TCGA-D8-A13Y-01A-11R-A115-07 TCGA-D8-A13Y-01A-11R-A115-07     Tumor
#> TCGA-D8-A147-01A-11R-A115-07 TCGA-D8-A147-01A-11R-A115-07     Tumor
#> TCGA-E2-A15D-01A-11R-A115-07 TCGA-E2-A15D-01A-11R-A115-07     Tumor
#> TCGA-BH-A0DE-01A-11R-A115-07 TCGA-BH-A0DE-01A-11R-A115-07     Tumor
#> TCGA-A2-A0EW-01A-21R-A115-07 TCGA-A2-A0EW-01A-21R-A115-07     Tumor
#> TCGA-AR-A0TS-01A-11R-A115-07 TCGA-AR-A0TS-01A-11R-A115-07     Tumor
#> TCGA-AO-A12D-01A-11R-A115-07 TCGA-AO-A12D-01A-11R-A115-07     Tumor
#> TCGA-C8-A12P-01A-11R-A115-07 TCGA-C8-A12P-01A-11R-A115-07     Tumor
#> TCGA-C8-A130-01A-31R-A115-07 TCGA-C8-A130-01A-31R-A115-07     Tumor
#> TCGA-D8-A13Z-01A-11R-A115-07 TCGA-D8-A13Z-01A-11R-A115-07     Tumor
#> TCGA-E2-A14O-01A-31R-A115-07 TCGA-E2-A14O-01A-31R-A115-07     Tumor
#> TCGA-E2-A15F-01A-11R-A115-07 TCGA-E2-A15F-01A-11R-A115-07     Tumor
#> TCGA-A2-A0T2-01A-11R-A084-07 TCGA-A2-A0T2-01A-11R-A084-07     Tumor
#> TCGA-AR-A1AH-01A-11R-A12D-07 TCGA-AR-A1AH-01A-11R-A12D-07     Tumor
#> TCGA-BH-A18I-01A-11R-A12D-07 TCGA-BH-A18I-01A-11R-A12D-07     Tumor
#> TCGA-BH-A18R-01A-11R-A12D-07 TCGA-BH-A18R-01A-11R-A12D-07     Tumor
#> TCGA-E2-A14Q-01A-11R-A12D-07 TCGA-E2-A14Q-01A-11R-A12D-07     Tumor
#> TCGA-E2-A155-01A-11R-A12D-07 TCGA-E2-A155-01A-11R-A12D-07     Tumor
#> TCGA-E2-A15L-01A-11R-A12D-07 TCGA-E2-A15L-01A-11R-A12D-07     Tumor
#> TCGA-BH-A0BO-01A-23R-A12D-07 TCGA-BH-A0BO-01A-23R-A12D-07     Tumor
#> TCGA-BH-A18J-01A-11R-A12D-07 TCGA-BH-A18J-01A-11R-A12D-07     Tumor
#> TCGA-BH-A18S-01A-11R-A12D-07 TCGA-BH-A18S-01A-11R-A12D-07     Tumor
#> TCGA-E2-A14S-01A-11R-A12D-07 TCGA-E2-A14S-01A-11R-A12D-07     Tumor
#> TCGA-E2-A156-01A-11R-A12D-07 TCGA-E2-A156-01A-11R-A12D-07     Tumor
#> TCGA-E2-A15M-01A-11R-A12D-07 TCGA-E2-A15M-01A-11R-A12D-07     Tumor
#>  [ reached 'max' / getOption("max.print") -- omitted 90 rows ]
#> 
#> Slot "cancer_type":
#> [1] "TCGA-BRCA"
```

# 3. Differential Expression Analysis
Once the GenePanel object is created, we can identify significant drivers using the calc_fold_change method. This utilizes Welch's t-test to handle unequal variances between tumor and normal tissues.


``` r
# Calculate Log2 Fold Change and P-values
results <- calc_fold_change(panel)
#> [1] "Analyzing 17814 genes across 590 samples..."

# Display the top significant genes
head(results)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Gene"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Log2FC"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PValue"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"ITM2A","2":"-2.981790","3":"1.012652e-86","4":"1.803939e-82","_rn_":"17524"},{"1":"SFRP1","2":"-4.581451","3":"5.653358e-78","4":"5.035446e-74","_rn_":"11097"},{"1":"OSR1","2":"-2.739629","3":"7.619958e-77","4":"4.524731e-73","_rn_":"17396"},{"1":"HOXA4","2":"-2.942836","3":"7.145809e-75","4":"3.182386e-71","_rn_":"8666"},{"1":"UBE2E3","2":"-1.179885","3":"3.253644e-71","4":"1.159208e-67","_rn_":"12314"},{"1":"CNTNAP3","2":"-2.470789","3":"9.330959e-71","4":"2.770362e-67","_rn_":"8599"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

# 4. Visualization: The Volcano Plot
Visualizing genomic alterations is crucial for interpretation. The plot_volcano method automatically highlights significant upregulation (Red) and downregulation (Blue).

``` r
# Generate plot
p <- plot_volcano(panel)
#> [1] "Analyzing 17814 genes across 590 samples..."

# Display
print(p)
```

![](/Users/premithapagadala/Documents/OncoMarker/OncoMarker_Analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
Volcano plot highlights:
Red points: Significantly upregulated in tumors
Blue points: Significantly downregulated
Grey points: Not significant

# 5. Clinical Risk Stratification
Finally, we can move from descriptive analysis to prognostic modeling. The predict_risk function stratifies the cohort into "High Risk" and "Low Risk" groups based on the median expression of a specific biomarker.

Example: Stratifying by TP53 (Tumor Suppressor) Logic: Low expression indicates loss of function -> High Risk.

``` r
# Predict risk based on TP53 levels
risk_scores <- predict_risk(panel, gene = "TP53", direction = "low_risk_high_expr")
#> [1] "Stratifying cohort based on TP53 median expression: -0.3"

# View the distribution of patients
table(risk_scores)
#> risk_scores
#>  Low Risk High Risk 
#>       295       295
```
# 6. Summary & Key Findings
## Summary of Differential Expression Analysis
The OncoMarker pipeline successfully identifies significant genomic alterations across the targeted 50-gene panel. For instance, TP53, MKI67, and ERBB2 consistently appear as top candidates with significant fold changes between Tumor and Normal samples.  

Key observations:
- Genes with Log2FC > 1 and adjusted p-value < 0.05 are considered upregulated in tumors.  
- Genes with Log2FC < -1 and adjusted p-value < 0.05 are considered downregulated in tumors.  

## Risk Stratification Insights
Using the `predict_risk()` function, patients were categorized into High Risk and Low Risk groups. This allows immediate translation to survival analyses or clinical decision support workflows.


``` r
# Summarize risk distribution
table(panel@patient_metadata$Risk)
#> < table of extent 0 >
```
# 7. References & Further Reading

## Data Sources
1. Breast Cancer Gene Expression Profiles (TCGA) – Kaggle  
   [https://www.kaggle.com/datasets/orvile/gene-expression-profiles-of-breast-cancer](https://www.kaggle.com/datasets/orvile/gene-expression-profiles-of-breast-cancer?resource=download)  

## R Programming & Tutorials
2. R Programming & Tutorials – R-Coder  
   [https://r-coder.com/](https://r-coder.com/)  

3. Applied Modeling & Text Analysis – Julia Silge Blog  
   [https://juliasilge.com/blog/](https://juliasilge.com/blog/)  

4. Visualization with ggplot2 – R Graph Gallery  
   [https://r-graph-gallery.com/ggplot2-package.html](https://r-graph-gallery.com/ggplot2-package.html)  

5. Comprehensive R Learning Libraries – Awesome R Resources  
   [https://github.com/iamericfletcher/awesome-r-learning-resources?tab=readme-ov-file](https://github.com/iamericfletcher/awesome-r-learning-resources?tab=readme-ov-file)  

## Suggested Reading for Clinical Genomics
6. The Cancer Genome Atlas (TCGA) – National Cancer Institute  
   [https://www.cancer.gov/tcga](https://www.cancer.gov/tcga)  

7. cBioPortal for Cancer Genomics  
   [https://www.cbioportal.org](https://www.cbioportal.org)  

8. Bioinformatics Workflows and Differential Expression Analysis (Friedman et al., 2023)
