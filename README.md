# Murine-Sleep-PSG-Power-Spectral-Analysis

This repository contains the R-based data processing and statistical modeling pipeline for analyzing laboratory rodent polysomnography (PSG) and quantitative EEG (qEEG) data. 

The methodologies executed in this code support the analyses described in our methods paper: *Laboratory Rodent Polysomnography and the Discovery of Hypnotic Orexin Receptor Antagonists: Methodological Considerations and Approaches* (Jacobson, Metha, et al.).

## Overview

High-resolution PSG recordings (e.g., 4-second epochs over 22+ hours) generate massive datasets. Traditional in-memory R workflows often bottleneck modern laptops[cite: 5]. This pipeline leverages **Apache Arrow** and partitioned Parquet datasets to perform "out-of-core" data processing. By keeping memory usage strictly constrained, this pipeline allows for complex, parallelized mixed-effects modeling across hundreds of frequency bands on standard desktop hardware.

The pipeline is split into three modular scripts:
1.  `01_preprocess_arrow_partitioned.R`: Cleans raw output, removes AC artifacts, normalizes power per epoch, and writes data to disk.
2.  `02_modeling_arrow_unified.R`: Lazily loads target data into a parallel cluster and fits `lme4` mixed-effects models.
3.  `03_unified_heatmaps.R`: Generates relative power spectral heatmaps with significance contours and vigilance bout plots.

## Data Naming Convention

To guarantee strict within-subject experimental control (ensuring a drug is only compared to the baseline vehicle from its exact experimental cohort), raw data files **must** follow this underscore-separated naming convention:

`[ProjectID]_[SubjectID]_[Treatment]_[Cohort].xlsx`

**Examples:**
*   `STUDY1_M001_Vehicle_Cohort1.xlsx`
*   `STUDY1_M001_Seltorexant40mg_Cohort1.xlsx`

## System Requirements

*   **R** (>= 4.1.0)
*   **Operating System:** Windows, macOS, or Linux. 
*   **Hardware:** Multi-core processor (parallelization uses `detectCores() - 1`). 16GB RAM recommended.

### R Dependencies
Ensure the following packages are installed:
```R
install.packages(c("tidyverse", "readxl", "arrow", "lme4", "emmeans", "parallel", "broom", "scico", "metR"))
