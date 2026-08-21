# 01_preprocess_arrow_partitioned.R
# Purpose: Aggregate raw spectral data, remove AC artifact, normalize per epoch,
#          and write to a partitioned Parquet dataset for low-RAM out-of-core processing.
setwd("~/Downloads/zhao data raw")
library(tidyverse)
library(readxl)
library(arrow) # Required for writing Parquet datasets

# --- USER DEFINED VARIABLES ---
input_dir <- "data"
output_dir <- "parquet_dataset" # Directory where the Parquet structure will live
binsize <- 60 # Epoch bin size in minutes for stratified analysis
recordingslength <- 22 # Total recording duration in hours
upperlimit <- 100 # Maximum frequency (Hz) to retain
ac_power <- 50 # Mains power frequency to remove (50 Hz for Australia/Europe, 60 Hz for Americas)

# Calculate derived variables
shrinkage <- 4 / 60 / binsize
dbmax <- recordingslength / 4 * 60 * 60

dir.create(output_dir, showWarnings = FALSE)

# Processing loop
folders <- list.dirs(input_dir)[-1]

for (folder in folders) {
  files <- list.files(folder, full.names = TRUE, pattern = "\\.xlsx?$")

  for (file in files) {

    # 1. Parse the underscore-separated filename
    # Expected format: "PROJECT_M001_Vehicle_Cohort1.xlsx"
    base_name <- tools::file_path_sans_ext(basename(file))
    filename_parts <- strsplit(base_name, split = "_")[[1]]

    if(length(filename_parts) < 4) {
      warning(paste("Filename does not match expected format:", base_name))
      next
    }

    project_id     <- filename_parts[1]
    subject_id     <- filename_parts[2]
    treatment_cond <- filename_parts[3]
    cohort_id      <- filename_parts[4]

    # 2. Read, Filter, and Normalize
    df <- read_excel(file) %>%
      pivot_longer(!1:7, names_to = "frequency", values_to = "count") %>%
      mutate(frequency = as.numeric(frequency)) %>%
      filter(frequency <= upperlimit) %>%

      # Remove AC power artifact prior to normalization
      filter(frequency < (ac_power - 1) | frequency > (ac_power + 1)) %>%

      # Normalize power WITHIN each 4-second epoch
      group_by(`Epoch Index`) %>%
      mutate(totalpower = sum(count)) %>%
      ungroup() %>%

      # Calculate variables and append metadata
      mutate(
        relpower = count / totalpower,
        id = subject_id,
        treatment = treatment_cond,
        cohort = cohort_id,
        bin = ceiling(`Epoch Index` * shrinkage)
      ) %>%

      filter(`Epoch Index` <= dbmax) %>%
      select(Epoch.Index = `Epoch Index`, bin, frequency, relpower, id, treatment, cohort, `Main Score`)

    # 3. Write directly to disk using Nested Partitioning
    # This creates a folder structure like: cohort=Cohort1/treatment=Vehicle/frequency=0.195/
    write_dataset(
      df,
      path = output_dir,
      format = "parquet",
      partitioning = c("cohort", "treatment", "frequency"),
      basename_template = paste0(subject_id, "-{i}.parquet") # <--- THE FIX
    )
    # 4. Strictly clear memory for laptop performance
    rm(df)
    gc()
  }
}


