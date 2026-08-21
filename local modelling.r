# 02_modeling_arrow_unified.R
# Purpose: Unified out-of-core modeling script. Analyzes specific treatment vs. control pairs
#          by lazily loading targeted partitions from the Parquet dataset into parallel cores.

library(broom)
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(parallel)
library(arrow)

# --- 1. USER DEFINED VARIABLES & EXPERIMENTAL DESIGN MAP ---

dataset_dir <- "parquet_dataset"
output_dir <- "statistical/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# TOGGLE ANALYSIS TYPE HERE:
# "stratified" = State-stratified hour-by-hour (21 hours)
# "1hr"        = Continuous non-stratified epoch-by-epoch (first hour)
# "22hr"       = Continuous non-stratified binned (full recording)
analysis_type <- "stratified"

epoch_duration_sec <- 4
target_hour <- 1       # Used if analysis_type == "1hr"
bin_size_minutes <- 15 # Used if analysis_type == "22hr"

# Lookup Table: Defines which drug maps to which control and restricts analysis to that specific cohort.
comparisons_map <- tribble(
  ~comparison_name,      ~target_drug,       ~target_control, ~target_cohort,
  "MK1064_v_Veh",         "MK1064",            "MC",       "1",
#  "Seltorexant40_v_Veh",         "Seltorexant40mg",            "MC",       "2",
#  "Seltorexant60_v_Veh",         "Seltorexant60mg",            "MC",       "3",
#  "Tiaabine_v_Veh",         "Tiagabine",            "MC",       "4"
)

# --- 2. EXTRACT FREQUENCIES ---
# Query the dataset to get all unique frequencies present on the hard drive
ds <- open_dataset(dataset_dir)
unique_freq_bands <- ds %>% select(frequency) %>% distinct() %>% collect() %>% pull(frequency)
unique_freq_bands <- sort(unique_freq_bands)

# --- 3. SET UP CLUSTER ---
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  library(lme4)
  library(emmeans)
  library(dplyr)
  library(arrow)
})

# --- 4. EXECUTE MODELING LOOP ---
for (i in 1:nrow(comparisons_map)) {

  current_run <- comparisons_map[i, ]
  comp_name   <- current_run$comparison_name
  drug_name   <- current_run$target_drug
  ctrl_name   <- current_run$target_control
  cohort_name <- current_run$target_cohort

  print(paste("Processing:", comp_name))



  # Define treatments as a single vector
  target_treatments <- c(drug_name, ctrl_name)

  clusterExport(cl, varlist = c(
    "dataset_dir", "analysis_type", "epoch_duration_sec",
    "target_hour", "bin_size_minutes", "drug_name", "ctrl_name",
    "cohort_name", "target_treatments"
  ))

  results_list <- parLapply(cl, unique_freq_bands, function(band) {

    # 4a. Out-of-core bypass: Construct physical file paths to skip Arrow's filter engine
    drug_path <- file.path(dataset_dir, paste0("cohort=", cohort_name),
                           paste0("treatment=", drug_name), paste0("frequency=", band))

    ctrl_path <- file.path(dataset_dir, paste0("cohort=", cohort_name),
                           paste0("treatment=", ctrl_name), paste0("frequency=", band))

    band_data <- data.frame()

    # Lazily load and append Drug Data
    if (dir.exists(drug_path)) {
      d_data <- open_dataset(drug_path) %>% collect() %>%
        mutate(treatment = drug_name, frequency = band)
      band_data <- bind_rows(band_data, d_data)
    }

    # Lazily load and append Control Data
    if (dir.exists(ctrl_path)) {
      c_data <- open_dataset(ctrl_path) %>% collect() %>%
        mutate(treatment = ctrl_name, frequency = band)
      band_data <- bind_rows(band_data, c_data)
    }

    if(nrow(band_data) == 0) return(NULL)

    # Guarantee Control is the reference level for contrasts
    band_data$treatment <- factor(band_data$treatment, levels = c(ctrl_name, drug_name))

    # 4b. Format data based on chosen methodology
    bouts_df <- NULL

    if (analysis_type == "stratified") {
      band_data <- band_data %>%
        filter(bin <= 21, !is.na(`Main Score`), `Main Score` %in% c("Wake", "NREM", "REM")) %>%
        select(-Epoch.Index) %>%   # <--- THE FIX: Drop the raw Epoch.Index before renaming bin
        rename(score = `Main Score`, Epoch.Index = bin) %>%
        group_by(id, Epoch.Index, treatment, score) %>%
        summarise(relpower = mean(relpower), count = n(), .groups = 'drop') %>%
        mutate(Epoch.Index = as.factor(Epoch.Index))

      # Calculate vigilance state bouts if processing the 0Hz band
      if(band == 0) {
        bouts_df <- band_data %>%
          group_by(Epoch.Index, score, treatment) %>%
          summarise(Epochs = sum(count), .groups = 'drop')
      }

    } else if (analysis_type == "1hr") {
      epochs_per_hour <- 3600 / epoch_duration_sec
      min_epoch <- (target_hour - 1) * epochs_per_hour + 1
      max_epoch <- target_hour * epochs_per_hour

      band_data <- band_data %>%
        filter(Epoch.Index >= min_epoch, Epoch.Index <= max_epoch) %>%
        group_by(id, Epoch.Index, treatment) %>%
        summarise(relpower = mean(relpower), .groups = 'drop') %>%
        mutate(Epoch.Index = as.factor(Epoch.Index))

    } else if (analysis_type == "22hr") {
      epochs_per_bin <- (bin_size_minutes * 60) / epoch_duration_sec

      band_data <- band_data %>%
        mutate(Time.Bin = ceiling(Epoch.Index / epochs_per_bin)) %>%
        group_by(id, Time.Bin, treatment) %>%
        summarise(relpower = mean(relpower), .groups = 'drop') %>%
        mutate(Time.Bin = as.factor(Time.Bin))
    }

    # 4c. Fit Models
    tryCatch({
      if (analysis_type == "stratified") {
        model <- lmer(relpower ~ treatment * Epoch.Index * score + (1|id), data = band_data)
        emm <- emmeans(model, ~ treatment | Epoch.Index | score)
      } else if (analysis_type == "1hr") {
        model <- lmer(relpower ~ treatment * Epoch.Index + (1|id), data = band_data)
        emm <- emmeans(model, ~ treatment | Epoch.Index)
      } else {
        model <- lmer(relpower ~ treatment * Time.Bin + (1|id), data = band_data)
        emm <- emmeans(model, ~ treatment | Time.Bin)
      }

      # 4d. Package results and clear worker memory
      res <- list(
        emmeans_df = as.data.frame(emm),
        pairs_df   = broom::tidy(pairs(emm)),
        band       = band,
        bouts_df   = bouts_df
      )

      rm(model, emm, band_data)
      gc()

      return(res)

    }, error = function(e) { return(NULL) })
  })

  # --- 5. EXTRACT & SAVE ---
  pvals <- data.frame()
  fold_change <- data.frame()
  master_bouts <- data.frame()

  for(res in results_list){
    if(is.null(res)) next

    band <- res$band

    if (analysis_type == "stratified") {
      fold <- res$emmeans_df %>% select(treatment, emmean, Epoch.Index, score) %>%
        pivot_wider(names_from = Epoch.Index, values_from = emmean) %>% mutate(band = band)
      pairs <- res$pairs_df %>% select(Epoch.Index, score, contrast, p.value) %>% mutate(band = band)
      if(!is.null(res$bouts_df)) master_bouts <- res$bouts_df

    } else if (analysis_type == "1hr") {
      fold <- res$emmeans_df %>% select(treatment, emmean, Epoch.Index) %>%
        pivot_wider(names_from = Epoch.Index, values_from = emmean) %>% mutate(band = band)
      pairs <- res$pairs_df %>% select(Epoch.Index, contrast, p.value) %>% mutate(band = band)

    } else {
      fold <- res$emmeans_df %>% select(treatment, emmean, Time.Bin) %>%
        pivot_wider(names_from = Time.Bin, values_from = emmean) %>% mutate(band = band)
      pairs <- res$pairs_df %>% select(Time.Bin, contrast, p.value) %>% mutate(band = band)
    }

    fold_change <- bind_rows(fold_change, fold)
    pvals <- bind_rows(pvals, pairs)
  }

  # Write outputs
  if (analysis_type == "stratified") {
    write_csv(fold_change, file.path(output_dir, paste0(comp_name, "_stratified_fold.csv")))
    write_csv(pvals, file.path(output_dir, paste0(comp_name, "_stratified_pvals.csv")))
    if(nrow(master_bouts) > 0) {
      write_csv(master_bouts, file.path(output_dir, paste0(comp_name, "_stratified_bouts.csv")))
    }
  } else {
    write_csv(fold_change, file.path(output_dir, paste0(comp_name, "_continuous_", analysis_type, "_fold.csv")))
    write_csv(pvals, file.path(output_dir, paste0(comp_name, "_continuous_", analysis_type, "_pvals.csv")))
  }
}

stopCluster(cl)
print("Processing complete.")
