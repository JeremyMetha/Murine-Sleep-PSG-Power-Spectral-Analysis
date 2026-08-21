# 03_unified_heatmaps.R
# Purpose: Generate spectral ratio heatmaps and bout plots from modeled CSV outputs.
# Supports state-stratified, 1-hour epoch-by-epoch, and 22-hour binned methodologies.

library(tidyverse)
library(scico)
library(metR)

# --- 1. USER DEFINED VARIABLES & EXPERIMENTAL DESIGN MAP ---

input_dir <- "statistical/outputs"
output_dir <- "statistical/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# TOGGLE ANALYSIS TYPE HERE:
# "stratified" = State-stratified hour-by-hour (21 hours)
# "1hr"        = Continuous non-stratified epoch-by-epoch (first hour)
# "22hr"       = Continuous non-stratified binned (full recording)
analysis_type <- "stratified"

epoch_duration_sec <- 4
bin_size_minutes <- 15
states <- c("Wake", "NREM", "REM")

# Lookup Table: Must match the modeling script exactly to dynamically pull the right names
comparisons_map <- tribble(
  ~comparison_name,      ~target_drug,       ~target_control, ~target_cohort,
  "MK1064_v_Veh",         "MK1064",            "MC",       "1",
  #  "Seltorexant40_v_Veh",         "Seltorexant40mg",            "MC",       "2",
  #  "Seltorexant60_v_Veh",         "Seltorexant60mg",            "MC",       "3",
  #  "Tiaabine_v_Veh",         "Tiagabine",            "MC",       "4"
)

setwd(input_dir)

# --- 2. EXECUTE VISUALIZATION LOOP ---
for (i in 1:nrow(comparisons_map)) {

  current_run <- comparisons_map[i, ]
  comp_name   <- current_run$comparison_name
  drug_name   <- current_run$target_drug
  ctrl_name   <- current_run$target_control

  print(paste("Plotting:", comp_name))

  # Load Data
  if(analysis_type == "stratified") {
    fold_file <- paste0(comp_name, "_stratified_fold.csv")
    pval_file <- paste0(comp_name, "_stratified_pvals.csv")
    bout_file <- paste0(comp_name, "_stratified_bouts.csv")
  } else {
    fold_file <- paste0(comp_name, "_continuous_", analysis_type, "_fold.csv")
    pval_file <- paste0(comp_name, "_continuous_", analysis_type, "_pvals.csv")
    bout_file <- NULL
  }

  # Skip if files haven't been generated yet
  if(!file.exists(fold_file) | !file.exists(pval_file)) {
    print(paste("Files missing for", comp_name, "- skipping."))
    next
  }

  fold  <- read_csv(fold_file, show_col_types = FALSE)
  pvals <- read_csv(pval_file, show_col_types = FALSE)
  bouts <- if(!is.null(bout_file) && file.exists(bout_file)) read_csv(bout_file, show_col_types = FALSE) else NULL

  # Identify time columns dynamically
  exclude_cols <- if(analysis_type == "stratified") c("treatment", "score", "band") else c("treatment", "band")
  time_cols <- setdiff(colnames(fold), exclude_cols)

  # Format Spectral Data
  df <- fold %>%
    pivot_longer(cols = all_of(time_cols), names_to = "TimePoint", values_to = "value") %>%
    mutate(TimePoint = as.numeric(TimePoint))

  treatmentpower <- df %>% filter(treatment == drug_name)
  vehiclepower <- df %>% filter(treatment == ctrl_name)

  # Calculate Percentage Ratio
  if(analysis_type == "stratified"){
    ratiopower <- treatmentpower %>%
      left_join(vehiclepower %>% select(TimePoint, band, score, veh_value = value),
                by = c("TimePoint", "band", "score")) %>%
      mutate(ratio = (value / veh_value) * 100 - 100)
  } else {
    ratiopower <- treatmentpower %>%
      left_join(vehiclepower %>% select(TimePoint, band, veh_value = value),
                by = c("TimePoint", "band")) %>%
      mutate(ratio = (value / veh_value) * 100 - 100)
  }

  # Scale X-axis for 1hr epoch data to minutes
  if(analysis_type == "1hr"){
    ratiopower <- ratiopower %>% mutate(TimePoint = (TimePoint * epoch_duration_sec) / 60)
  }

  # Format P-values
  time_col_name <- if(analysis_type == "22hr") "Time.Bin" else "Epoch.Index"
  plot_data_sig <- pvals %>%
    filter(as.numeric(band) > 0) %>%
    mutate(
      TimePoint = as.numeric(!!sym(time_col_name)),
      band = as.numeric(band),
      log_p = -log10(p.value)
    )

  if(analysis_type == "1hr"){
    plot_data_sig <- plot_data_sig %>% mutate(TimePoint = (TimePoint * epoch_duration_sec) / 60)
  }

  # --- 3. PLOTTING FUNCTION ---
  plot_heatmap <- function(rp_data, pval_data, title, x_lab) {
    ggplot(rp_data, aes(x = TimePoint, y = as.numeric(band))) +
      geom_contour_fill(bins = 2000, aes(z = ratio)) +
      scale_fill_scico(palette = "vik", limits = c(-50, 100), values = c(0, 0.333, 1),
                       oob = scales::squish, name = "Spectral Ratio\n(% change)") +
      geom_contour(data = pval_data, bins = 2000,
                   aes(x = TimePoint, y = band, z = p.value, color = as.character(after_stat(level))),
                   breaks = c(0.001, 0.01, 0.05)) +
      scale_color_manual(values = c("0.05" = "black", "0.01" = "darkgreen", "0.001" = "yellow"),
                         name = "Significance (p)") +
      # Capped at 50Hz with pseudo-log transformation
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                         breaks = c(0.5, 1, 2, 5, 10, 20, 30, 50), limits = c(0, 50)) +
      scale_x_continuous(expand = c(0, 0)) +
      coord_cartesian(expand = FALSE) +
      theme_minimal() +
      labs(x = x_lab, y = "Frequency (Hz)", title = title)
  }

  # --- 4. RENDER & SAVE ---
  if(analysis_type == "stratified"){
    for(state in states){

      rp_sub <- ratiopower %>% filter(score == state)
      pval_sub <- plot_data_sig %>% filter(score == state)

      heatmap_title <- paste("Spectral Differences:", drug_name, "vs", ctrl_name, "-", state)
      hm_plot <- plot_heatmap(rp_sub, pval_sub, heatmap_title, "Time (h)")

      # Plot Bouts
      if(!is.null(bouts)) {
        bouts2 <- bouts %>% filter(score == state) %>% mutate(Epoch = as.numeric(Epoch.Index))
        bout_title <- paste("# Bouts of", state, "-", drug_name, "vs", ctrl_name)

        boutplot <- ggplot(bouts2, aes(fill = treatment, x = Epoch, y = Epochs)) +
          geom_area(alpha = 0.7, position = "identity") +
          theme_classic() +
          scale_fill_manual(labels = c(drug_name, ctrl_name), values = c("red", "blue")) +
          scale_y_continuous(expand = expansion(0)) +
          scale_x_continuous(expand = expansion(0), breaks = 0:21) +
          labs(x = "Time (h)", y = "Bout Count", fill = "Treatment", title = bout_title)

        ggsave(plot = boutplot, filename = file.path(output_dir, paste0(comp_name, "_", state, "_bouts.tiff")),
               width = 2400, height = 1800, units = "px", dpi = 300)
      }

      ggsave(plot = hm_plot, filename = file.path(output_dir, paste0(comp_name, "_", state, "_heatmap.tiff")),
             width = 2400, height = 1800, units = "px", dpi = 300)
    }

  } else {

    x_label <- if(analysis_type == "1hr") "Time (Minutes)" else paste0("Time Bins (", bin_size_minutes, " min)")
    heatmap_title <- paste("Continuous Spectral Differences:", drug_name, "vs", ctrl_name)

    hm_plot <- plot_heatmap(ratiopower, plot_data_sig, heatmap_title, x_label)

    ggsave(plot = hm_plot, filename = file.path(output_dir, paste0(comp_name, "_continuous_", analysis_type, "_heatmap.tiff")),
           width = 2400, height = 1800, units = "px", dpi = 300)
  }
}
