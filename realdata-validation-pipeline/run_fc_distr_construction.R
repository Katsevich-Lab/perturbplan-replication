# This is a script for fold change distribution construction for enhancer and gene analysis
library(ggplot2)
library(dplyr)

# load the estimated effect size
results_dir <- paste0("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/", experiment)

# Load positive control results including both gene_tss and known_enhancer
positive_cntrl_results <- readRDS(paste0(results_dir, "/discovery_es_at_scale.rds")) |> 
  dplyr::filter(target_type %in% c("gene_tss", "known_enhancer"))

# load the positive results for fold_change
positive_results <- readRDS(paste0(results_dir, "/discovery_es_at_scale_positive.rds"))

# Column information:
# positive_results: grna_target, response_id, fold_change_at_scale, se_fold_change_at_scale
# positive_cntrl_results: relative_expression_at_scale, expression_size_at_scale, etc.

# Step 1: Identify target-gene pairs with exactly 2 gRNAs in the main dataset
pairs_with_some_grnas <- positive_cntrl_results |>
  dplyr::group_by(grna_target, response_id) |>
  dplyr::summarise(num_grna = dplyr::n(), .groups = "drop") |>
  dplyr::filter(num_grna == guides_per_target) |>
  dplyr::select(grna_target, response_id)

# Step 2: Filter positive results for analysis based on window_def
cat("Using window definition:", window_def, "\n")

if (window_def == "enhancer_analysis") {
  # For enhancer analysis: focus on effect size 0.7 to 0.85
  positive_results_filtered <- positive_results |>
    dplyr::inner_join(pairs_with_some_grnas, by = c("grna_target", "response_id")) |>
    dplyr::filter(fold_change_at_scale >= 0.7 & fold_change_at_scale <= 0.85)
    
  # Define three ranges for new window specification
  range1 <- c(0.7, 0.75)   # 0.7 - 0.75
  range2 <- c(0.75, 0.8)   # 0.75 - 0.8  
  range3 <- c(0.8, 0.85)   # 0.8 - 0.85
  
} else if (window_def == "gene_analysis") {
  # For gene analysis: focus on effect size 0.6 to 0.75
  positive_results_filtered <- positive_results |>
    dplyr::inner_join(pairs_with_some_grnas, by = c("grna_target", "response_id")) |>
    dplyr::filter(fold_change_at_scale >= 0.6 & fold_change_at_scale <= 0.75)
    
  # Define three ranges for gene analysis
  range1 <- c(0.6, 0.65)   # 0.6 - 0.65
  range2 <- c(0.65, 0.7)   # 0.65 - 0.7
  range3 <- c(0.7, 0.75)   # 0.7 - 0.75
  
} else {
  stop("Invalid window_def. Must be 'enhancer_analysis' or 'gene_analysis'")
}

# Step 3: Count pairs in each range
range1_pairs <- positive_results_filtered |>
  dplyr::filter(fold_change_at_scale >= range1[1] & fold_change_at_scale <= range1[2]) |>
  dplyr::select(grna_target, response_id)

range2_pairs <- positive_results_filtered |>
  dplyr::filter(fold_change_at_scale >= range2[1] & fold_change_at_scale <= range2[2]) |>
  dplyr::select(grna_target, response_id)

range3_pairs <- positive_results_filtered |>
  dplyr::filter(fold_change_at_scale >= range3[1] & fold_change_at_scale <= range3[2]) |>
  dplyr::select(grna_target, response_id)

# Count available pairs in each range
n_range1 <- nrow(range1_pairs)
n_range2 <- nrow(range2_pairs)
n_range3 <- nrow(range3_pairs)

cat("Available pairs in each range:\n")
if (window_def == "enhancer_analysis") {
  cat("Range 1 (0.7-0.75):", n_range1, "\n")
  cat("Range 2 (0.75-0.8):", n_range2, "\n")
  cat("Range 3 (0.8-0.85):", n_range3, "\n")
} else if (window_def == "gene_analysis") {
  cat("Range 1 (0.6-0.65):", n_range1, "\n")
  cat("Range 2 (0.65-0.7):", n_range2, "\n")
  cat("Range 3 (0.7-0.75):", n_range3, "\n")
}

# Step 4: Include ALL pairs cumulatively for each tail type
# Light-tail: ALL pairs from range3 only (highest fold change range)
# Medium-tail: ALL pairs from range2 + range3 (medium to high fold changes)  
# Heavy-tail: ALL pairs from range1 + range2 + range3 (all fold change ranges)

# Create samples that include cumulative ranges
light_tail_samples <- list(
  sampled_1 = data.frame(),  # No range1 pairs for light tail
  sampled_2 = data.frame(),  # No range2 pairs for light tail
  sampled_3 = range3_pairs,  # ALL range3 pairs for light tail
  counts = c(0, 0, n_range3)
)

medium_tail_samples <- list(
  sampled_1 = data.frame(),  # No range1 pairs for medium tail
  sampled_2 = range2_pairs,  # ALL range2 pairs for medium tail
  sampled_3 = range3_pairs,  # ALL range3 pairs for medium tail
  counts = c(0, n_range2, n_range3)
)

heavy_tail_samples <- list(
  sampled_1 = range1_pairs,  # ALL range1 pairs for heavy tail
  sampled_2 = range2_pairs,  # ALL range2 pairs for heavy tail
  sampled_3 = range3_pairs,  # ALL range3 pairs for heavy tail
  counts = c(n_range1, n_range2, n_range3)
)

cat("Light-tail (range3 only):", light_tail_samples$counts, "Total:", sum(light_tail_samples$counts), "\n")
cat("Medium-tail (range2+3):", medium_tail_samples$counts, "Total:", sum(medium_tail_samples$counts), "\n")
cat("Heavy-tail (range1+2+3):", heavy_tail_samples$counts, "Total:", sum(heavy_tail_samples$counts), "\n")

# Step 5: Create dataframes for each tail type
create_tail_dataframe <- function(samples, tail_name) {
  # Combine all sampled pairs
  combined_pairs <- dplyr::bind_rows(samples$sampled_1, samples$sampled_2, samples$sampled_3)
  
  # Join with positive results and target-level info
  target_level_info <- positive_cntrl_results |>
    dplyr::group_by(grna_target, response_id) |>
    dplyr::summarise(
      relative_expression = first(relative_expression_at_scale),
      expression_size = first(expression_size_at_scale),
      num_trt_cells_full = sum(num_oracle_cells),
      min_cells_per_grna = min(num_oracle_cells),
      grna_per_target = dplyr::n(),
      num_total_plan_cells = first(num_total_plan_cells),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      num_trt_cells_conservative = min_cells_per_grna * grna_per_target
    ) |>
    dplyr::select(-min_cells_per_grna)
  
  # compute the gRNA_sd when experiment is Ray
  gRNA_sd <- switch (experiment,
    Ray = {
      positive_cntrl_results |> 
        dplyr::group_by(grna_target, response_id) |> 
        dplyr::summarise(var_per_target = var(fold_change_at_scale, na.rm = TRUE), .groups = "drop") |> 
        dplyr::ungroup() |> 
        dplyr::summarise(grna_sd = sqrt(mean(var_per_target, na.rm = TRUE))) |> 
        dplyr::pull()
    },
    Gasperini = {
      0
    }
  )

  # Create final dataframe
  result_df <- combined_pairs |>
    dplyr::left_join(positive_results |> dplyr::select(grna_target, response_id, fold_change_at_scale, se_fold_change_at_scale), 
                     by = c("grna_target", "response_id")) |>
    dplyr::left_join(target_level_info, by = c("grna_target", "response_id")) |>
    dplyr::rename(fold_change = fold_change_at_scale, se_fold_change = se_fold_change_at_scale) |>
    dplyr::mutate(gRNA_sd = gRNA_sd) |>  # Set gRNA_sd to 0 as per previous instructions
    dplyr::filter(!is.na(fold_change) & !is.na(relative_expression) & !is.na(expression_size) & 
                  !is.na(num_trt_cells_full) & !is.na(num_trt_cells_conservative) & !is.na(num_total_plan_cells))
  
  return(result_df)
}

# Create dataframes for each tail type
fc_distr_light_tail_df <- create_tail_dataframe(light_tail_samples, "light_tail")
fc_distr_medium_tail_df <- create_tail_dataframe(medium_tail_samples, "medium_tail")
fc_distr_heavy_tail_df <- create_tail_dataframe(heavy_tail_samples, "heavy_tail")

# Step 6: Load negative control data (non-targeting only)
negative_cntrl_results <- readRDS(paste0(results_dir, "/discovery_es_at_scale.rds")) |> 
  dplyr::filter(target_type %in% c("non-targeting"))

# Get target-level summary info for negative controls
negative_target_level_info <- negative_cntrl_results |>
  dplyr::group_by(grna_target, response_id) |>
  dplyr::summarise(
    relative_expression = first(relative_expression_at_scale),
    expression_size = first(expression_size_at_scale),
    num_trt_cells_full = sum(num_oracle_cells),
    num_total_plan_cells = first(num_total_plan_cells),
    .groups = "drop"
  ) |>
  dplyr::filter(!is.na(relative_expression) & !is.na(expression_size) & 
                !is.na(num_trt_cells_full) & !is.na(num_total_plan_cells))

# Function to create discovery_pairs
create_discovery_pairs <- function(positive_df, tail_name) {
  # Calculate required sample sizes
  n_positive <- nrow(positive_df)
  n_total <- ceiling(n_positive / positive_proportion)
  n_negative <- n_total - n_positive
  
  cat("Creating", tail_name, "pairs: ", n_positive, "positive +", n_negative, "negative =", n_total, "total pairs\n")
  
  # Sample negative control pairs
  set.seed(123)  # For reproducibility
  sampled_negative_pairs <- negative_target_level_info |>
    dplyr::slice_sample(n = min(n_negative, nrow(negative_target_level_info))) |>
    dplyr::mutate(
      pair_id = paste(grna_target, response_id, sep = "_"),
      pair_type = "negative"
    ) |>
    dplyr::select(pair_id, pair_type)
  
  # Create positive pairs
  positive_pairs <- positive_df |>
    dplyr::mutate(
      pair_id = paste(grna_target, response_id, sep = "_"),
      pair_type = "positive"
    ) |>
    dplyr::select(pair_id, pair_type)
  
  # Combine positive and negative pairs
  discovery_pairs <- dplyr::bind_rows(positive_pairs, sampled_negative_pairs)
  
  return(discovery_pairs)
}

# Create discovery_pairs for each tail type
discovery_pairs_light_tail <- create_discovery_pairs(fc_distr_light_tail_df, "light_tail")
discovery_pairs_medium_tail <- create_discovery_pairs(fc_distr_medium_tail_df, "medium_tail")
discovery_pairs_heavy_tail <- create_discovery_pairs(fc_distr_heavy_tail_df, "heavy_tail")

# Step 7: Save results
# Create output directory if it doesn't exist
output_dir <- paste0("reproduction-code-prep/realdata-validation-pipeline/intermediate-files/", experiment)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the three dataframes
saveRDS(fc_distr_light_tail_df, paste0(output_dir, "/fc_distr_light_tail_", window_def, ".rds"))
saveRDS(fc_distr_medium_tail_df, paste0(output_dir, "/fc_distr_medium_tail_", window_def, ".rds"))
saveRDS(fc_distr_heavy_tail_df, paste0(output_dir, "/fc_distr_heavy_tail_", window_def, ".rds"))

# Save the discovery_pairs
saveRDS(discovery_pairs_light_tail, paste0(output_dir, "/discovery_pairs_light_tail_", window_def, ".rds"))
saveRDS(discovery_pairs_medium_tail, paste0(output_dir, "/discovery_pairs_medium_tail_", window_def, ".rds"))
saveRDS(discovery_pairs_heavy_tail, paste0(output_dir, "/discovery_pairs_heavy_tail_", window_def, ".rds"))

cat("Fold change distribution construction completed for", window_def, "\n")