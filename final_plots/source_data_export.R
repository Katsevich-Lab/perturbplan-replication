# Export source data for every figure produced by an R script in this folder
# into a single .xlsx file with one sheet per panel (sheets prefixed by figure).
#
# Skipped (Keynote-built, no R source data): Figure 1, Figure 2, Extended Table 1.
library(writexl)
library(dplyr)

results_dir   <- "reproduction-code-prep"
plots_dir     <- "reproduction-code-prep/final_plots/figures"
out_xlsx      <- file.path(plots_dir, "source_data.xlsx")

# Accumulator for all sheets
all_sheets <- list()
add_sheets <- function(prefix, sheets) {
  sheets <- sheets[!vapply(sheets, function(df) is.null(df) || nrow(df) == 0, logical(1))]
  names(sheets) <- paste0(prefix, "_", names(sheets))
  all_sheets <<- c(all_sheets, sheets)
}

# Helper: pull $data from a (possibly patchworked) ggplot object
plot_data <- function(p) p$data

# Helper: pull $data from a (possibly patchworked) ggplot object
plot_data <- function(p) p$data

# ============================================================================
# Figure 3
# ============================================================================
fig3_a <- readRDS(sprintf("%s/simulation-benchmarking/benchmarking_plot.rds", results_dir))
fig3_b <- readRDS(sprintf("%s/saturation-curve-fitting/results/saturation_curve_plot.rds", results_dir))
fig3_c <- readRDS(sprintf("%s/realdata-validation-pipeline/PerturbPlan/results/Gasperini_gene_analysis_plot.rds", results_dir))
fig3_d <- readRDS(sprintf("%s/realdata-validation-pipeline/PerturbPlan/results/Ray_enhancer_analysis_plot.rds", results_dir))

add_sheets("fig3", list(
  a_sim_benchmark   = plot_data(fig3_a),
  b_read_umi_valid  = plot_data(fig3_b),
  c_Gasperini_valid = plot_data(fig3_c),
  d_Ray_valid       = plot_data(fig3_d)
))

# ============================================================================
# Figure 4 (Perturb-seq case study, three datasets)
# ============================================================================
f4_dir <- "reproduction-code-prep/case-study/perturb-seq-case-study/results/"
load_three <- function(rel_path) {
  dplyr::bind_rows(
    readRDS(paste0(f4_dir, "Morris/",    rel_path)) |> dplyr::mutate(assay = "Morris"),
    readRDS(paste0(f4_dir, "Gasperini/", rel_path)) |> dplyr::mutate(assay = "Gasperini"),
    readRDS(paste0(f4_dir, "Replogle/",  rel_path)) |> dplyr::mutate(assay = "Replogle")
  )
}

power_over_tpm  <- load_three("power_over_tpm_df.rds")
power_over_fc   <- load_three("power_over_fc_df.rds") |>
  dplyr::mutate(minimum_effect_size_percent = round((1 - minimum_fold_change) * 100, 1))

opt_power_tpm   <- load_three("power_opt_tpm_df.rds")
opt_cost_tpm_df <- opt_power_tpm |> dplyr::filter(total_cost == minimum_cost)

opt_power_fc    <- load_three("power_opt_fc_df.rds")
opt_cost_fc_df  <- opt_power_fc |> dplyr::filter(total_cost == minimum_cost)

add_sheets("fig4", list(
  b_power_vs_tpm        = power_over_tpm  |> dplyr::filter(TPM_threshold <= 1000),
  c_power_vs_eff_size   = power_over_fc   |> dplyr::filter(minimum_effect_size_percent <= 50),
  d_opt_cells_reads_tpm = opt_cost_tpm_df |> dplyr::select(assay, TPM_threshold, num_captured_cells, sequenced_reads_per_cell, minimum_cost),
  e_opt_cost_vs_tpm     = opt_cost_tpm_df |> dplyr::select(assay, TPM_threshold, minimum_cost),
  f_opt_cells_reads_fc  = opt_cost_fc_df  |> dplyr::select(assay, minimum_fold_change, num_captured_cells, sequenced_reads_per_cell, minimum_cost),
  g_opt_cost_vs_fc      = opt_cost_fc_df  |> dplyr::select(assay, minimum_fold_change, minimum_cost)
))

# ============================================================================
# Figure 5 (TAP-seq vs Perturb-seq case study)
# ============================================================================
f5_dir <- "reproduction-code-prep/case-study/TAP-seq-case-study/results/"

sat_plot <- readRDS(paste0(f5_dir, "saturation_plotting.rds"))
# Patchwork stored the saturation curves in the "inset" slot and the
# mapping-efficiency barplot at the top level (see patchwork inset_element layout).
panel_5a_curves <- sat_plot$patches$plots[[1]]$data |>
  dplyr::select(assay, reads_per_cell_sequenced, reads_per_cell_mapped, library_size, umis_per_cell_goi)
panel_5a_mapping_eff <- sat_plot$data

umi_plot <- readRDS(paste0(f5_dir, "UMI_comparison_plotting.rds"))
panel_5b <- plot_data(umi_plot) |>
  dplyr::select(response_id, TAP_seq, Perturb_seq) |>
  dplyr::rename(TAP_seq_UMI_per_cell = TAP_seq, Perturb_seq_UMI_per_cell = Perturb_seq)

power_tpm_plot <- readRDS(paste0(f5_dir, "power_versus_TPM_plotting.rds"))
panel_5c <- plot_data(power_tpm_plot) |>
  dplyr::select(assay, UMI_threshold, overall_power, dplyr::any_of("minimum_fold_change"))

power_fc_plot <- readRDS(paste0(f5_dir, "power_versus_FC_plotting.rds"))
panel_5d <- plot_data(power_fc_plot) |>
  dplyr::filter(minimum_fold_change >= 0.75) |>
  dplyr::mutate(minimum_effect_size_percent = round((1 - minimum_fold_change) * 100, 1)) |>
  dplyr::filter(minimum_effect_size_percent <= 50) |>
  dplyr::select(assay, minimum_effect_size_percent, overall_power)

eff_plot <- readRDS(paste0(f5_dir, "cell_read_cost_primer_eff_plotting.rds"))
# patchwork swapped main and inset: cost-savings is in patches$plots[[1]]
panel_5e <- eff_plot$patches$plots[[1]]$data |>
  dplyr::filter(arm_primer_efficiency_bound) |>
  dplyr::mutate(TAP_seq_cost_savings = 1 - 1 / Perturb_TAP_cost_ratio) |>
  dplyr::select(primer_efficiency_bound, Perturb_seq, TAP_seq, Perturb_TAP_cost_ratio, TAP_seq_cost_savings)

ratio_plot <- readRDS(paste0(f5_dir, "cell_read_cost_cost_ratio_plotting.rds"))
panel_5f <- plot_data(ratio_plot) |>
  dplyr::filter(arm_ratio) |>
  dplyr::mutate(seq_to_lib_cost_ratio = 1 / cost_ratio,
                TAP_seq_cost_savings  = 1 - 1 / Perturb_TAP_cost_ratio) |>
  dplyr::select(cost_ratio, seq_to_lib_cost_ratio, primer_efficiency_bound,
                Perturb_seq, TAP_seq, Perturb_TAP_cost_ratio, TAP_seq_cost_savings)

add_sheets("fig5", list(
  a_saturation_curves   = panel_5a_curves,
  a_mapping_eff_inset   = panel_5a_mapping_eff,
  b_umi_at_saturation   = panel_5b,
  c_power_vs_umi_thresh = panel_5c,
  d_power_vs_eff_size   = panel_5d,
  e_savings_primer_eff  = panel_5e,
  f_savings_cost_ratio  = panel_5f
))

# ============================================================================
# Extended Figure 1 (dataset metrics bar plots)
# ============================================================================
ext1_summary <- readRDS("reproduction-code-prep/final_plots/Extended_figure_1/cells_and_reads_summary.rds")
ext1_df <- do.call(rbind, lapply(ext1_summary, function(x) {
  data.frame(
    dataset                  = x$dataset_name %||% NA_character_,
    num_cells                = as.numeric(x$num_cells),
    umi_per_cell             = as.numeric(x$umi_per_cell),
    sequenced_reads_per_cell = as.numeric(x$sequenced_reads_per_cell),
    stringsAsFactors = FALSE
  )
}))
add_sheets("extfig1", list(dataset_metrics = ext1_df))

# ============================================================================
# Extended Figure 3 (saturation curve validation + downsampling)
# ============================================================================
sat_dir <- "reproduction-code-prep/saturation-curve-fitting/results"
add_sheets("extfig3", list(
  ab_validation_curves = readRDS(sprintf("%s/saturation_curve_comparison.rds", sat_dir)),
  c_K562_downsampled   = readRDS(sprintf("%s/saturation_curve_comparison_varying_sequencing_depth.rds", sat_dir)),
  d_TCD8_downsampled   = readRDS(sprintf("%s/saturation_curve_comparison_varying_sequencing_depth_tcd8.rds", sat_dir)),
  e_K562_true_curve    = readRDS(sprintf("%s/read_umi_df_k562.rds", sat_dir)),
  f_TCD8_true_curve    = readRDS(sprintf("%s/read_umi_df_tcd8.rds", sat_dir)),
  K562_downsamp_train  = readRDS(sprintf("%s/read_umi_df_varying_sd.rds", sat_dir)),
  TCD8_downsamp_train  = readRDS(sprintf("%s/read_umi_df_varying_sd_tcd8.rds", sat_dir))
))

# ============================================================================
# Extended Figure 4 (K562 Gasperini vs K562 10x case study)
# ============================================================================
load_two_k562 <- function(rel_path) {
  dplyr::bind_rows(
    readRDS(paste0(f4_dir, "Gasperini/", rel_path)) |> dplyr::mutate(assay = "K562 (Gasperini)"),
    readRDS(paste0(f4_dir, "K562_10x/",  rel_path)) |> dplyr::mutate(assay = "K562 (10x)")
  )
}

ext4_power_tpm <- load_two_k562("power_over_tpm_df.rds") |> dplyr::filter(TPM_threshold <= 1000)
ext4_power_fc  <- load_two_k562("power_over_fc_df.rds") |>
  dplyr::mutate(minimum_effect_size_percent = round((1 - minimum_fold_change) * 100, 1)) |>
  dplyr::filter(minimum_effect_size_percent <= 50)

ext4_opt_tpm <- load_two_k562("power_opt_tpm_df.rds") |> dplyr::filter(total_cost == minimum_cost)
ext4_opt_fc  <- load_two_k562("power_opt_fc_df.rds")  |> dplyr::filter(total_cost == minimum_cost)

# Panels a, b: per-gene TPM and dispersion comparison K562 (Gasperini) vs K562 (10x)
data("K562_10x", package = "perturbplan")
default_expr <- perturbplan:::get_pilot_data_from_package("K562")$baseline_expression_stats |>
  dplyr::filter(!is.na(relative_expression)) |>
  dplyr::rename(default_expr = relative_expression, default_size = expression_size)
custom_expr  <- K562_10x$baseline_expression_stats |>
  dplyr::filter(!is.na(relative_expression)) |>
  dplyr::rename(custom_expr = relative_expression, custom_size = expression_size)
ext4_ab <- dplyr::inner_join(default_expr, custom_expr, by = "response_id") |>
  dplyr::transmute(
    response_id,
    K562_Gasperini_TPM   = default_expr * 1e6,
    K562_10x_TPM         = custom_expr * 1e6,
    K562_Gasperini_size  = default_size,
    K562_10x_size        = custom_size
  )

add_sheets("extfig4", list(
  ab_tpm_and_dispersion = ext4_ab,
  c_power_vs_tpm        = ext4_power_tpm,
  d_power_vs_eff_size   = ext4_power_fc,
  e_opt_cells_reads_tpm = ext4_opt_tpm |> dplyr::select(assay, TPM_threshold, num_captured_cells, sequenced_reads_per_cell, minimum_cost),
  f_opt_cells_reads_fc  = ext4_opt_fc  |> dplyr::select(assay, minimum_fold_change, num_captured_cells, sequenced_reads_per_cell, minimum_cost),
  g_opt_cost_vs_tpm     = ext4_opt_tpm |> dplyr::select(assay, TPM_threshold, minimum_cost),
  h_opt_cost_vs_fc      = ext4_opt_fc  |> dplyr::select(assay, minimum_fold_change, minimum_cost)
))

# ============================================================================
# Cost minimization demonstration (recompute from scratch into a child env)
# ============================================================================
source("reproduction-code-prep/final_plots/cost_minimization_demonstration.R")
add_sheets("fig1", list(
  f_equi_power_curve = power_data,
  f_equi_cost_curve  = cost_data,
  f_optimal_design   = optimal_point
))

# ============================================================================
# Supplementary FDP curve figure (recompute from scratch)
# ============================================================================
source("reproduction-code-prep/final_plots/Supplementary-figure-FDP-curve.R")
fdp_df <- analytic_results |>
  dplyr::filter(num_total_cells == 37888, reads_per_cell %in% c(12287, 22117))
add_sheets("supp_fdp", list(FDP_vs_cutoff = fdp_df))

# ============================================================================
# Write the consolidated workbook
# ============================================================================
prefix_order <- c("fig1", "fig3", "fig4", "fig5", "extfig1", "extfig3", "extfig4", "supp_fdp")
sheet_prefix <- sub("_.*", "", names(all_sheets))
all_sheets <- all_sheets[order(match(sheet_prefix, prefix_order))]
all_sheets <- lapply(all_sheets, as.data.frame)
names(all_sheets) <- substr(names(all_sheets), 1, 31)
writexl::write_xlsx(all_sheets, path = out_xlsx)
message("wrote ", out_xlsx, " (", length(all_sheets), " sheets)")


`%||%` <- function(a, b) if (is.null(a)) b else a
