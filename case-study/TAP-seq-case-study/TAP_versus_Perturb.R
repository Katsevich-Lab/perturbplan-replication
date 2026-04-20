# This is a script comparing TAP-seq and Perturb-seq
library(dplyr)
library(perturbplan)
set.seed(1)

# source necessary files
source("reproduction-code-prep/case-study/helper.R")

# load sceptre object
sceptre_object <- readRDS(paste0(.get_config_path("LOCAL_RAY_2025_RAW_DATA_DIR"), "raw/sceptre_output/final_sceptre_object.rds"))

# load TAP-seq pilot data
pilot_data <- readRDS("reproduction-code-prep/case-study/auxiliary-files/discovery_es_at_scale.rds")

# load discovery table
discovery_table <- obtain_discovery_table(dataset_name = "Ray")

# path to power results
dir_to_results <- "reproduction-code-prep/case-study/TAP-seq-case-study/results/"
if(!dir.exists(dir_to_results)){
  dir.create(dir_to_results, recursive = TRUE)
}

############################## Obtain fixed parameters #########################
fixed_parameters <- readRDS("reproduction-code-prep/case-study/fixed_parameters.rds")
fdr_target <- fixed_parameters$fdr_target
precision <- fixed_parameters$precision
mapping_efficiency_perturb <- fixed_parameters$mapping_efficiency
mapping_efficiency_tap <- fixed_parameters$mapping_efficiency_tap
gRNA_variability <- fixed_parameters$gRNA_variability
prop_non_null <- fixed_parameters$prop_non_null
control_group <- fixed_parameters$control_group
side <- fixed_parameters$side
grid_size <- fixed_parameters$grid_size
cost_per_captured_cell <- fixed_parameters$cost_per_captured_cell
cost_per_million_reads <- fixed_parameters$cost_per_million_reads
power_target <- fixed_parameters$power_target
UMI_threshold_default <- 1
minimum_fold_change <- 0.9

# specify the parameter grid for retrospective analyses
UMI_threshold_retrospective_list <- 10^{seq(-1, 1, length.out = 100)}
minimum_fold_change_retrospective_list <- fixed_parameters$minimum_fold_change_retrospective_list

############################ import TAP-seq pilot data #########################
# Extract reads per cell
reads_per_cell <- round((readRDS("reproduction-code-prep/case-study/auxiliary-files/summary_persample.rds") |> 
                           as.data.frame() |>
                           dplyr::summarise(reads_per_cell = sum(num_cells * avg_reads) / sum(num_cells)) |> 
                           dplyr::pull()) / mapping_efficiency_perturb)

# specify gRNAs per target
gRNAs_per_target <- round(sceptre_object@grna_target_data_frame |> 
  dplyr::group_by(grna_target) |>
  dplyr::summarise(grna_per_target = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean(grna_per_target)) |>
  dplyr::pull())

# specify non-targeting gRNAs
non_targeting_gRNAs <- 51

# specify number of targets
num_guides <- length(unique(sceptre_object@grna_target_data_frame$grna_id))

# specify number of target
num_targets <- length(unique(sceptre_object@grna_target_data_frame$grna_target))

# MOI parameter
MOI <- mean(sceptre_object@grnas_per_cell)

# number of total cells
num_total_cells <- unique(pilot_data$num_total_plan_cells)

# extract cells per target
cells_per_target <- round(num_total_cells * MOI * gRNAs_per_target / num_guides)

##################### Obtain baseline expression statistics ####################
# specify the primer efficiency tolerance (lower bound on share-normalized
# TAP/Perturb ratio; 0 keeps all target genes)
primer_efficiency_threshold <- 0

# filter based on the UMI_threshold_default
baseline_info <- get_baseline_expression(primer_efficiency_threshold = primer_efficiency_threshold, 
                                         pilot_data = pilot_data)
baseline_expression_stats_perturb <- baseline_info$baseline_expression_stats_perturb
baseline_expression_stats_tap <- baseline_info$baseline_expression_stats_tap
library_parameters_perturb <- baseline_info$library_parameters_perturb
library_parameters_tap <- baseline_info$library_parameters_tap

#################### Perform the retrospective analyses ########################
# Perform power analysis
source("reproduction-code-prep/case-study/TAP-seq-case-study/retrospective_analysis.R")

# plot for retrospective analysis
source("reproduction-code-prep/case-study/TAP-seq-case-study/plot_retrospective_analysis.R")

#################### Perform the prospective analyses ##########################
# Perform cost-power analysis
source("reproduction-code-prep/case-study/TAP-seq-case-study/prospective_analysis.R")

# plot for retrospective analysis
source("reproduction-code-prep/case-study/TAP-seq-case-study/plot_prospective_analysis.R")
