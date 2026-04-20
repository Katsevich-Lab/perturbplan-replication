# This is a Rscript preparing for data preprocessing for Gasperini data
# source the helper function
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")

# create preprocessing directory
preprocessing_dir <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Gasperini"
dir.create(preprocessing_dir, recursive = TRUE)

##################### 1. Preprocess the at-scale data ##########################
# preprocess the molecule_info.h5 for gene expression
path_to_plan_data <- obtain_file_path(experiment = "Gasperini")
SRR_dirs_list <- list.files(path = path_to_plan_data["raw"], pattern = "^SRR", full.names = TRUE)
summary_perSRR <- matrix(0, nrow = length(SRR_dirs_list), ncol = 2,
                         dimnames = list(SRR = basename(SRR_dirs_list),Summary = c("num_cells", "avg_reads")))

# obtain srr_labels
srr_labels <- readRDS(sprintf("%s/srr_labels.rds", path_to_plan_data["raw"]))
cell_id_list <- NULL

# loop over different SRRs
for (SRR_dir in SRR_dirs_list) {
  
  # make a director under this path
  path_to_save <- sprintf("%s/processed", SRR_dir)
  if(!dir.exists(path_to_save)){
    dir.create(path_to_save)
  }
  
  # if the file is already there then the QC should be performed or not
  if(file.exists(sprintf("%s/QC_data.rds", path_to_save))){
    QC_data <- readRDS(sprintf("%s/QC_data.rds", path_to_save))
  }else{
    # process the raw data
    QC_data <- obtain_qc_data(path_to_outs_folder = sprintf("%s/outs", SRR_dir),
                              library_index = srr_labels[basename(SRR_dir), "gene_expression_label"])
    
    # save the QC-ed data
    saveRDS(QC_data, file = sprintf("%s/QC_data.rds", path_to_save))
  }
                            
  # extract the summary statistics
  summary_perSRR[basename(SRR_dir), ] <- summary_data(QC_data)
  
  # extract cell id
  cell_id_list <- c(cell_id_list, unique(QC_data$cell_id))
  
  # delete the variable
  rm(QC_data)
}

# save the summary information
summary_SRR_dir <- sprintf("%s/summary_SRR.rds", preprocessing_dir)
if(!file.exists(summary_SRR_dir)){
  summary_SRR <- list(
    minimum_num_cells = min(summary_perSRR[, "num_cells"]),
    minimum_num_reads_per_cell = min(summary_perSRR[, "avg_reads"]),
    num_SRRs = length(SRR_dirs_list)
  )
  saveRDS(summary_SRR, summary_SRR_dir)
}

###################### 2. Construct the discovery pairs ########################
# get the grna group information
grna_at_scale_table <- readr::read_tsv(file = paste0(path_to_plan_data["processed"], "raw/GSE120861_all_deg_results.at_scale.txt"))
grna_id_to_group_at_scale <- readr::read_tsv(file = paste0(path_to_plan_data["processed"], "raw/GSE120861_grna_groups.at_scale.txt"),
                                             col_types = "cc", col_names = c("grna_group", "barcode"))

# preprocess at-scale discovery
discovery_pairs <- grna_at_scale_table |> 
  dplyr::select(gRNA_group, ENSG, site_type) |>
  dplyr::distinct() |>
  dplyr::rename(grna_group = gRNA_group) |>
  dplyr::rename(target_type = site_type) |>
  dplyr::filter(target_type != "TSS") |>
  dplyr::mutate(target_type = factor(target_type),
                target_type = forcats::fct_recode(target_type, 
                                                  gene_tss = "selfTSS",
                                                  candidate_enhancer = "DHS",
                                                  known_enhancer = "positive_ctrl",
                                                  "non-targeting" = "NTC")) |>
  dplyr::left_join(grna_id_to_group_at_scale, by = "grna_group", relationship = "many-to-many") |>
  dplyr::rename(grna_id = barcode) |>
  dplyr::distinct() |>
  dplyr::rename(response_id = ENSG, grna_target = grna_group)

# intersect with the processed data for at_scale data
grna_counts <- ondisc::read_odm(odm_fp = sprintf("%s/processed/grna_assignment/matrix.odm", path_to_plan_data["processed"]), 
                                metadata_fp = sprintf("%s/processed/grna_assignment/metadata.rds", path_to_plan_data["processed"]))
gene_counts <- ondisc::read_odm(odm_fp = sprintf("%s/processed/gene/matrix.odm", path_to_plan_data["processed"]), 
                                metadata_fp = sprintf("%s/processed/gene/metadata.rds", path_to_plan_data["processed"]))
discovery_pairs_at_scale <- discovery_pairs |>
  dplyr::filter(grna_id %in% rownames(grna_counts@feature_covariates)) |>
  dplyr::filter(response_id %in% rownames(gene_counts@feature_covariates)) |>
  dplyr::mutate(target_type = factor(target_type))

# compute cells per guide
grna_matrix <- grna_counts[[unique(discovery_pairs$grna_id), ]]
rownames(grna_matrix) <- unique(discovery_pairs$grna_id)
colnames(grna_matrix) <- stringi::stri_sub(rownames(grna_counts@cell_covariates), 1, -10)
cells_per_guide <- apply(grna_matrix[, intersect(cell_id_list, colnames(grna_matrix))], 1, sum)
cells_per_guide_df <- data.frame(
  grna_id = names(cells_per_guide),
  num_oracle_cells = cells_per_guide,
  num_total_plan_cells = length(intersect(cell_id_list, colnames(grna_matrix)))
)
discovery_pairs_extended <- discovery_pairs_at_scale |> dplyr::left_join(cells_per_guide_df, by = "grna_id")

# save the discovery pairs
saveRDS(discovery_pairs_extended, "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Gasperini/discovery_pairs.rds")
