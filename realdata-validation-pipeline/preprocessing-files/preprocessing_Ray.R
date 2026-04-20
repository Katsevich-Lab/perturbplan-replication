# This is a Rscript preparing for data preprocessing for Gasperini data
library(sceptre)

# source the helper function
source("~/.Rprofile")
source("reproduction-code-prep/realdata-validation-pipeline/helper-preprocess.R")

# create preprocessing directory
preprocessing_dir <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Ray"
dir.create(preprocessing_dir, recursive = TRUE)

##################### 1. Preprocess the at-scale data ##########################
# preprocess the molecule_info.h5 for gene expression
path_to_plan_data <- obtain_file_path(experiment = "Ray")
sample_dirs_list <- list.files(path = sprintf("%s/cell_ranger", path_to_plan_data["raw"]), 
                               pattern = "230327",
                               full.names = TRUE)
summary_persample <- matrix(0, nrow = length(sample_dirs_list), ncol = 2,
                            dimnames = list(
                              sample = basename(sample_dirs_list),
                              Summary = c("num_cells", "avg_reads")
                            ))

# loop over different SRRs
for (sample_dir in sample_dirs_list) {
  
  # make a director under this path
  path_to_save <- sprintf("%s/processed", sample_dir)
  if(!dir.exists(path_to_save)){
    dir.create(path_to_save, recursive = TRUE)
  }
  
  # if the file is already there then the QC should be performed or not
  if(file.exists(sprintf("%s/QC_data.rds", path_to_save))){
    QC_data <- readRDS(sprintf("%s/QC_data.rds", path_to_save))
  }else{
    # process the raw data
    QC_data <- obtain_qc_data(path_to_outs_folder = sprintf("%s", sample_dir),
                              library_index = basename(sample_dir))
    
    # save the QC-ed data
    saveRDS(QC_data, file = sprintf("%s/QC_data.rds", path_to_save))
  }
  
  # extract the summary statistics
  summary_persample[basename(sample_dir), ] <- summary_data(QC_data[startsWith(QC_data$response_id, "ENSG"), ])
  
  # print the summary
  print(summary_persample[basename(sample_dir), ])
  
  # delete the variable
  rm(QC_data)
}

# save the summary information
summary_sample_path <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Ray/summary_SRR.rds"
summary_persample_path <- "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Ray/summary_persample.rds"
if(!file.exists(summary_sample_path)){
  summary_sample <- list(
    minimum_num_cells = min(summary_persample[, "num_cells"]),
    minimum_num_reads_per_cell = min(summary_persample[, "avg_reads"]),
    num_SRRs = length(sample_dirs_list)
  )
  saveRDS(summary_sample, summary_sample_path)
  saveRDS(summary_persample, summary_persample_path)
}

###################### 2. Construct the discovery pairs ########################
set.seed(1)

# compute the median of guides per target
guides_per_target <- 15

# load sceptre object
sceptre_object <- readRDS(paste0(path_to_plan_data["raw"], "/sceptre_output/final_sceptre_object.rds"))

# load the reference gene information
feature_file_dir <- paste0(sample_dirs_list[1], "/features.tsv")
if(!file.exists(feature_file_dir)){
  R.utils::gunzip(paste0(sample_dirs_list[1], "/features.tsv.gz"), remove = FALSE)  
}
gene_name_to_id <- readr::read_tsv(feature_file_dir, col_names = FALSE)
colnames(gene_name_to_id) <- c("response_id", "gene_name", "type")

# get the gRNA threshold
threshold <- sceptre_object@grna_assignment_hyperparameters$threshold
grna_matrix <- sceptre_object@grna_matrix[[1]]
response_list <- rownames(sceptre_object@response_matrix[[1]])

########################### 2.1 get the positive pairs #########################
# filter the target which has exactly 15 guides
tss_target_guides <- sceptre_object@grna_target_data_frame[grepl("TSS", sceptre_object@grna_target_data_frame$grna_target), ]
filterd_target <- tss_target_guides |> 
  dplyr::group_by(grna_target) |> 
  dplyr::summarise(num_guides = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::filter(num_guides == guides_per_target) |> 
  dplyr::select(grna_target) |>
  dplyr::pull()
tss_target_guides <- tss_target_guides |> 
  dplyr::filter(grna_target %in% filterd_target) |>
  dplyr::rowwise() |>
  dplyr::mutate(num_oracle_cells = sum(grna_matrix[grna_id, ] >= threshold)) |> 
  dplyr::ungroup()

# extract the target-gene information
positive_target_gene <- sceptre_object@discovery_pairs |> 
  dplyr::filter(grna_target %in% unique(tss_target_guides$grna_target)) |> 
  dplyr::mutate(target_gene = stringr::str_replace(grna_group, "_.*", "")) |>
  dplyr::left_join(gene_name_to_id, by = "response_id") |>
  dplyr::filter(target_gene == gene_name) |> 
  dplyr::select(grna_target, response_id) |> 
  dplyr::mutate(target_type = "gene_tss")

# combine information 
positive_pairs <- positive_target_gene |> 
  dplyr::left_join(tss_target_guides, by = "grna_target") |> 
  dplyr::select(target_type, response_id, grna_id, grna_target, num_oracle_cells)

########################### 2.2 get the negative pairs #########################
# extract the non-tageting guides
non_target_guides <- sceptre_object@grna_target_data_frame[grepl("non-targeting", sceptre_object@grna_target_data_frame$grna_target), ] 

# construct non-targeting dataframe
genes_per_target <- 100
num_targets <- round(nrow(non_target_guides) / guides_per_target)
negative_target_gene <- data.frame(
  response_id = unlist(lapply(1:num_targets, function(x) sample(response_list, genes_per_target))),
  grna_target = rep(sprintf("random_%d", 1:num_targets), each = genes_per_target),
  target_type = "non-targeting"
)

# extract the number of cells per guide
non_target_guides <- non_target_guides |> 
  dplyr::rowwise() |>
  dplyr::mutate(num_oracle_cells = sum(grna_matrix[grna_id, ] >= threshold)) |>
  dplyr::ungroup() |>
  dplyr::mutate(grna_target = sprintf("random_%d", pmin(ceiling(seq_len(nrow(non_target_guides)) / guides_per_target), num_targets)))

# combine the information
negative_pairs <- negative_target_gene |> 
  dplyr::left_join(non_target_guides, by = "grna_target") |> 
  dplyr::select(target_type, response_id, grna_id, grna_target, num_oracle_cells)

# combine positive and negative control pairs
discovery_pairs <- dplyr::bind_rows(positive_pairs, negative_pairs) |> 
  dplyr::mutate(num_total_plan_cells = ncol(grna_matrix))

# save the discovery pairs
saveRDS(discovery_pairs, "reproduction-code-prep/realdata-validation-pipeline/intermediate-files/Ray/discovery_pairs.rds")
