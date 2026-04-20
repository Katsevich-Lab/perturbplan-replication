# This is a Rscript for helper functions backing up the subsample analysis in run_subsample_analysis.sh
library(Matrix)
library(ondisc)

# obtain subsample for Gasperini data
construct_subsample_Gasperini <- function(plan_experiment, R, N, 
                                          gene_of_interest,
                                          grna_of_interest){
  
  # obtain the data folder
  path_to_data <- obtain_file_path(experiment = plan_experiment)
  
  # load gRNA count odm file
  grna_counts <- ondisc::read_odm(odm_fp = sprintf("%s/processed/grna_assignment/matrix.odm", path_to_data["processed"]), 
                                  metadata_fp = sprintf("%s/processed/grna_assignment/metadata.rds", path_to_data["processed"]))
  
  # extract the cells QC-ed
  cells_in_grna_counts <- stringi::stri_sub(rownames(grna_counts@cell_covariates), 1, -10)
  
  ############### sample over different SRR files ##############################
  # List all files in the directory with prefix "SRR"
  SRR_dirs_list <- list.files(path = path_to_data["raw"], pattern = "^SRR", full.names = TRUE)
  subsample_perSRR <- N / length(SRR_dirs_list)
  sampled_QC_data <- data.frame(NULL)
  gene_covariate <- data.frame(NULL)
  
  # sample N distinct cells
  for (SRR_dir in SRR_dirs_list) {
    
    # load the QC-ed data frame
    QC_data <- readRDS(sprintf("%s/processed/QC_data.rds", SRR_dir))
    
    # filter the cell in QC_data to ensure a subset of cells in processed data
    QC_data <- QC_data |> dplyr::filter(cell_id %in% cells_in_grna_counts)
    
    # sample subsample_perSRR cells and R * subsample_perSRR from QC_data
    sampled_cell <- sample(unique(QC_data$cell_id), subsample_perSRR)
    
    # subsample the reads
    sampled_reads <- QC_data |> 
      # sample over cells
      dplyr::filter(cell_id %in% sampled_cell) |>
      dplyr::mutate(SRR = basename(SRR_dir)) |>
      # sample over reads
      tidyr::uncount(num_reads) |>
      dplyr::slice_sample(n = subsample_perSRR * R) |>
      dplyr::distinct()
    
    # obtain the covariate from gene expression
    gene_covariate <- sampled_reads |>
      dplyr::group_by(cell_id) |>
      dplyr::summarise(
        library_size = dplyr::n(),
        num_expressed_gene = length(unique(response_id))
      ) |>
      dplyr::ungroup() |>
      dplyr::bind_rows(gene_covariate)
    
    # obtain the subsampled cell-gene data
    sampled_QC_data <- sampled_reads |>
      # form observed count df
      dplyr::filter(response_id %in% gene_of_interest) |>
      dplyr::group_by(cell_id, response_id) |>
      dplyr::summarise(obs_count = dplyr::n()) |>
      dplyr::ungroup() |>
      # bind the sampled data in previous rounds
      dplyr::bind_rows(sampled_QC_data)
    
    # release the memory by deleting QC_data
    rm(QC_data, sampled_reads, sampled_cell)
  }
  
  #################### construct gene expression matrix ########################
  # Define gene and cell order
  gene_levels <- gene_of_interest  
  cell_levels <- unique(sampled_QC_data$cell_id) 

  # Create index mappings
  cell_idx <- match(sampled_QC_data$cell_id, cell_levels)
  gene_idx <- match(sampled_QC_data$response_id, gene_levels)
  
  # Create sparse gene expression matrix
  response_matrix <- Matrix::sparseMatrix(
    i = gene_idx, 
    j = cell_idx, 
    x = sampled_QC_data$obs_count, 
    dims = c(length(gene_levels), length(cell_levels)), 
    dimnames = list(gene_levels, cell_levels)
  )
  rm(sampled_QC_data)
  
  ####################### construct gRNA count matrix ##########################
  # obtain the filtered grna barcode
  grna_idx <- match(grna_of_interest, rownames(grna_counts@feature_covariates)) 
  cell_idx <- match(cell_levels, cells_in_grna_counts)
  grna_matrix <- as(as(grna_counts[[grna_idx, cell_idx]], "generalMatrix"), "CsparseMatrix")
  rownames(grna_matrix) <- grna_of_interest
  colnames(grna_matrix) <- cell_levels
  
  # return the subsampled data
  return(list(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    gene_covariate = gene_covariate
  ))
}


# obtain subsample for Ray data
construct_subsample_Ray <- function(plan_experiment, R, N, 
                                    gene_of_interest,
                                    grna_of_interest,
                                    grna_threshold = 5,
                                    sample_summary){
  
  # obtain the data folder
  path_to_data <- obtain_file_path(experiment = plan_experiment)
  
  ############### sample over different SRR files ##############################
  # List all files in the directory with prefix "SRR"
  sample_dirs_list <- list.files(path = paste0(path_to_data["raw"], "/cell_ranger/"), pattern = "230327", full.names = TRUE)
  sampled_gene_data <- data.frame(NULL)
  sampled_grna_data <- data.frame(NULL)
  gene_covariate <- data.frame(NULL)
  
  # load target_panel
  target_panel <- readr::read_csv(sprintf("%s/target_panel.csv", sample_dirs_list[1]), 
                                  skip = 5, col_names = TRUE) |> 
    dplyr::select(gene_id) |> 
    dplyr::distinct() |> 
    dplyr::pull()
  
  # decide the sampling budge for cells
  sample_summary <- as.data.frame(sample_summary) |> 
    dplyr::mutate(num_total_samples = sum(num_cells)) |>
    dplyr::mutate(sampled_cells = round(num_cells * N / num_total_samples))
  
  # sample N distinct cells
  for (sample_dir in sample_dirs_list) {
    
    # load the QC-ed data frame
    QC_data <- readRDS(sprintf("%s/processed/QC_data.rds", sample_dir))
    
    # extract the number of subsamples
    num_subsample <- sample_summary[basename(sample_dir), "sampled_cells"]
    
    # sample num_subsample cells and R * num_subsample from QC_data
    sampled_cell <- sample(unique(QC_data$cell_id), num_subsample)
    
    # sample cells 
    sampled_cells <- QC_data |> dplyr::filter(cell_id %in% sampled_cell)
    
    # extract the gRNA
    sampled_grna_data <- sampled_cells |>
      # form observed count df
      dplyr::filter(response_id %in% grna_of_interest) |>
      dplyr::group_by(cell_id, response_id) |>
      dplyr::summarise(obs_count = dplyr::n()) |>
      dplyr::ungroup() |>
      # bind the sampled data in previous rounds
      dplyr::bind_rows(sampled_grna_data)
    
    # subsample the reads from sampled_cells
    sampled_reads <- sampled_cells[startsWith(sampled_cells$response_id, "ENSG"), ] |> 
      # sample over reads
      tidyr::uncount(num_reads) |>
      dplyr::slice_sample(n = num_subsample * R) |>
      dplyr::distinct()
    
    # obtain the covariate from gene expression
    gene_covariate <- sampled_reads |>
      dplyr::filter(response_id %in% target_panel) |>
      dplyr::group_by(cell_id) |>
      dplyr::summarise(
        library_size = dplyr::n(),
        num_expressed_gene = length(unique(response_id))
      ) |>
      dplyr::ungroup() |>
      dplyr::bind_rows(gene_covariate)
    
    # obtain the subsampled cell-gene data
    sampled_gene_data <- sampled_reads |>
      # form observed count df
      dplyr::filter(response_id %in% gene_of_interest) |>
      dplyr::group_by(cell_id, response_id) |>
      dplyr::summarise(obs_count = dplyr::n()) |>
      dplyr::ungroup() |>
      # bind the sampled data in previous rounds
      dplyr::bind_rows(sampled_gene_data)
    
    # release the memory by deleting QC_data
    rm(QC_data, sampled_reads, sampled_cell)
  }
  
  # extract the common cell indices
  common_cell_idx <- intersect(unique(sampled_grna_data$cell_id), unique(sampled_gene_data$cell_id))
  
  #################### Create sparse gene expression matrix ####################
  # filter the sampled gene data according to the common cell_idx
  sampled_gene_data <- sampled_gene_data |> dplyr::filter(cell_id %in% common_cell_idx)
  
  # create x and y axes mappings
  cell_levels <-  unique(sampled_gene_data$cell_id)
  gene_levels <- gene_of_interest
  cell_idx <- match(sampled_gene_data$cell_id, cell_levels)
  gene_idx <- match(sampled_gene_data$response_id, gene_levels)
  
  # create matrix
  response_matrix <- Matrix::sparseMatrix(
    i = gene_idx, 
    j = cell_idx, 
    x = sampled_gene_data$obs_count, 
    dims = c(length(gene_levels), length(cell_levels)), 
    dimnames = list(gene_levels, cell_levels)
  )
  rm(sampled_gene_data)
  
  ################### Create sparse grna assignment matrix #####################
  # filter the sampled grna data according to the common cell_idx
  sampled_grna_data <- sampled_grna_data |> dplyr::filter(cell_id %in% common_cell_idx)
  
  # create x and y axes mappings
  cell_levels <-  unique(sampled_grna_data$cell_id)
  grna_levels <- grna_of_interest
  cell_idx <- match(sampled_grna_data$cell_id, cell_levels)
  grna_idx <- match(sampled_grna_data$response_id, grna_levels)
  
  # do gRNA assignment
  sampled_grna_data <- sampled_grna_data |> 
    dplyr::mutate(grna_assignment = ifelse(obs_count >= grna_threshold, 1, 0))
  
  # create sparse matrix
  grna_matrix <- Matrix::sparseMatrix(
    i = grna_idx, 
    j = cell_idx, 
    x = sampled_grna_data$grna_assignment, 
    dims = c(length(grna_levels), length(cell_levels)), 
    dimnames = list(grna_levels, cell_levels)
  )
  rm(sampled_grna_data)
  
  # return the subsampled data
  return(list(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    gene_covariate = gene_covariate |> dplyr::filter(cell_id %in% common_cell_idx)
  ))
}


