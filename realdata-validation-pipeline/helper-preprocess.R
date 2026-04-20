# function for obtaining QC-ed data
obtain_qc_data <- function(path_to_outs_folder, library_index){
  
  # specify the particular h5 file of interest
  raw_count_file_path <- sprintf("%s/molecule_info.h5", path_to_outs_folder)
  qc_info_file_path <- sprintf("%s/filtered_feature_bc_matrix.h5", path_to_outs_folder)
  
  ###################### construct the raw data frame ##########################
  raw_count_file <- rhdf5::h5read(raw_count_file_path, "count")
  umi_idx <- rhdf5::h5read(raw_count_file_path, "umi")
  
  # obtain cell index
  barcode_idx <- rhdf5::h5read(raw_count_file_path, "barcode_idx")
  cell_barcodes <- rhdf5::h5read(raw_count_file_path, "barcodes")
  cell_idx <- cell_barcodes[barcode_idx + 1]
  
  # append the gem_group
  gem_group <- rhdf5::h5read(raw_count_file_path, "gem_group")
  cell_id_with_gem <- paste(cell_idx, gem_group, sep = "-")
  
  # obtain gene index for each RNA
  RNA_idx <- rhdf5::h5read(raw_count_file_path, "feature_idx")
  gene_reference <- rhdf5::h5read(raw_count_file_path, "features")
  gene_idx <- gene_reference$id[RNA_idx + 1]
  
  # store the data frame
  raw_data_frame <- data.frame(
    num_reads = raw_count_file,
    UMI_id = umi_idx + 1,
    cell_id = cell_id_with_gem,
    response_id = gene_idx
  )
  
  ############################ QC the raw data #################################
  qc_cell <- rhdf5::h5read(qc_info_file_path, "matrix/barcodes")
  
  # QC the raw data
  qc_df <- raw_data_frame |> dplyr::filter(cell_id %in% qc_cell)
  
  ############ Append the barcode with gem_group and library index #############
  qc_df <- qc_df |> dplyr::mutate(cell_id = stringi::stri_sub(paste(cell_id, library_index, sep = "_"), 1, -10))
    
  # return the reads vector
  return(qc_df)
}


# obtain the SRR file path
obtain_file_path <- function(experiment, type = "at_scale"){
  
  # obtain the folder path
  folder_path <- switch(type,
                        at_scale = {
                          switch (experiment,
                                  Gasperini = {
                                    setNames(c(paste0(.get_config_path("LOCAL_GASPERINI_2019_RAW_DATA_DIR"),
                                                      "processed/at_scale/run2/", sep = "/"),
                                               paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"),
                                                      "at-scale", sep = "/")),
                                             c("raw", "processed"))
                                  },
                                  Ray = {
                                    setNames(c(paste0(.get_config_path("LOCAL_RAY_2025_RAW_DATA_DIR"),
                                                      "raw", sep = "/"),
                                               paste0(.get_config_path("LOCAL_RAY_2025_RAW_DATA_DIR"),
                                                      "processed", sep = "/")),
                                             c("raw", "processed"))
                                  }
                          )
                        },
                        pilot = {
                          switch (experiment,
                                  Gasperini = {
                                    setNames(c(paste0(.get_config_path("LOCAL_GASPERINI_2019_RAW_DATA_DIR"),
                                                      "processed/pilot/", sep = "/"),
                                               paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"),
                                                      "high-MOI-pilot", sep = "/")),
                                             c("raw", "processed"))
                                  }
                          )
                        }
  )
  
  # return the folder path
  return(folder_path)
}

# function for computing the summary statistics
summary_data <- function(QC_data){
  
  # extract the number of total cells
  num_cells <- length(unique(QC_data$cell_id))
  
  # extract the number of reads per cell
  total_reads <- sum(QC_data$num_reads)
  num_reads_per_cell <- total_reads / num_cells
  
  # output the summary statistics
  return(
    setNames(c(num_cells, num_reads_per_cell), c("num_cells", "avg_reads"))
  )
}

# Use apply for submatrix divided into chunks
apply_in_chunks <- function(mat, chunk_size, margin, FUN, ...) {
  n <- if (margin == 1) nrow(mat) else ncol(mat)  # Total elements along the margin
  result <- vector("list", ceiling(n / chunk_size))  # Store results
  
  for (i in seq(1, n, by = chunk_size)) {
    end <- min(i + chunk_size - 1, n)
    submat <- if (margin == 1) mat[i:end, , drop = FALSE] else mat[, i:end, drop = FALSE]
    result[[ceiling(i / chunk_size)]] <- apply(submat, margin, FUN, ...)
  }
  
  do.call(c, result)  # Combine results
}