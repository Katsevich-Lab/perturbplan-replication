library(data.table)
library(rhdf5)
library(Seurat)

get_number_of_genes <- function(obj) {
  length(obj$baseline_expression_stats$response_id)
}

get_cells_and_reads <- function(name) {
  config_name <- reference_expression_datasets |>
    dplyr::filter(dataset_name == name) |>
    dplyr::pull(config_name)
  process_function <- reference_expression_datasets |>
    dplyr::filter(dataset_name == name) |>
    dplyr::pull(process_function)
  path <- .get_config_path(config_name)
  get(process_function)(path)
}

get_UMI_saturation <- function(obj) {
  obj$library_parameters$UMI_per_cell_at_saturation
}

get_UMI_variation <- function(obj) {
  obj$library_parameters$variation
}

get_mapping_efficiency <- function(obj) {
  obj$mapping_efficiency
}

# ---- one unified "compiler" for ALL strings ----
canon <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)   # non-breaking space -> space
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x)              # drop spaces/punct/underscores etc.
  x
}

compile_to_actual <- function(desired, available, what = "strings") {
  desired <- as.character(desired)
  available <- as.character(available)

  amap <- setNames(available, canon(available))
  out <- unname(amap[canon(desired)])

  if (anyNA(out)) {
    miss <- desired[is.na(out)]
    stop(
      sprintf("Unmatched %s: %s", what, paste(miss, collapse = ", ")),
      "\nAvailable options are: ", paste(available, collapse = ", ")
    )
  }
  out
}

# ---- plotting helper (one panel per metric) ----
make_metric_plot <- function(metric_name) {
  ggplot(long_dt[metric == metric_name], aes(x = `Cell type`, y = value, fill = `Cell type`)) +
    geom_col(width = 0.8) +
    labs(x = NULL, y = NULL, title = metric_name, fill = "Cell type") +
    theme_classic() +
    scale_fill_manual(values = manual_colors) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.position = "bottom"
    )
}

# ---- Processing functions ----
process_k562_gasperini <- function(path_to_dataset) {
  message("Start processing K562_Gasperini")

  path_to_run1 <- file.path(path_to_dataset, "processed", "at_scale", "run1")
  path_to_run2 <- file.path(path_to_dataset, "processed", "at_scale", "run2")

  k562_data_1 <- perturbplan::reference_data_preprocessing_10x(path_to_run1)
  response_matrix_1 <- k562_data_1[[1]]
  read_umi_table    <- k562_data_1[[2]]

  k562_data_2 <- perturbplan::reference_data_preprocessing_10x(path_to_run2)
  response_matrix_2 <- k562_data_2[[1]]

  response_matrix <- cbind(response_matrix_1, response_matrix_2)

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}

process_k562_10x <- function(path_to_dataset) {
  message("Start processing K562_10x")

  k562_data <- perturbplan::reference_data_preprocessing_10x(path_to_dataset)
  response_matrix <- k562_data[[1]]
  read_umi_table  <- k562_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}

process_t_cd8_10x <- function(path_to_dataset) {
  message("Start processing T_CD8_Shifrut")

  path_to_runs <- file.path(path_to_dataset, "processed")
  t_cd8_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)

  response_matrix <- t_cd8_data[[1]]
  read_umi_table  <- t_cd8_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}

process_a549_10x <- function(path_to_dataset) {
  message("Start processing A549_10x")

  path_to_runs <- file.path(path_to_dataset, "processed")
  a549_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)

  response_matrix <- a549_data[[1]]
  read_umi_table  <- a549_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}

process_ipsc_10x <- function(path_to_dataset) {
  message("Start processing iPSC_10x")

  ipsc_data <- perturbplan::reference_data_preprocessing_10x(path_to_dataset)
  response_matrix <- ipsc_data[[1]]
  read_umi_table  <- ipsc_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}
process_ipsc_neuron_10x <- function(path_to_dataset) {
  message("Start processing iPSC_neuron_10x")

  ipsc_neuron_data <- perturbplan::reference_data_preprocessing_10x(path_to_dataset)
  response_matrix <- ipsc_neuron_data[[1]]
  read_umi_table  <- ipsc_neuron_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}

process_thp1_10x <- function(path_to_dataset) {
  message("Start processing THP1_Yao")

  path_to_runs <- file.path(path_to_dataset, "processed")
  dir_srrs <- list.dirs(path_to_runs, recursive = FALSE, full.names = TRUE)

  seurat_obj <- readRDS(file.path(path_to_runs, "GSM6858447_KO_conventional.rds"))

  nt_cells <- colnames(seurat_obj[["perturbations"]])[which(
    grepl("non-targeting|safe-targeting", seurat_obj@meta.data$Guides) &
      seurat_obj@meta.data$Total_number_of_guides == 1
  )]

  nt_cells_once <- nt_cells[which(sub(".*-(.*)$", "\\1", nt_cells) == 1)]

  response_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  num_all_cells <- ncol(response_matrix)
  response_matrix <- response_matrix[, nt_cells_once]

  if (any(response_matrix < 0)) {
    stop("Response matrix contains negative values.")
  }

  genes <- fread(file.path(dir_srrs[1], "outs", "filtered_feature_bc_matrix", "features.tsv.gz"),
                 header = FALSE)
  colnames(genes) <- c("gene_id", "gene_name", "gene_type")

  gene_map <- genes[!duplicated(gene_name), .(gene_name, gene_id)]
  response_matrix <- response_matrix[rownames(response_matrix) %in% gene_map$gene_name, ]

  setkey(gene_map, gene_name)
  gene_ids <- gene_map[rownames(response_matrix), gene_id]
  rownames(response_matrix) <- gene_ids

  srr_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)
  read_umi_table <- srr_data[[2]]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads)/num_all_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}


process_k562_tap <- function(path_to_dataset) {
  message("Start processing K562_TAP")

  path_to_runs <- file.path(path_to_dataset, "processed")
  k562_tap_data <- perturbplan::reference_data_preprocessing_10x(path_to_runs)

  response_matrix <- k562_tap_data[[1]]
  read_umi_table  <- k562_tap_data[[2]]

  panel_path <- file.path(path_to_runs, "perturbplan-demo", "outs", "target_panel.csv")
  if (!file.exists(panel_path)) {
    stop("target_panel.csv not found at: ", panel_path)
  }

  .trim <- function(x) gsub("^\\s+|\\s+$", "", x)

  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(panel_path, header = FALSE, skip = 6, select = 1L, data.table = FALSE)
    gene_list <- unique(.trim(dt[[1]]))
  } else {
    con <- file(panel_path, open = "r"); on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    first_col <- vapply(strsplit(lines[-(1:6)], ",", fixed = TRUE),
                        function(z) z[[1]], character(1L))
    gene_list <- unique(.trim(first_col))
  }

  gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]
  if (length(gene_list) == 0L) {
    stop("Parsed gene_list is empty.")
  }

  response_matrix <- response_matrix[rownames(response_matrix) %in% gene_list, , drop = FALSE]

  num_cells <- ncol(response_matrix)
  umi_per_cell <- mean(colSums(response_matrix))
  sequenced_reads_per_cell <- sum(read_umi_table$num_reads[read_umi_table$response_id %in% gene_list])/num_cells

  list(
    num_cells = num_cells,
    umi_per_cell = umi_per_cell,
    sequenced_reads_per_cell = sequenced_reads_per_cell
  )
}
