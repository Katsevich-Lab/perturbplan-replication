library(perturbplan)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)
library(gridExtra)
source("reproduction-code-prep/final_plots/Extended_figure_1/helper.R")

# load data
data("K562_Gasperini")
data("K562_10x")
data("K562_Ray")
data("A549_Sakellaropoulos")
data("THP1_Yao")
data("T_CD8_Shifrut")
data("iPSC_Tian")
data("iPSC_neuron_Tian")
data("reference_expression_datasets")

## 1. Create dataset information table for later use

# Create full list of all datasets for generating cells_and_reads_summary
objs <- list(
  K562_Gasperini       = K562_Gasperini,
  K562_10x             = K562_10x,
  K562_Ray             = K562_Ray,
  A549_Sakellaropoulos = A549_Sakellaropoulos,
  THP1_Yao             = THP1_Yao,
  T_CD8_Shifrut        = T_CD8_Shifrut,
  iPSC_Tian            = iPSC_Tian,
  iPSC_neuron_Tian     = iPSC_neuron_Tian
)

# Generate or load the cells and reads summary
rds_file <- "reproduction-code-prep/final_plots/Extended_figure_1/cells_and_reads_summary.rds"
if (file.exists(rds_file)) {
  message("Loading existing cells_and_reads_summary.rds")
  cells_and_reads <- readRDS(rds_file)
} else {
  message("Generating cells_and_reads_summary.rds (this may take a while)...")
  cells_and_reads <- lapply(
    names(objs),
    get_cells_and_reads
  )
  saveRDS(cells_and_reads, rds_file)
  message("Successfully saved cells_and_reads_summary.rds")
}

# Create list of datasets used in paper
objs_paper <- list(
  K562_Gasperini       = K562_Gasperini,
  A549_Sakellaropoulos = A549_Sakellaropoulos,
  THP1_Yao             = THP1_Yao,
  T_CD8_Shifrut        = T_CD8_Shifrut,
  iPSC_Tian            = iPSC_Tian,
  iPSC_neuron_Tian     = iPSC_neuron_Tian
)

# extract num_cells
num_cells <- vapply(
  cells_and_reads,
  function(x) as.numeric(x$num_cells),
  numeric(1)
)

# extract umi_per_cell
umi_per_cell <- vapply(
  cells_and_reads,
  function(x) as.numeric(x$umi_per_cell),
  numeric(1)
)

# extract sequenced_reads_per_cell
sequenced_reads_per_cell <- vapply(
  cells_and_reads,
  function(x) as.numeric(x$sequenced_reads_per_cell),
  numeric(1)
)
# construct dataset
dataset_information_df_paper <- data.frame(
  `Cell type` = c("K562","A549","THP-1","CD8+ T cells","iPSC","iPSC-derived neurons"),
  `Number of Genes Captured` = vapply(objs_paper, get_number_of_genes, numeric(1)),
  `Number of Cells in Learning` = num_cells[c(1,4:length(num_cells))],
  `Sequenced Reads per Cell` = sequenced_reads_per_cell[c(1,4:length(num_cells))],
  `Mapping Efficiency` = vapply(objs_paper, get_mapping_efficiency, numeric(1)),
  `UMIs per Cell` = umi_per_cell[c(1,4:length(num_cells))],
  `UMIs per Cell at Saturation` = vapply(objs_paper, get_UMI_saturation, numeric(1)),

  stringsAsFactors = FALSE,
  row.names = NULL
)

## 2. Visualize paper information as bar plots ----
manual_colors <- c(
  "A549" = "#b2182b",
  "iPSC"            = "#d6604d",
  "iPSC-derived neurons"     = "#f4a582",
  "K562"       = "#d1e5f0",
  "CD8+ T cells"        = "#2166ac",
  "THP-1"             = "#fddbc7"
)

# ---- input ----
df <- dataset_information_df_paper
dt <- as.data.table(df)

# ---- robustly find the "Cell type" column (compiled) ----
cell_col <- compile_to_actual("Cell type", names(dt), what = "id column")[1]

# ---- user-specified metric order (human-readable) ----
metric_order_raw <- c(
  "Number of genes captured",
  "Number of cells",
  "Sequenced reads per cell",
  "Mapping efficiency",
  "UMIs per cell",
  "UMIs per cell at saturation"
)

# rename dt
old_names <- names(dt)
metrics_old <- setdiff(old_names, "Cell.type")
stopifnot(length(metrics_old) == length(metric_order_raw))
setnames(
  dt,
  old = metrics_old,
  new = metric_order_raw
)

# compile desired metric order to the *actual* column names in dt
metric_order <- compile_to_actual(metric_order_raw, setdiff(names(dt), cell_col), what = "metrics")

# ---- melt to long ----
long_dt <- melt(
  dt,
  id.vars = cell_col,
  measure.vars = metric_order,
  variable.name = "metric",
  value.name = "value"
)

# apply the same compiler to cell type values (for stable factor levels)
cell_levels <- unique(long_dt[[cell_col]])
cell_levels_compiled <- canon(cell_levels)
# keep original display labels, but order by compiled keys for stability
cell_levels <- cell_levels[order(cell_levels_compiled)]

setnames(long_dt, cell_col, "Cell type")
long_dt[, `Cell type` := factor(`Cell type`, levels = cell_levels)]
long_dt[, metric := factor(as.character(metric), levels = metric_order)]

# plot different metrics
metric_plots <- lapply(metric_order, make_metric_plot)

# ---- combine panels in the specified order ----
Figure_metrics <- wrap_plots(metric_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18)
  )

# ---- save (set plots_dir if you use it elsewhere) ----
plots_dir <- "reproduction-code-prep/final_plots/figures/"
ggsave(file.path(plots_dir, "extended_figure_1.pdf"), Figure_metrics, width = 14, height = 16)
