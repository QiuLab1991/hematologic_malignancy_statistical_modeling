suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(scTenifoldKnk)
})

seed <- 2025
qc_mtThreshold <- 0.1
qc_minLSize <- 50
nc_nNet <- 10
nc_nCells <- 500
nc_nComp <- 3

set.seed(seed)

input_file <- "data/scRNA_seq_data.rds"
output_dir <- "results/scrna"

candidate_genes <- c(
  "TCL1A", "CD7", "FAS", "CTSZ", "ICAM1",
  "PDCD1", "LTBR", "RNASET2", "ADM"
)

ko_genes <- c("FAS", "ICAM1", "LTBR")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scRNA <- readRDS(input_file)

required_metadata <- c("Sample.type", "cell_type")
missing_metadata <- setdiff(required_metadata, colnames(scRNA@meta.data))

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (!identical(rownames(scRNA@meta.data), Cells(scRNA))) {
  meta_fixed <- scRNA@meta.data

  if (nrow(meta_fixed) != length(Cells(scRNA))) {
    stop("Cell metadata and expression matrix contain different numbers of cells.")
  }

  prefix_cells <- sub("_.*", "", Cells(scRNA))
  prefix_meta <- sub("_.*", "", rownames(meta_fixed))

  if (!identical(prefix_cells, prefix_meta)) {
    stop("Cell metadata order does not match the expression matrix.")
  }

  rownames(meta_fixed) <- Cells(scRNA)
  scRNA@meta.data <- meta_fixed
}

scRNA@meta.data$cell_type <- recode(
  scRNA@meta.data$cell_type,
  "CD4+ NCM" = "Naive CD4+ T",
  "CD8+ CT" = "Cyt CD8+ T",
  "CD8+ NCM" = "Naive CD8+ T",
  "Cycling B" = "Cycling B",
  "Cycling T-cells" = "Cycling T",
  "Exhausted CD8+" = "Exh CD8+ T",
  "FCRL4+ B" = "FCRL4+ B",
  "FDC" = "FDC",
  "FL" = "Tumor B",
  "GC B" = "GC B",
  "Macrophages" = "Macrophages",
  "Memory B" = "Memory B",
  "NK" = "NK",
  "Naive B" = "Naive B",
  "Plasmablast" = "Plasma_cell",
  "TfH" = "Tfh",
  "Treg" = "Treg",
  "cDC1" = "cDC1",
  "mDC" = "mDC",
  "pDC" = "pDC"
)

candidate_genes <- candidate_genes[
  candidate_genes %in% rownames(scRNA)
]

if (length(candidate_genes) == 0) {
  stop("None of the candidate genes were found in the Seurat object.")
}

if ("umap" %in% names(scRNA@reductions)) {
  p_umap <- DimPlot(
    scRNA,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE
  )

  ggsave(
    file.path(output_dir, "umap_cell_types.pdf"),
    p_umap,
    width = 6,
    height = 5
  )
}

reference_object <- subset(
  scRNA,
  subset = Sample.type %in% c("FL", "Normal", "Reactive")
)

reference_object$Comparison <- ifelse(
  reference_object$Sample.type == "FL",
  "FL",
  "Reference"
)

expr_data <- FetchData(
  reference_object,
  vars = c("Comparison", candidate_genes)
)

comparison_results <- rbindlist(
  lapply(candidate_genes, function(gene) {

    x <- expr_data[[gene]]
    group <- expr_data$Comparison

    p_value <- tryCatch(
      wilcox.test(
        x ~ group,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )

    data.table(
      Gene = gene,
      FL_mean = mean(
        x[group == "FL"],
        na.rm = TRUE
      ),
      Reference_mean = mean(
        x[group == "Reference"],
        na.rm = TRUE
      ),
      P = p_value
    )
  })
)

comparison_results[, FDR := p.adjust(P, method = "BH")]

fwrite(
  comparison_results,
  file.path(
    output_dir,
    "FL_vs_reference_candidate_genes.csv"
  )
)

b_cell_types <- c(
  "Tumor B",
  "Naive B",
  "Memory B",
  "GC B",
  "FCRL4+ B",
  "Cycling B",
  "Plasma_cell"
)

fl_b <- subset(
  scRNA,
  subset = Sample.type == "FL" &
    cell_type %in% b_cell_types
)

Idents(fl_b) <- "cell_type"

b_cell_results <- FindAllMarkers(
  object = fl_b,
  features = candidate_genes,
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0
)

fwrite(
  as.data.table(b_cell_results),
  file.path(
    output_dir,
    "FL_B_cell_candidate_gene_expression.csv"
  )
)

tumor_b <- subset(
  scRNA,
  subset = Sample.type == "FL" &
    cell_type == "Tumor B"
)

countMatrix <- GetAssayData(
  tumor_b,
  assay = "RNA",
  slot = "counts"
)

missing_ko_genes <- setdiff(
  ko_genes,
  rownames(countMatrix)
)

if (length(missing_ko_genes) > 0) {
  stop(
    "KO genes not detected in the count matrix: ",
    paste(missing_ko_genes, collapse = ", ")
  )
}

if (ncol(countMatrix) < nc_nCells) {
  stop(
    "The number of tumor B cells is smaller than nc_nCells."
  )
}

for (gene in ko_genes) {

  set.seed(seed)

  result <- scTenifoldKnk(
    countMatrix = countMatrix,
    gKO = gene,
    qc = TRUE,
    qc_mtThreshold = qc_mtThreshold,
    qc_minLSize = qc_minLSize,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp
  )

  saveRDS(
    result,
    file.path(
      output_dir,
      paste0("scTenifoldKnk_", gene, ".rds")
    )
  )

  fwrite(
    as.data.table(result$diffRegulation),
    file.path(
      output_dir,
      paste0(
        "scTenifoldKnk_",
        gene,
        "_diffRegulation.csv"
      )
    )
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo.txt"
  )
)
