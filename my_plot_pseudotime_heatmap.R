my_plot_pseudotime_heatmap <- function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2", 
  num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
  add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
  norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3, 
  trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE, 
  cores = 1) 
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
    max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
    relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
    FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
    cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
    clustering_distance_rows = row_dist, clustering_method = hclust_method, 
    cutree_rows = num_clusters, silent = TRUE, filename = NA, 
    breaks = bks, border_color = NA, color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
      num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
      ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
        "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
        "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if (!is.null(annotation_row)) 
      row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row)) 
    row.names(annotation_row) <- row_ann_labels
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
    cluster_rows = cluster_rows, show_rownames = show_rownames, 
    show_colnames = F, clustering_distance_rows = row_dist, 
    clustering_method = hclust_method, cutree_rows = num_clusters, 
    annotation_row = annotation_row, annotation_col = annotation_col, 
    treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
    border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  resaa=list()
  resaa$annotation_row = annotation_row
  resaa$ph_res = ph_res
  if (return_heatmap) {
    return(resaa)
  }
}