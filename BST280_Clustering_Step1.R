library(Biobase)
library(pheatmap)
library(SummarizedExperiment)
library(dendextend)

filtered_eset <- readRDS("luad_gtex_eset_1207.rds")

### 1. Extract expression matrix and gene names
expr_mat <- assayData(filtered_eset)[["norm_logtpm"]]
expr_mat <- expr_mat[apply(expr_mat, 1, var) > 0, ]        # drop zero-variance genes
gene_names <- fData(filtered_eset)$gene_name
rownames(expr_mat) <- gene_names

### 2. Scale genes (Z-score across samples)
expr_scaled <- t(scale(t(expr_mat)))

### 3. Sample + gene distances for clustering
dist_samples <- as.dist(1 - cor(expr_scaled, method = "spearman"))
dist_genes   <- as.dist(1 - cor(t(expr_scaled), method = "pearson"))

hc_samples <- hclust(dist_samples, method = "complete")
hc_genes   <- hclust(dist_genes,   method = "complete")

### 4. Heatmap (no annotation)
pheatmap(
  expr_scaled,
  cluster_rows = hc_genes,
  cluster_cols = hc_samples,
  show_rownames = FALSE,
  fontsize_col = 8,
  show_colnames = FALSE,
  main = "Hierarchical Clustering: LUAD vs GTEx (All Genes)"
)

### 5. Annotation (Dataset: TCGA_LUAD vs GTEX_LUNG)
annotation <- data.frame(
  Dataset = pData(filtered_eset)$dataset
)
rownames(annotation) <- colnames(expr_scaled)

ann_colors <- list(
  Dataset = c("TCGA_LUAD" = "red", "GTEx_LUNG" = "blue")
)

### 6. Heatmap with annotations
pheatmap(
  expr_scaled,
  cluster_rows = hc_genes,
  cluster_cols = hc_samples,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Hierarchical Clustering with Dataset Annotation"
)

geneTree <- hclust(dist_genes, method = "average")
geneDend <- as.dendrogram(geneTree)
hclusth0.5 <- cutree(geneTree, h = 5)
hclusth1.0 <- cutree(geneTree, h = 1)
hclusth1.5 <- cutree(geneTree, h = 1.5)

plot(
  geneDend,
  leaflab = "none",
  main = "Gene Clustering Dendrogram",
  ylab = "Height"
)

# Add the colored cluster bars under the dendrogram
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)
colored_bars(
  the_bars,
  geneDend,
  sort_by_labels_order = TRUE,   # ensure alignment with dendrogram order
  y_shift = -0.1,
  rowLabels = c("h = 10", "h = 5", "h = -5","h = -10"),
  cex.rowLabels = 0.7
)
abline(h = 1.5, lty = 2, col = "grey")
abline(h = 1.0, lty = 2, col = "grey")
abline(h = 0.5, lty = 2, col = "grey")
