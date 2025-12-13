library(SummarizedExperiment)
library(Biobase)
library(limma)
library(biomaRt)
library(ggfortify)


list.files("data")

#Load TCGA lung adenocarcinoma file and create an eset dataset
#Here we will set logtpm to be the exprs level for now to run differential expression later
#We will also append raw counts and normalized value into the eset dataset
luad <- readRDS("data/tcga_luad.rds")
raw_count_luad <- assay(luad,"raw_counts")
normalize_luad <- assay(luad,"logtpm")
clin_luad <-as.data.frame(colData(luad))
gene_annot_luad <- as.data.frame(rowData(luad))

eset_luad <- ExpressionSet(
  assayData = assayDataNew(
    exprs   = as.matrix(normalize_luad),
    counts  = as.matrix(raw_count_luad),
    log_tpm = as.matrix(normalize_luad)
  ),
  phenoData    = AnnotatedDataFrame(clin_luad),
  featureData  = AnnotatedDataFrame(gene_annot_luad)
)

#Load GTEX normal lung tissue file and create an eset dataset
gtex_lung <- readRDS("data/gtex_lung.rds")
raw_count_gtex_lung  <- assay(gtex_lung ,"raw_counts")
normalize_gtex_lung <- assay(gtex_lung ,"logtpm")
clin_gtex_lung <-as.data.frame(colData(gtex_lung ))
gene_annot_gtex_lung <- as.data.frame(rowData(gtex_lung ))

eset_gtex_lung <- ExpressionSet(
  assayData = assayDataNew(
    exprs   = as.matrix(normalize_gtex_lung),
    counts  = as.matrix(raw_count_gtex_lung),
    log_tpm = as.matrix(normalize_gtex_lung)
  ),
  phenoData    = AnnotatedDataFrame(clin_gtex_lung),
  featureData  = AnnotatedDataFrame(gene_annot_gtex_lung)
)

#Harmonizing the two datasets tcga_luad and gtex_lung
#We will extract exprs, assayData, pData(clinical data) and fData(gene annotation)
#We will then find overlaping columns and combine these two dataset by outer join
#We will duplicate the new individual eset files first
eset_luad_H  <- eset_luad
eset_gtex_H  <- eset_gtex_lung
mat_expr_luad   <- exprs(eset_luad)
mat_expr_gtex   <- exprs(eset_gtex_lung)
mat_counts_luad <- assayData(eset_luad)$counts
mat_counts_gtex <- assayData(eset_gtex_lung)$counts
mat_logtpm_luad <- assayData(eset_luad)$log_tpm
mat_logtpm_gtex <- assayData(eset_gtex_lung)$log_tpm
pd_luad <- pData(eset_luad)
pd_gtex <- pData(eset_gtex_lung)
fd_luad <- fData(eset_luad)
fd_gtex <- fData(eset_gtex_lung)

#Find intersect fData and manipulate individual eset accordingly
common_fd_cols <- intersect(colnames(fd_luad), colnames(fd_gtex))
unique_fd_luad <- setdiff(colnames(fd_luad), colnames(fd_gtex))
unique_fd_gtex <- setdiff(colnames(fd_gtex), colnames(fd_luad))
for (col in unique_fd_luad) fd_gtex[[col]] <- NA
for (col in unique_fd_gtex) fd_luad[[col]] <- NA
fd_luad <- fd_luad[, sort(colnames(fd_luad))]
fd_gtex <- fd_gtex[, sort(colnames(fd_gtex))]
fData(eset_luad_H) <- fd_luad
fData(eset_gtex_H) <- fd_gtex

#Find intersect pData and manipulate individual eset accordingly
common_pd_cols <- intersect(colnames(pd_luad), colnames(pd_gtex))
unique_pd_luad <- setdiff(colnames(pd_luad), colnames(pd_gtex))
unique_pd_gtex <- setdiff(colnames(pd_gtex), colnames(pd_luad))
for (col in unique_pd_luad) pd_gtex[[col]] <- NA
for (col in unique_pd_gtex) pd_luad[[col]] <- NA
pd_luad <- pd_luad[, sort(colnames(pd_luad))]
pd_gtex <- pd_gtex[, sort(colnames(pd_gtex))]
pData(eset_luad_H) <- pd_luad
pData(eset_gtex_H) <- pd_gtex

#find interecting expressed gene and subset the eset files 
common_genes <- intersect(rownames(eset_luad_H), rownames(eset_gtex_H))
esetH_luad <- eset_luad_H[common_genes, ]
esetH_gtex <- eset_gtex_H[common_genes, ]

#combine the assayData, pdata, fdata into a new combined dataset
mat_expr_combined <- cbind(exprs(esetH_luad), exprs(esetH_gtex))
mat_counts_combined <- cbind(assayData(esetH_luad)$counts,assayData(esetH_gtex)$counts)
mat_logtpm_combined <- cbind(assayData(esetH_luad)$log_tpm,assayData(esetH_gtex)$log_tpm)
pData(esetH_luad)$dataset <- "TCGA_LUAD"
pData(esetH_gtex)$dataset <- "GTEx_LUNG"
pd_combined <- rbind(pData(esetH_luad), pData(esetH_gtex))
fd_combined <- fData(esetH_luad)

combined_eset <- ExpressionSet(
  assayData = assayDataNew(
    exprs   = mat_expr_combined,
    counts  = mat_counts_combined,
    log_tpm = mat_logtpm_combined
  ),
  phenoData   = AnnotatedDataFrame(pd_combined),
  featureData = AnnotatedDataFrame(fd_combined)
)

validObject(combined_eset)

#Create counts per million assayData element 'tpm' and filter accordingly
#Filter criteria 1: average row mean for tpm greater than 1 is > 0.2
#Filter criteria 2: average row mean for raw counts greater than 10 is > 50%
counts_luad <- assayData(eset_luad)$counts
counts_gtex <- assayData(eset_gtex_lung)$counts
tpm_luad <- 2^(assayData(eset_luad)$log_tpm) - 1
tpm_gtex <- 2^(assayData(eset_gtex_lung)$log_tpm) - 1
common_genes <- intersect(rownames(counts_luad), rownames(counts_gtex))
counts_luad2 <- counts_luad[common_genes, ]
counts_gtex2 <- counts_gtex[common_genes, ]
tpm_luad2 <- tpm_luad[common_genes, ]
tpm_gtex2 <- tpm_gtex[common_genes, ]
counts_all <- cbind(counts_luad2, counts_gtex2)
tpm_all    <- cbind(tpm_luad2,   tpm_gtex2)
tpm_filter   <- rowMeans(tpm_all > 1) >= 0.20
count_filter <- rowMeans(counts_all > 10) >= 0.50
keep_genes <- tpm_filter & count_filter
sum(keep_genes)
filtered_eset <- combined_eset[keep_genes, ]
validObject(filtered_eset)
dim(filtered_eset)

#Correct for batch effect, we will skip quantile normalization for this step
logtpm <- exprs(filtered_eset)
batch <- pData(filtered_eset)$dataset   # “TCGA_LUAD” / “GTEx_Lung”
logtpm_corrected <- removeBatchEffect(
  logtpm,                   # or logtpm if not quantile normalized
  batch = batch
)
assayData(filtered_eset) <- assayDataNew(
  exprs      = logtpm_corrected,                # replace main slot
  counts     = assayData(filtered_eset)$counts,
  log_tpm    = assayData(filtered_eset)$log_tpm,
  norm_logtpm = logtpm_corrected                # new slot added
)

dim(exprs(filtered_eset))
dim(assayData(filtered_eset)$norm_logtpm)
ls(assayData(filtered_eset))
colnames(pData(filtered_eset))
colnames(fData(filtered_eset))
table(pData(filtered_eset)$study)

library(Biobase) 

# Assuming filtered_eset and logtpm_corrected have been correctly loaded into the environment

# -----------------------------------------------
# Step A: Resolve the issue of "logical subscript is too long"
# -----------------------------------------------
corrected_sample_ids <- colnames(logtpm_corrected)
metadata_matched <- pData(filtered_eset)[corrected_sample_ids, ]

luad_samples_final <- metadata_matched$dataset == "TCGA_LUAD"
gtex_samples_final <- metadata_matched$dataset == "GTEx_LUNG"

luad_corrected_expr <- logtpm_corrected[, luad_samples_final]
gtex_corrected_expr <- logtpm_corrected[, gtex_samples_final]


# -----------------------------------------------
# Step B: Resolve the Ensembl ID version number issue (PANDA error cause)
# -----------------------------------------------

# Use the simplest sub function to remove the dot and everything after it
rownames(luad_corrected_expr) <- sub("\\..*", "", rownames(luad_corrected_expr))
rownames(gtex_corrected_expr) <- sub("\\..*", "", rownames(gtex_corrected_expr))


# -----------------------------------------------
# Step C: Format and resave the files
# -----------------------------------------------

# Format and save the LUAD file
luad_panda_df <- as.data.frame(luad_corrected_expr)
luad_panda_df$ID <- rownames(luad_panda_df)
luad_panda_df <- luad_panda_df[, c(ncol(luad_panda_df), 1:(ncol(luad_panda_df)-1))] # Move ID column to the first position

write.table(
  luad_panda_df,
  file = "luad_corrected_expression_final.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE # Ensure R's internal row names are not written
)

# Format and save the GTEx file
gtex_panda_df <- as.data.frame(gtex_corrected_expr)
gtex_panda_df$ID <- rownames(gtex_panda_df)
gtex_panda_df <- gtex_panda_df[, c(ncol(gtex_panda_df), 1:(ncol(gtex_panda_df)-1))]

write.table(
  gtex_panda_df,
  file = "gtex_corrected_expression_final.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --- Files have been resaved ---


# --- Specific fix for GTEx file writing (using the most robust data.frame method) ---

# 1. Reformat the GTEx data frame (using explicit column definition, no index reordering needed)
# Ensure the ID column is the first column and only appears once.
gtex_panda_df_fixed <- data.frame(
  ID = rownames(gtex_corrected_expr), # Explicitly define the ID column
  as.data.frame(gtex_corrected_expr) # Followed by expression values for all samples
)

# 2. Overwrite and save the GTEx file
write.table(
  gtex_panda_df_fixed,
  file = "gtex_corrected_expression_final.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)