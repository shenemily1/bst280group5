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
eset2_luad <- eset_luad
eset2_gtex_lung <- eset_gtex_lung
expr_luad     <- exprs(eset_luad)
expr_gtex     <- exprs(eset_gtex_lung)
counts_luad   <- assayData(eset_luad)$counts
counts_gtex   <- assayData(eset_gtex_lung)$counts
logtpm_luad   <- assayData(eset_luad)$log_tpm
logtpm_gtex   <- assayData(eset_gtex_lung)$log_tpm
pdata_luad    <- pData(eset_luad)
pdata_gtex    <- pData(eset_gtex_lung)
fdata_luad    <- fData(eset_luad)
fdata_gtex    <- fData(eset_gtex_lung)

#Find intersect fData and manipulate individual eset accordingly
common_f_cols <- intersect(colnames(fdata_luad), colnames(fdata_gtex))
unique_luad   <- setdiff(colnames(fdata_luad), colnames(fdata_gtex))
unique_gtex   <- setdiff(colnames(fdata_gtex), colnames(fdata_luad))
for (col in unique_luad) {fdata_gtex[[col]] <- NA}
for (col in unique_gtex) {fdata_luad[[col]] <- NA}
fdata_luad <- fdata_luad[, sort(colnames(fdata_luad))]
fdata_gtex <- fdata_gtex[, sort(colnames(fdata_gtex))]
fData(eset2_luad)      <- fdata_luad
fData(eset2_gtex_lung) <- fdata_gtex

#Find intersect pData and manipulate individual eset accordingly
common_p_cols <- intersect(colnames(pdata_luad), colnames(pdata_gtex))
unique_luad_p <- setdiff(colnames(pdata_luad), colnames(pdata_gtex))
unique_gtex_p <- setdiff(colnames(pdata_gtex), colnames(pdata_luad))
for (col in unique_luad_p) {pdata_gtex[[col]] <- NA}
for (col in unique_gtex_p) {pdata_luad[[col]] <- NA}
pdata_luad <- pdata_luad[, sort(colnames(pdata_luad))]
pdata_gtex <- pdata_gtex[, sort(colnames(pdata_gtex))]
pData(eset2_luad)      <- pdata_luad
pData(eset2_gtex_lung) <- pdata_gtex

#find interecting expressed gene and subset the eset files 
common_genes <- intersect(rownames(eset_luad), rownames(eset_gtex_lung))
eset3_luad      <- eset2_luad[common_genes, ]
eset3_gtex_lung <- eset2_gtex_lung[common_genes, ]

#combine the assayData, pdata, fdata into a new combined dataset
expr_combined <- cbind(exprs(eset3_luad), exprs(eset3_gtex_lung))
counts_combined <- cbind(assayData(eset3_luad)$counts, assayData(eset3_gtex_lung)$counts)
logtpm_combined <- cbind(assayData(eset3_luad)$log_tpm,assayData(eset3_gtex_lung)$log_tpm)
pData(eset2_luad)$dataset      <- "tcga_luad"
pData(eset2_gtex_lung)$dataset <- "gtex_lung"
pdata_combined <- rbind(pData(eset2_luad),pData(eset2_gtex_lung))
fdata_combined <- fData(eset2_luad) #fData should be the same for luad and lung now 


combined_eset <- ExpressionSet(
  assayData = assayDataNew(
    exprs   = expr_combined,
    counts  = counts_combined,
    log_tpm = logtpm_combined
  ),
  phenoData   = AnnotatedDataFrame(pdata_combined),
  featureData = AnnotatedDataFrame(fdata_combined)
)
validObject(combined_eset)

#Create counts per million assayData element 'tpm' and filter accordingly
#Filter criteria 1: average row mean for tpm greater than 1 is > 0.2
#Filter criteria 2: average row mean for raw counts greater than 10 is > 50%
counts_luad <- assayData(eset_luad)$counts
counts_gtex <- assayData(eset_gtex_lung)$counts
tpm_luad <- 2^(assayData(eset_luad)$log_tpm) - 1
tpm_gtex <- 2^(assayData(eset_gtex_lung)$log_tpm) - 1
counts_all <- cbind(counts_luad, counts_gtex)
tpm_all    <- cbind(tpm_luad, tpm_gtex)
tpm_filter <- rowMeans(tpm_all > 1) >= 0.20
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
table(pData(filtered_eset)$study)