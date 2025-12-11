#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("rnaseqGene", version = "3.8")
#BiocManager::install("DESeq2")
library("DESeq2") # for RNA-seq analysis
library("airway") # for loading the data set
library("dplyr") # for data processing
library("magrittr")
library("apeglm") # Used by DESeq2 to fit a linear model
library("ggplot2") # for visualization
library("vsn") # for variance stabilization
library("hexbin") # for visualization

#Load
filtered_eset <- readRDS("luad_gtex_eset_1207.rds")
#investigate
filtered_eset

# Extract counts and sample metadata
countdata <- assayData(filtered_eset)$counts   # genes x samples
coldata   <- pData(filtered_eset)             # samples x variables

# Create a group variable: normal (GTEx) vs tumor (TCGA)
coldata$group <- ifelse(coldata$dataset == "TCGA_LUAD", "tumor", "normal")
coldata$group <- factor(coldata$group, levels = c("normal", "tumor"))
table(coldata$group)

dds <- DESeqDataSetFromMatrix(
  countData = round(countdata),  # integers
  colData   = coldata,
  design    = ~ group
)

dds
nrow(dds)  # number of genes

#colnames(colData(dds))

#handle covariates - sex
coldata$sex <- ifelse(
  !is.na(coldata$gtex.sex), 
  coldata$gtex.sex, 
  coldata$tcga.cgc_case_gender
)

coldata <- as.data.frame(coldata)
coldata$sex <- as.character(coldata$sex)
# Recode GTEx: 1 -> MALE, 2 -> FEMALE
is_gtex <- coldata$dataset == "GTEx_LUNG"
coldata$sex[is_gtex & coldata$sex %in% c("1", 1)] <- "MALE"
coldata$sex[is_gtex & coldata$sex %in% c("2", 2)] <- "FEMALE"

coldata$sex <- factor(coldata$sex, levels = c("FEMALE", "MALE"))
table(coldata$sex, coldata$dataset)

#handle covariates - age
#colnames(pData(filtered_eset)[grep("age", colnames(pData(filtered_eset)), ignore.case = TRUE)])
#table(coldata$gtex.age)
#summary(coldata$tcga.cgc_case_age_at_diagnosis)

# TCGA ages (numeric)
tcga_age <- coldata$tcga.cgc_case_age_at_diagnosis

# Create comparable age categories
tcga_age_bins <- cut(
  tcga_age,
  breaks = c(20, 30, 40, 50, 60, 70, 80, 90),
  labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89"),
  right = FALSE
)

# Assign only to TCGA samples
tcga_age_bins[ coldata$dataset != "TCGA_LUAD" ] <- NA

coldata$age_bin <- ifelse(
  coldata$dataset == "GTEx_LUNG",
  coldata$gtex.age,
  as.character(tcga_age_bins)
)

coldata$age_bin <- factor(coldata$age_bin,
                          levels = c("20-29","30-39","40-49","50-59","60-69","70-79", "80-89"))
table(coldata$age_bin, coldata$dataset)
length(coldata$age_bin)


#no variables first
design(dds) <- ~ group
dds

dds <- DESeq(dds)

# Contrast: tumor vs normal (group_tumor / group_normal)
res <- results(dds, contrast = c("group", "tumor", "normal"))
head(res)
summary(res)

# Add gene symbols if they’re in fData
gene_annot <- fData(filtered_eset)
res_df <- as.data.frame(res)
res_df$gene_id    <- rownames(res_df)
res_df$gene_symbol <- gene_annot[rownames(res_df), "gene_symbol"]  # adjust col name as needed

head(res_df[order(res_df$padj), ], 10)  # top 10 DE genes

# Library sizes
round(colSums(countdata) / 1e6, 1)

# PCA on rlog or vst-transformed counts
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "group")

