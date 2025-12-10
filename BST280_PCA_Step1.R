library(ggfortify)
library(limma)
# PCA before batch correction
PCA_filter_prebatch <- prcomp(t(assayData(filtered_eset)[["log_tpm"]]))
# PCA after batch correction
PCA_filter_postbatch <- prcomp(t(assayData(filtered_eset)[["norm_logtpm"]]))


#autoplot(PCA_filter_prebatch, data=pData(filtered_eset), colour="study",
#         main="Before batch correction")
study <- pData(filtered_eset)$dataset 

#plotting PCA pre and post batch correction for tcga / gtex dataset with base r
plot(
  PCA_filter_prebatch$x[,1], PCA_filter_prebatch$x[,2],
  col = as.factor(study),
  pch = 19,
  xlab = paste0("PC1 (", round(summary(PCA_filter_prebatch)$importance[2,1] * 100, 1), "%)"),
  ylab = paste0("PC2 (", round(summary(PCA_filter_prebatch)$importance[2,2] * 100, 1), "%)"),
  main = "PCA Before Batch Correction"
)
legend("topright", legend = levels(as.factor(study)),
       col = 1:length(levels(as.factor(study))), pch = 19)

plot(
  PCA_filter_postbatch$x[,1], PCA_filter_postbatch$x[,2],
  col = as.factor(study),
  pch = 19,
  xlab = paste0("PC1 (", round(summary(PCA_filter_postbatch)$importance[2,1] * 100, 1), "%)"),
  ylab = paste0("PC2 (", round(summary(PCA_filter_postbatch)$importance[2,2] * 100, 1), "%)"),
  main = "PCA After Batch Correction"
)
legend("topright", legend = levels(as.factor(study)),
       col = 1:length(levels(as.factor(study))), pch = 19)

saveRDS(filtered_eset, file = "luad_gtex_eset_1207.rds")
##
filtered_eset <- readRDS("luad_gtex_eset_1207.rds")
study <- pData(filtered_eset)$dataset 
PCA_filter_postbatch <- prcomp(t(assayData(filtered_eset)[["norm_logtpm"]]))

pc1 <- PCA_filter_postbatch$x[,1]
plot(
  density(pc1[study == "GTEX_LUNG"]),
  col = "blue", lwd = 3,
  main = "Distribution of PC1: GTEx (homogeneous) vs LUAD (heterogeneous)",
  xlab = "PC1"
)
lines(
  density(pc1[study == "TCGA_LUAD"]),
  col = "red", lwd = 3
)
legend("topleft",
       legend = c("GTEx Normal Lung", "TCGA-LUAD Tumor"),
       col = c("blue", "red"),
       lwd = 3)

#To do list, derive common gender/sex covariate, race/ethnicity covariate, smoking status covariate
#histological subtype covariate for the clustering analysis from the tcga and gtex unique variables
#Visualize varaince difference along the PC1 axis