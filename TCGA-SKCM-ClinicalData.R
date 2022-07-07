library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)

BiocManager::install("S4Vectors")
library(S4Vectors)
#####CLinical Data
colnames(colData(TCGA_SKCM_Transcriptome))

skcm_clinical_whole_data = colData(TCGA_SKCM_Transcriptome)  

# for (p in 1:58) {
#   if (isS4(skcm_clinical_whole_data[p]) == TRUE) {
#     print(p)
#   }
#   
# }
# 
# for (p in 1:58) {
#   if (nrow(skcm_clinical_whole_data[p]) == 473) {
#     print(p)
#   } else if (nrow(skcm_clinical_whole_data < 473)) {
#     print("not", p)
#   }
#   
# }
skcm_clinical = data.frame(skcm_clinical_whole_data[1:58])
colsremoved = c("barcode","patient","sample","sample_submitter_id", "sample_id", "submitter_id", "exposure_id", "demographic_id")
skcm_clinical = select(skcm_clinical, -colsremoved)





