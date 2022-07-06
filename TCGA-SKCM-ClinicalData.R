library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)

#####CLinical Data
colnames(colData(TCGA_SKCM_Transcriptome))

skcm_clinical_whole_data = colData(TCGA_SKCM_Transcriptome) ##will be analyzed and processed. 