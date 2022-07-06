if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install("TCGAbiolinks")
BiocManager::install("dplyr")
BiocManager::install("SummarizedExperiment")
BiocManager::install("edgeR")
BiocManager::install("tidyverse")
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)
TCGAbiolinks:::getProjectSummary("TCGA-SKCM")
########RNASeq SKCM
query_RNAseq = GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification"
)
query_infoRNASeq = getResults(query_RNAseq) #query info is compatible with the TCGA portal
colnames(query_infoRNASeq)  ##for clinical data---> review code block.

GDCdownload(
  query = query_RNAseq,
  method = "client",
  directory = "TCGA_TranscriptomeProfiling_SKCM"
)
TCGA_SKCM_Transcriptome = GDCprepare(
  query = query_RNAseq,
  save = TRUE,
  save.filename = "prepdataSKCMRNAseq",
  directory = "TCGA_TranscriptomeProfiling_SKCM"
)
RNASeq_genes_SKCM = rowData(TCGA_SKCM_Transcriptome)
SKCM_RNASeqraw = assay(TCGA_SKCM_Transcriptome)

######## miRNASeq SKCM
query_miRSeq = GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling",
  data.type = "miRNA Expression Quantification"
  
)
query_infomiRSeq = getResults(query_miRSeq) #query data is compatible with the TCGA portal.
colnames(query_miRSeq)
GDCdownload(
  query = query_miRSeq,
  method = "client",
  directory = "miRNASeqraw.SKCM"
)
TCGA_SKCM_miRSeq = GDCprepare(
  query = query_miRSeq,
  save = TRUE,
  save.filename = "prepdatamiRSeq",
  directory = "miRNASeqraw.SKCM"
)
#RNASeq Normalization & gene name append
which(colSums(is.na(SKCM_RNASeqraw)) > 0)    #named integer(0)

RNAseq_normalized = cpm(SKCM_RNASeqraw, log = FALSE, prior.count = 1)

gene_ID_RNASeq = c(RNASeq_genes_SKCM[,"gene_name"])
RNAseq_normalized_wgenename = RNAseq_normalized
rownames(RNAseq_normalized_wgenename) = gene_ID_RNASeq

#miRSeq process-normalization
SKCM_miRSeqraw = select(TCGA_SKCM_miRSeq, c(starts_with("miRNA_ID"),starts_with("read_count") ))

miRNA_names = select(SKCM_miRSeqraw, c(starts_with("miRNA_ID")))

miRNA_counts = select(SKCM_miRSeqraw, starts_with("read_count"))

miRNA_counts_normalized = cpm(miRNA_counts, log = FALSE, prior.count = 1)

miRSeq_normalized = cbind(miRNA_names,miRNA_counts_normalized)

colnames(miRSeq_normalized) = gsub("read_count_","",colnames(miRSeq_normalized))###Look this later. 
colnames(miRNA_counts_normalized) = gsub("read_count_","",colnames(miRNA_counts_normalized))


#####CLinical Data
colnames(colData(TCGA_SKCM_Transcriptome))

skcm_clinical_whole_data = colData(TCGA_SKCM_Transcriptome)
