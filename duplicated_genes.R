# library(TCGAbiolinks)
# library(dplyr)
# library(SummarizedExperiment)
# library(edgeR)
# library(tidyverse)
#To identify duplicated rows (genes) in RNASeq_normalized_wgenenames data and create another matrix with them. 

sum(duplicated(RNAseq_normalized_wgenename))
duplicate_logicals = duplicated(RNAseq_normalized_wgenename)

df_duplicated_genes = data.frame()
i = 1
while(i <= length(duplicate_logicals)) {
  if (duplicate_logicals[i] == TRUE) {
    df_duplicated_genes=rbind(df_duplicated_genes,df_RNASeq_normalized_wgenename[i,])
  }
  i = i + 1
}
sum(df_duplicated_genes["AL358613.3",])
sum(RNAseq_normalized_wgenename["AL358613.3",])

x = df_duplicated_genes[rowSums(df_duplicated_genes[]) > 0, ]
y = df_duplicated_genes[rowSums(df_duplicated_genes[]) ==  0, ]












