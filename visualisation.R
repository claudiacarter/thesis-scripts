# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# BiocManager::install(organism, character.only = TRUE)
# devtools::install_github('btmonier/ggDESeq')

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
# options for data vis:
library(ggDESeq)
library(ggpubr)

#==Set Rattus norvegicus as organism annotation to load==
library("org.Rn.eg.db", character.only = TRUE)  # Load annotation for Rattus norvegicus
keytypes(org.Rn.eg.db)  # display available datasets to use

#==Generate visualisations of DESeq2 Output==
# using ggDESeq:
#ggMA(data=dds,padj=0.01,lfc=1)
ggVolcano(data = dds, padj = 0.1)
ggMA(data=dds, padj=0.05, lfc=2)

# using DESeq's own plots:
plotMA(dds, ylim = c(-2, 2))

# using ggpubr:
ggmaplot(
  data = res,
  fdr = 0.05,
  fc = 4,
  genenames = NULL,
  detection_call = NULL,
  size = 1,
  alpha = 0.5,
  seed = 42,
  font.label = c(8, "plain", "black"),
  label.rectangle = FALSE,
  palette = c("#4195F5", "#D21359", "#E4DEE0"),
  top = 10,
  select.top.method = c("padj", "fc"),
  label.select = NULL,
  main = "Differential Gene Expression MA Plot",
  xlab = expression(log[2]~mean~expression),
  ylab = expression(log[2]~fold~change),
  ggtheme = theme_classic())


#==Prep data for clusterProfiler==
# read.csv("results.csv", header = TRUE)  # read in CSV if necessary
og_genes <- res$log2FoldChange
genes <- na.omit(og_genes)
genes = sort(genes, decreasing = TRUE)

#==