library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

#==Set Rattus norvegicus as organism annotation to load==
library("org.Rn.eg.db", character.only = TRUE)  # Load annotation for Rattus norvegicus
keytypes(org.Rn.eg.db)  # display available datasets to use

#==Generate visualisations of DESeq2 Output==
# Principal Component Analysis
rld <- rlog(dds)
plotPCA(rld,intgroup=c("sample_status"), ntop = 500, returnData=FALSE) + 
  geom_point(aes(color=group)) + 
  scale_color_manual(values = c('#244061','#FFC107')) +
  ggtitle("PCA of Differential Gene Expression") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  colour="#3a414a",
                                  face = "bold"))

# MA plot
# sort results most to least significant (by Padj) so all NS are plotted first

ggmaplot(
  data = res[order(res$padj),],
  fdr = 0.05,
  fc = 4, # a fold change of 4 is log2 fold change 2
  genenames = NULL,
  detection_call = NULL,
  size = 1.25,
  alpha = 0.7,
  seed = 42,
  font.label = c(8, "plain", "black"),
  label.rectangle = FALSE,
  palette = c("#40B0A6", "#B32240", "#E4DEE0"),
  top = 10,
  select.top.method = c("padj", "fc"),
  label.select = NULL,
  xlab = expression(log[2]~mean~expression),
  ylab = expression(log[2]~fold~change),
  ggtheme = theme_classic()) + 
  ggtitle("Differentially Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  colour="#3a414a",
                                  face = "bold"))


#==Prep data for clusterProfiler==
# read.csv("results.csv", header = TRUE)  # read in CSV if necessary
og_genes <- res$log2FoldChange
genes <- na.omit(og_genes)
genes = sort(genes, decreasing = TRUE)

#==