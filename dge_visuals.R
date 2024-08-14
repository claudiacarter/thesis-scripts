library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

#==Set Rattus norvegicus as organism annotation to load== # check still needed
organism = "org.Rn.eg.db"
library(organism, character.only = TRUE)  # Load annotation for Rattus norvegicus
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
ma_data <- na.omit(res[order(res$padj),])

ggmaplot(
  data = ma_data,
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

#==Gene Set Enrichment Plots==
# DotPlot [not used in final report]
require(DOSE) # may need to cite this in plot legend/ somewhere
dotplot(gse, showCategory=12, split=".sign") + 
  facet_grid(.~.sign) + 
  scale_fill_continuous(name = "FDR", low='#40B0A6', high="#244061") +
  theme(text = element_text(size=10)) + 
  theme(axis.title.x=element_text(size=10)) +
  theme(axis.text.x=element_text(size=rel(0.9))) + 
  theme(axis.text.y=element_text(size=rel(0.9))) + 
  ggtitle("Enrichment of Biological Process GO Terms in DEGs") +
  theme(plot.title = element_text(size=10, 
                                  hjust = 0.5, 
                                  colour="#3a414a",
                                  face = "bold"))

# Ridge Plot [saving dimensions = w785*h900 for report]
ridgeplot(gse) + labs(x = "Enrichment Distribution") +
  scale_fill_continuous(name = "FDR", low='#40B0A6', high="#244061") +
  geom_vline(aes(xintercept=0), colour="#808080", linetype="dashed") +
  theme(text = element_text(size=10)) + 
  theme(axis.title.x=element_text(size=10)) + 
  theme(axis.text.x=element_text(size=rel(0.9))) + 
  theme(axis.text.y=element_text(size=rel(0.9))) + 
  ggtitle("Enrichment of Biological Process GO Terms in DEGs") +
            theme(plot.title = element_text(size=10, 
                                            hjust = 0.5, 
                                            colour="#3a414a",
                                            face = "bold"))
