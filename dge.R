# Note: featureCounts output needed to be converted to CSV from TSV prior to this R script

#===load libraries===
library(DESeq2)
library(tidyverse)
library(clusterProfiler)

#==Set Rattus norvegicus as organism annotation to load==
organism = "org.Rn.eg.db"
library(organism, character.only = TRUE)  # Load annotation for Rattus norvegicus


#===read in counts data===
featureCounts_out <- read.csv('counts_data.csv')
head(featureCounts_out) # contains extraneous columns, new minimal dataframe
count_data <- dplyr::select(featureCounts_out, "Geneid", "C1", "C2", "C3", "I1", "I2", "I3")
rownames(count_data) <- count_data[,1]    # make geneid the rownames
count_data[,1] <- NULL

#===assigning samples status===
col_data <- data.frame(
  row.names = colnames(count_data),
  sample_status = c("control", "control", "control", "infected", "infected", "infected")
)

all(colnames(count_data) %in% rownames(col_data)) # checks
all(colnames(count_data) == rownames(col_data))

#===construct DESeqDataSet object===
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ sample_status)

dds
pre_filter_len <- nrow(dds)

#===filtering===
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] # removing rows with under 10 read counts, total

post_filter_len <- nrow(dds)
rows_rmvd <- pre_filter_len - post_filter_len
rows_rmvd

#===set control status as norm to be compared to===
dds$sample_status <- relevel(dds$sample_status, ref = "control")

#===Run DESeq===
dds <- DESeq(dds)
res <- results(dds)

res

write.csv(res, "result.csv")

#==Gene Set Enrichment==

# Prep data for clusterProfiler
og_genes <- res$log2FoldChange
names(og_genes) <- rownames(res)
genes <- na.omit(og_genes)
genes = sort(genes, decreasing = TRUE)

length(genes)

# Gene Set Enrichment
keytypes(org.Rn.eg.db)  # display available datasets to use
gse <- gseGO(geneList=genes, 
             ont ="BP",  #Biological Process, can also take "CC", "MS" or "ALL"
             keyType = "SYMBOL", 
             minGSSize = 3, # allowing group sizes of 3-800
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr")

