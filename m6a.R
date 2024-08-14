if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager", version='3.18')
 
BiocManager::install(version='3.18')
BiocManager::install("biomaRt")

library(biomaRt)
library(tidyverse)

# read in differential gene expression (DGE) data
all_genes <- as.data.frame(res)
all_genes <- all_genes[c("log2FoldChange", "padj")]
colnames(all_genes) <- c("DGE_log2FC", "DGE_Padj")

# table freq.s of methylation sites per represented transcript (non-redundant)
# note: data pre-filtered for >0.9 probability

transcripts <- sort(unique(c(c2_m6A_site$transcript_id, 
                             c3_m6A_site$transcript_id, 
                             i2_m6A_site$transcript_id, 
                             i3_m6A_site$transcript_id)))

c2 <- factor(c2_m6A_site$transcript_id, transcripts)
c3 <- factor(c3_m6A_site$transcript_id, transcripts)
i2 <- factor(i2_m6A_site$transcript_id, transcripts)
i3 <- factor(i3_m6A_site$transcript_id, transcripts)


sites_per_transcript <- cbind(freq_c2=table(c2), 
                              freq_c3=table(c3), 
                              freq_i2=table(i2), 
                              freq_i3=table(i3))

head(sites_per_transcript)
per_transcript_df <- as.data.frame(sites_per_transcript, header=TRUE)

# Add control average (baseMean) and infected average and log2FC columns

per_transcript_df$DM_baseMean <- (per_transcript_df$freq_c2 + per_transcript_df$freq_c3) / 2
per_transcript_df$infectedMean <- (per_transcript_df$freq_i2 + per_transcript_df$freq_i3) / 2
per_transcript_df$DM_log2FC <- log2(per_transcript_df$infectedMean/per_transcript_df$DM_baseMean)
head(per_transcript_df)

# translate refseq IDs into gene symbols (submitted to NCBI site)
refseq_ids <- as.list(rownames(per_transcript_df))

write.table(refseq_ids, 
            append = FALSE, 
            "m6A_genes.txt",
            sep = "\n", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)



# merge gene symbol column into dataframe
cbind(gene_symbol=0, per_transcript_df)

# still need to differentiate between transcripts not expressed and transcripts
# with no methylation?
