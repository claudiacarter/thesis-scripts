library(tidyverse)

#==load in site level m6Anet results==
c2_m6A_site <- read.csv('c2_data.site_proba_tra_filtered.csv')
c3_m6A_site <- read.csv('c3_data.site_proba_tra_filtered.csv')
i2_m6A_site <- read.csv('i2_data.site_proba_tra_filtered.csv')
i3_m6A_site <- read.csv('i3_data.site_proba_tra_filtered.csv')

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

# translate refseq IDs into gene symbols
refseq_ids <- as.list(rownames(per_transcript_df))
write.table(refseq_ids, 
            append = FALSE, 
            "m6A_genes.txt",
            sep = "\n", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

#~~~submitted txt file to web tool NCBI databases for gene names and symbols~~~#

gene_table <- read.csv("transcripts_to_genes.csv")
head(gene_table)

# checking all transcripts are present and in the same order as methylation data
all(rownames(per_transcript_df) %in% gene_table$RefSeq.Transcript.ID)
all(rownames(per_transcript_df) == gene_table$RefSeq.Transcript.ID)

# add gene symbol column into dataframe
per_transcript_df <- cbind(gene_symbol=gene_table$Symbol, per_transcript_df)
head(per_transcript_df)



# Subset genes that didn't have 0 as mean of control or infected

# This is to make sure I'm only dealing with differentially methylated genes 
# and not just ones that aren't expressed in control/infectious state.
# I'm sure there's a better way to do this but for the purposes of the plot 
# it's okay.

filtered_per_transcript_df <- subset(per_transcript_df,
                                      DM_log2FC != Inf & DM_log2FC != -Inf,
                                      select=c(gene_symbol, DM_log2FC))

length(rownames(filtered_per_transcript_df))
head(filtered_per_transcript_df)

# find duplicates (i.e. multiple transcripts -isoforms?- for one gene)
transcript_freq <- data.frame(table(filtered_per_transcript_df$gene_symbol))
transcript_freq <- subset(transcript_freq,
                          Freq != 1)


head(transcript_freq)

