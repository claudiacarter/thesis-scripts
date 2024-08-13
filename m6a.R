
# read in differential gene expression (DGE) data
all_genes <- as.data.frame(res)
all_genes <- all_genes[c("log2FoldChange", "padj")]
colnames(all_genes) <- c("DGE_log2FC", "DGE_Padj")

# read in site-level m6A data
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

rownames(sites_per_transcript)

# Convert IDs to Genes?
write.table(rownames(sites_per_transcript), 
            append = FALSE, 
            "m6A_genes.txt",
            sep = "\n", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)