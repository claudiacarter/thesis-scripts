library("ggseqlogo")

#==visualising the DRACH motifs methylated==
motifs_c2 <- as.vector(c2_m6A_site$kmer)
motifs_c3 <- as.vector(c3_m6A_site$kmer)
motifs_i2 <- as.vector(i2_m6A_site$kmer)
motifs_i3 <- as.vector(i3_m6A_site$kmer)

motifs_c <- as.vector(c(motifs_c2, motifs_c3))
motifs_i <- as.vector(c(motifs_i2, motifs_i3))

# convert U to T
motifs_c2 <- gsub(pattern = "T", replacement = "U", x = motifs_c2)
motifs_c3 <- gsub(pattern = "T", replacement = "U", x = motifs_c3)
motifs_i2 <- gsub(pattern = "T", replacement = "U", x = motifs_i2)
motifs_i3 <- gsub(pattern = "T", replacement = "U", x = motifs_i3)

motifs_c <- gsub(pattern = "T", replacement = "U", x = motifs_c)
motifs_i <- gsub(pattern = "T", replacement = "U", x = motifs_i)


# create custom colour scheme for colourblind accessibility
colours = make_col_scheme(chars=c('A', 'U', 'C', 'G'),
                      cols=c('#40B0A6', '#B32240', '#244061', '#FFC107'))
#B32213 - other 'red' option not used

# generate sequence logos and align plots
# plots per sample
c2_logo <- ggseqlogo(motifs_c2,
                     method = 'prob',
                     col_scheme = colours) + 
           ggtitle("m6A DRACH Motifs of Control 2") + 
           theme(plot.title = element_text(hjust = 0.5))
c3_logo <- ggseqlogo(motifs_c3,
                     method = 'prob',
                     col_scheme = colours) + 
           ggtitle("m6A DRACH Motifs of Control 3") + 
           theme(plot.title = element_text(hjust = 0.5))
i2_logo <- ggseqlogo(motifs_i2,
                     method = 'prob',
                     col_scheme = colours) + 
           ggtitle("m6A DRACH Motifs of Infected 2") + 
           theme(plot.title = element_text(hjust = 0.5))
i3_logo <- ggseqlogo(motifs_i3,
                     method = 'prob',
                     col_scheme = colours) + 
           ggtitle("m6A DRACH Motifs of Infected 3") + 
           theme(plot.title = element_text(hjust = 0.5))

c2_logo
c3_logo
i2_logo
i3_logo
plot_grid(c2_logo, c3_logo, i2_logo, i3_logo, ncol = 2, align = 'v')

# plots control group vs infected group
control_logo <- ggseqlogo(motifs_c, 
                          method = 'prob', 
                          col_scheme = colours) +
                ggtitle("m6A DRACH Motifs of Control Group") +
                theme(plot.title = element_text(hjust = 0.5, 
                                                colour="#3a414a",
                                                face = "bold"))
                                                
infected_logo <- ggseqlogo(motifs_i,
                           method = 'prob', 
                           col_scheme = colours) +
                 ggtitle("m6A DRACH Motifs of Infected Group") + 
                 theme(plot.title = element_text(hjust = 0.5, 
                                                 colour="#3a414a",
                                                 face = "bold"))

plot_grid(control_logo, infected_logo, ncol = 1, align = 'v')

#==compare transcript methylation levels==

dge <- filter(degs, gene_symbol %in% datapoints)
dm <- filter(filtered_per_transcript_df, gene_symbol %in% datapoints)
dge_dm <- merge(dge, dm)
dge_dm

dge_dm_plot <- ggplot(dge_dm, 
                      aes(DGE_log2FC, DM_log2FC, label = label),
                      scale="globalminmax") +
                     geom_vline(xintercept = 0, linetype = 1, colour = "#3a414a") +
                     geom_hline(yintercept = 0, linetype = 1, colour = "#3a414a") +
                     geom_point(colour = ifelse(dge_dm$DM_log2FC != 0, "#244061", "#707b7c")) +
                     geom_text_repel(max.overlaps = 5, size = 3) + 
                     theme_grey() +
                     ggtitle("Differential Methylation vs Differential Expression") +
                     theme(plot.title = element_text(hjust = 0.5, 
                                      colour="#3a414a",
                                      face = "bold"))