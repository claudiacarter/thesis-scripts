library("ggseqlogo")

#==load in site level m6Anet results==
c2_m6A_site <- read.csv('c2_data.site_proba_tra_filtered.csv')
c3_m6A_site <- read.csv('c3_data.site_proba_tra_filtered.csv')
i2_m6A_site <- read.csv('i2_data.site_proba_tra_filtered.csv')
i3_m6A_site <- read.csv('i3_data.site_proba_tra_filtered.csv')

motifs_c2 <- as.vector(c2_m6A_site$kmer)
motifs_c3 <- as.vector(c3_m6A_site$kmer)
motifs_i2 <- as.vector(i2_m6A_site$kmer)
motifs_i3 <- as.vector(i3_m6A_site$kmer)

motifs_c <- as.vector(c(motifs_c2, motifs_c3))
motifs_i <- as.vector(c(motifs_i2, motifs_i3))

#==convert Us to Ts for true motifs==
motifs_c2 <- gsub(pattern = "T", replacement = "U", x = motifs_c2)
motifs_c3 <- gsub(pattern = "T", replacement = "U", x = motifs_c3)
motifs_i2 <- gsub(pattern = "T", replacement = "U", x = motifs_i2)
motifs_i3 <- gsub(pattern = "T", replacement = "U", x = motifs_i3)

motifs_c <- gsub(pattern = "T", replacement = "U", x = motifs_c)
motifs_i <- gsub(pattern = "T", replacement = "U", x = motifs_i)


#==create custom colour scheme for colourblind accessibility==
colours = make_col_scheme(chars=c('A', 'U', 'C', 'G'),
                      cols=c('#40B0A6', '#B32240', '#244061', '#FFC107'))
#B32213 - other 'red' option not used

#==generate sequence logos and align plots==
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
