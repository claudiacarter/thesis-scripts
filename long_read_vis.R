library("ggseqlogo")

#==load in site level m6Anet results==
c2_m6A_site <- read.csv('c2_data.site_proba_tra_filtered.csv')
c3_m6A_site <- read.csv('c3_data.site_proba_tra_filtered.csv')
i2_m6A_site <- read.csv('i2_data.site_proba_tra_filtered.csv')
i3_m6A_site <- read.csv('i3_data.site_proba_tra_filtered.csv')

motifs_c <- as.vector(c(c2_m6A_site$kmer, c3_m6A_site$kmer))
motifs_i <- as.vector(c(i2_m6A_site$kmer, i3_m6A_site$kmer))

#==convert Us to Ts for true motifs==
motifs_c <- gsub(pattern = "T", replacement = "U", x = motifs_c)
motifs_i <- gsub(pattern = "T", replacement = "U", x = motifs_i)


#==create custom colour scheme for colourblind accessibility==
colours = make_col_scheme(chars=c('A', 'U', 'C', 'G'),
                      cols=c('#40B0A6', '#B32240', '#244061', '#FFC107'))
#B32213 - other 'red' option not used

#==generate sequence logos and align plots==
# plots per sample


# plots control group vs infected group
control_logo <- ggseqlogo(motifs_c, method = 'prob', col_scheme = colours) +
               ggtitle("a. m6A DRACH Motifs of Control Group") +
               theme(plot.title = element_text(hjust = 0.5))
infected_logo <- ggseqlogo(motifs_i,method = 'prob', col_scheme = colours) +
               ggtitle("b. m6A DRACH Motifs of Infected Group") +
               theme(plot.title = element_text(hjust = 0.5))
plot_grid(control_logo,infected_logo,  ncol = 1, align = 'v')
