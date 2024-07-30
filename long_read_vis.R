library("ggseqlogo")

#==load in site level m6Anet results==
m6a_site <- read.csv('data.site_proba_attempt1_b.csv')
head(m6a_site)

# replace w/ following once all results calc'ed
#c2_m6A_site <- read.csv('data.site_proba_c2.csv')
#c3_m6A_site <- read.csv('data.site_proba_c3.csv')
#i2_m6A_site <- read.csv('data.site_proba_i2.csv')
#i3_m6A_site <- read.csv('data.site_proba_i3.csv')

motifVector <- as.vector(m6a_site$kmer)
head(motifVector)

# replace w/ following when need to concatenate controls vectors
#motifs_control <- as.vector(c(c2_m6a_site$kmer, c3_m6a_site$kmer)
#motifs_infected <- as.vector(c(i2_m6a_site$kmer, i3_m6a_site$kmer)

#==convert Us to Ts for true motifs==
motifVector <- gsub(pattern = "T", replacement = "U", x = motifVector)
head(motifVector)

# replace w/ following when multiple motif vectors exist
#motifs_c <- gsub(pattern = "T", replacement = "U", x = motifs_c)
#motifs_i <- gsub(pattern = "T", replacement = "U", x = motifs_i)

#==create custom colour scheme for colourblind accessibility==
colours = make_col_scheme(chars=c('A', 'U', 'C', 'G'),
                      cols=c('#40B0A6', '#B32240', '#244061', '#FFC107'))
#B32213 - other 'red' option not used

#==generate sequence logo(s and align plots)==
ggseqlogo(motifVector, method = 'prob', col_scheme = colours)

# replace w/ following when multiple motif vectors exist
#control_logo <- ggseqlogo(motifs_c, method = 'prob', col_scheme = colours)
#infected_logo <- ggseqlogo(motifs_i,method = 'prob', col_scheme = colours)
#suppressMessages( require(cowplot) )
#plot_grid(control_logo,infected_logo,  ncol = 1, align = 'v')
