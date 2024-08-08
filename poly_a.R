library(readr)
library(dplyr)

# read in nanopolish output CSVs
np_out_c2 <- as.data.frame(read_tsv("c2_polya_out.tsv"))
np_pass_c2 <- filter(np_out_c2, qc_tag == "PASS")

np_out_c3 <- as.data.frame(read_tsv("c3_polya_out.tsv"))
np_pass_c3 <- filter(np_out_c3, qc_tag == "PASS")

np_out_i2 <- as.data.frame(read_tsv("i2_polya_out.tsv"))
np_pass_i2 <- filter(np_out_i2, qc_tag == "PASS")

np_out_i3 <- as.data.frame(read_tsv("i3_polya_out.tsv"))
np_pass_i3 <- filter(np_out_i3, qc_tag == "PASS")
