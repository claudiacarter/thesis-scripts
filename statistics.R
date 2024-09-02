# check significance of GDMT enrichment in DEGs with hypergeometric test
deg_group = 672
gdmt_group = 665
group_overlap = 16
rn_genes = 23194
phyper(group_overlap,deg_group, rn_genes-deg_group, gdmt_group, lower.tail=FALSE)