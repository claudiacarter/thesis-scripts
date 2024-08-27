# define control_percentmeth and infected_percent_meth
#read csv of each FULL m6anet result

# t tests
sink(file="welch_test.txt")
t.test(control_percentmeth, infected_percent_meth, var.equal = FALSE)
sink(file=NULL)

sink(file="t_test.txt")
t.test(control_percentmeth, infected_percent_meth, var.equal = TRUE)
sink(file=NULL)

# check significance of GDMT enrichment in DEGs with hypergeometric test
deg_group = 672
gdmt_group = 665
group_overlap = 16
rn_genes = 23194
phyper(group_overlap,deg_group, rn_genes-deg_group, gdmt_group, lower.tail=FALSE)