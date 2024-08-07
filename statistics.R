# define control_percentmeth and infected_percent_meth


# t tests
sink(file="welch_test.txt")
t.test(control_percentmeth, infected_percent_meth, var.equal = FALSE)
sink(file=NULL)

sink(file="t_test.txt")
t.test(control_percentmeth, infected_percent_meth, var.equal = TRUE)
sink(file=NULL)
