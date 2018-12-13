library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]

# file_name <- "data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes"

read_tsv(file_name) %>%
	select(group, axis1, axis2) %>%
	mutate(group = str_replace(group, "_PMM_FS", "")) %>%
	separate(group, into=c("rounds", "polymerase", "subject")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=axis1, y=axis2, color=polymerase, shape=rounds)) +
		stat_ellipse(aes(group=subject), type="norm", geom="polygon", color="gray", fill="lightgray") +
		geom_point() +
		ggsave("results/figures/stool_pcoa.pdf")
