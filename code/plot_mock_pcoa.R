library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]

# file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes"

read_tsv(file_name) %>%
	select(group, axis1, axis2) %>%
	separate(group, into=c("rounds", "polymerase", "mock")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=axis1, y=axis2, color=rounds, shape=polymerase)) +
		geom_point() +
		stat_ellipse(aes(group=rounds)) +
		ggsave("results/figures/mock_pcoa.pdf")
