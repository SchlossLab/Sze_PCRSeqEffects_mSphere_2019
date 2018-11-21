library(tidyverse)

pre_vsearch_file <- "data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.count_table"

pre_vsearch <- read_tsv(pre_vsearch_file) %>%
	select(-total) %>%
	gather(sample, total, -Representative_Sequence) %>%
	group_by(sample) %>%
	summarize(pre_vsearch = sum(total))

post_vsearch_file <- "data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table"

post_vsearch <- read_tsv(post_vsearch_file) %>%
	select(-total) %>%
	gather(sample, total, -Representative_Sequence) %>%
	group_by(sample) %>%
	summarize(post_vsearch = sum(total))

inner_join(pre_vsearch, post_vsearch, by="sample") %>%
	mutate(frac_chimera = 1 - post_vsearch / pre_vsearch) %>%
	separate(sample, into=c("rounds", "polymerase", "dummy1", "dummy2", "subject")) %>%
	mutate(rounds = as.numeric(str_replace(rounds, "x", ""))) %>%
	select(-dummy1, -dummy2) %>%
	write_tsv("data/process/stool_chimera.tsv")




read_tsv("data/process/stool_chimera.tsv") %>%
	ggplot(aes(x=rounds, y=frac_chimera, color=subject)) +
		geom_line() +
		facet_grid(.~polymerase) +
		ggsave("stool_chimera.pdf")
