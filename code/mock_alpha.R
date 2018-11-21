library(tidyverse)

perfect_chimera_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.groups.ave-std.summary"

perfect_chimera <- read_tsv(perfect_chimera_file) %>%
	separate(group, into=c("rounds", "polymerase", "mock")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(method == 'ave') %>%
	mutate(method="perfect") %>%
	select(method, rounds, polymerase, nseqs, coverage, invsimpson, shannon, sobs)


vsearch_chimera_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.groups.ave-std.summary"

vsearch_chimera <- read_tsv(vsearch_chimera_file) %>%
	separate(group, into=c("rounds", "polymerase", "mock")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(method == 'ave') %>%
	mutate(method="vsearch") %>%
	select(method, rounds, polymerase, nseqs, coverage, invsimpson, shannon, sobs)


ideal_list_file <- "data/mothur/zymo_mock.filter.pick.unique.precluster.opti_mcc.list"
ideal_list <- readLines(ideal_list_file)[2] %>% str_split("\t") %>% unlist() %>% str_subset("_")

ideal_count_table_file <- "data/mothur/zymo_mock.filter.pick.unique.precluster.count_table"
ideal_count_table <- read_tsv(ideal_count_table_file)

stopifnot(nrow(ideal_count_table) == length(ideal_list))

mock_rel_abund <-
		tibble(
			species = c("Pseudomonas_aeruginosa", "Escherichia_coli", "Salmonella_enterica",
				"Lactobacillus_fermentum", "Enterococcus_faecalis", "Staphylococcus_aureus",
				"Listeria_monocytogenes", "Bacillus_subtilis"),
			# ideal_abund = 12.5,
			actual_abund = c(10.4, 9.0, 11.8, 10.3, 14.1, 14.6, 13.2, 13.2),
			genome_size = c(6.792, 4.875, 4.760, 1.905, 2.845, 2.730, 2.992, 4.045),
			copy_number = c(4, 7, 7, 5, 4, 6, 6, 10),
			rrn_theor_abund = c(4.2, 10.1, 10.4, 18.4, 9.9, 15.5, 14.1, 17.4)/100,
			genome_gc = c(66.2, 46.7, 52.2, 52.4, 37.5, 32.9, 38, 43.9)
		) %>%
		mutate(ideal_rel_abund = copy_number / genome_size,
					ideal_rel_abund = ideal_rel_abund/sum(ideal_rel_abund),
					actual_rel_abund = actual_abund * copy_number / genome_size,
					actual_rel_abund = actual_rel_abund/sum(actual_rel_abund)) %>%
		select(-rrn_theor_abund, -ideal_rel_abund)

ideal <- mock_rel_abund %>%
	mutate(scaled_rel_abund = actual_rel_abund * 1000) %>%
	summarize(sobs = nrow(.),
				shannon = -sum(actual_rel_abund * log(actual_rel_abund)),
				simpson = (sum(scaled_rel_abund * (scaled_rel_abund-1)) /
									(sum(scaled_rel_abund) * (sum(scaled_rel_abund)-1))),
				invsimpson = 1/simpson) %>%
	select(-simpson) %>%
	gather(metric, ideal)



bind_rows(vsearch_chimera, perfect_chimera) %>%
	gather(key=metric, value=value, nseqs, coverage, invsimpson, shannon, sobs) %>%
	inner_join(., ideal, by="metric") %>%
	write_tsv("data/process/mock_alpha_diversity.tsv")




read_tsv("data/process/mock_alpha_diversity.tsv") %>%
	filter(metric== "invsimpson" | metric == "shannon" | metric == "sobs") %>%
	mutate(metric = factor(metric, levels = c("sobs", "shannon", "invsimpson"))) %>%
	ggplot(aes(x=rounds, y=value, group=method, color=method)) +
	geom_line() +
	geom_hline(aes(yintercept=ideal)) +
	facet_grid(metric~polymerase, scales="free_y") +
	ggsave("results/figures/mock_alpha.pdf")
