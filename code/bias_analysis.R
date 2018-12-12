library(tidyverse)

# see data/references/ZRC187325.pdf for certificate of analysis from Zymo
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


# mock_rel_abund <- read_tsv("data/mothur/zymo_mock.filter.pick.count_table") %>%
# 	rename(reference=Representative_Sequence) %>%
# 	mutate(species = str_replace(reference, "_16S_.*", "")) %>%
# 	group_by(species) %>%
# 	mutate(total_copies = sum(total)) %>%
# 	ungroup() %>%
# 	mutate(frac_genome_copies = total/total_copies) %>%
# 	inner_join(., measured_composition, by="species") %>%
# 	mutate(exp_rel_abund = actual_rel_abund * frac_genome_copies) %>%
# 	select(reference, species, exp_rel_abund)


fasta <- scan("data/mothur/zymo_mock.filter.pick.fasta", quiet=T, what=character())

mock_gc <- tibble(reference = str_replace(fasta[c(T,F)], ">", ""),
										sequences = str_replace_all(fasta[c(F,T)], "[-.]", "")) %>%
						mutate(species = str_replace(reference, "_16S_.*", ""),
									length = nchar(sequences),
									gc = nchar(sequences) - nchar(str_replace_all(sequences, "[GC]", "")),
									percent_gc = gc / length
									) %>%
						group_by(species) %>%
						summarize(percent_gc = round(100 * mean(percent_gc), digits=1))


mock <- inner_join(mock_rel_abund, mock_gc, by="species") %>%
				arrange(species)



contig_error_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary"
contig_errors <- read_tsv(contig_error_file_name) %>%
						mutate(species = str_replace(reference, "_16S_.*", "")) %>%
						select(query, reference, species, mismatches, total, numparents) %>%
						rename(seq_name = query)

contig_count_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table"
contig_counts <- read_tsv(contig_count_file_name) %>%
						select(-total) %>%
						rename(seq_name = Representative_Sequence) %>%
						gather(sample, count, -seq_name) %>%
						filter(count != 0)

# This dataframe contains the average error rate for the sequences after running pre.cluster on
# the data. We have pooled the four replicates to get a total average error rate.
rel_abundances <- inner_join(contig_errors, contig_counts, by="seq_name") %>%
			separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
			mutate(rounds = str_replace(rounds, "x", "")) %>%
			filter(numparents == 1, mismatches <= 20) %>%
			select(species, rounds, polymerase, dna, count) %>%
			group_by(polymerase, rounds) %>%
			mutate(n_seqs = sum(count)) %>%
			filter(n_seqs > 200) %>%
			group_by(polymerase, rounds, species) %>%
			summarize(obs_count = sum(count), obs_rel_abund = sum(count) / n_seqs[1]) %>%
			ungroup() %>%
			arrange(species)


rel_abundances %>%
	select(polymerase, rounds, species, obs_count) %>%
	spread(species, obs_count) %>%
	mutate(label="1", Group=paste(polymerase, rounds, sep="_"), numOtus=8) %>%
	select(label, Group, numOtus, everything(), -polymerase, -rounds, ) %>%
	write_tsv("data/mothur/taxa_mapping.shared")


inner_join(rel_abundances, mock, by="species") %>%
	mutate(species = str_replace(species, "_", " ")) %>%
	mutate(species = factor(species, levels= mock %>% arrange(percent_gc) %>% pull(species) %>% str_replace("_", " "))) %>%
	ggplot(aes(x=rounds, y=obs_rel_abund, group=polymerase, color=polymerase)) +
		facet_wrap(~species)+#, labeller = label_parsed) +
		geom_hline(aes(yintercept=actual_rel_abund), col="gray") +
		geom_line() +
		geom_text(aes(x=1, y=0.25, label=format(percent_gc, nsmall=1L)), size=3, color="black") +
		coord_cartesian(ylim=c(0,0.25)) +
		theme_classic() +
		theme(strip.text=element_text(face="italic"),
					strip.text.x=element_text(hjust=0, size=9),
					strip.background = element_rect(color = 'white', fill = 'white', size = 10)) +
		ggsave("results/figures/mock_bias.pdf")






salmonella_ratio <- inner_join(contig_errors, contig_counts, by="seq_name") %>%
		separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
		mutate(rounds = str_replace(rounds, "x", "")) %>%
		filter(numparents == 1, mismatches <= 20) %>%
		select(reference, rounds, polymerase, dna, count) %>%
		group_by(polymerase, rounds) %>%
		mutate(n_seqs = sum(count)) %>%
		filter(n_seqs > 200) %>%
		group_by(polymerase, rounds, reference) %>%
		summarize(obs_rel_abund = sum(count) / n_seqs[1]) %>%
		ungroup() %>%
		filter(str_detect(reference, "Salmonella_enterica")) %>%
		group_by(polymerase, rounds) %>%
		summarize(ratio = obs_rel_abund[2]/obs_rel_abund[1]) %>%
		ungroup()

ggplot(salmonella_ratio, aes(x=rounds, y=ratio, group=polymerase, color=polymerase)) +
	geom_hline(aes(yintercept=6), col="gray") +
	geom_line() +
	theme_classic() +
	ggsave("results/figures/mock_salmonella_ratio.pdf")
