library(tidyverse)

contig_error_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.error.summary"
contig_errors <- read_tsv(contig_error_file_name) %>%
						select(query, reference, mismatches, total, numparents) %>%
						rename(seq_name = query)

contig_count_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.pick.count_table"
contig_counts <- read_tsv(contig_count_file_name) %>%
						select(-total) %>%
						rename(seq_name = Representative_Sequence) %>%
						gather(sample, count, -seq_name) %>%
						filter(count != 0)

# This dataframe contains the average error rate for the sequences after forming contigs. We have
# pooled the four replicates to get a total average error rate. We also toss any polymerase/rounds
# that had fewer than 50 sequences on average
contig_error_summary <- inner_join(contig_errors, contig_counts, by="seq_name") %>%
	separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
	select(-mock, -dna) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(numparents == 1, mismatches <= 20) %>%
	group_by(rounds, polymerase) %>%
	summarize(contig_error = sum(count*mismatches)/sum(count*total),
						contig_n = sum(count)) %>%
	filter(contig_n > 200) %>%
	ungroup()


# This dataframe contains the fraction of sequences that were flagged as being chimeras after
# forming contigs. We have pooled the four replicates to get an overall chimera rate.
contig_chimera_summary <- inner_join(contig_errors, contig_counts, by="seq_name") %>%
	separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
	select(-mock, -dna) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(mismatches <= 20) %>%
	group_by(rounds, polymerase) %>%
	summarize(contig_frac_chimera = sum(count * (numparents == 2))/sum(count)) %>%
	ungroup()



pc_error_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.error.summary"
pc_errors <- read_tsv(pc_error_file_name) %>%
						select(query, reference, mismatches, total, numparents) %>%
						rename(seq_name = query)

pc_count_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table"
pc_counts <- read_tsv(pc_count_file_name) %>%
						select(-total) %>%
						rename(seq_name = Representative_Sequence) %>%
						gather(sample, count, -seq_name) %>%
						filter(count != 0)

# This dataframe contains the average error rate for the sequences after running pre.cluster on
# the data. We have pooled the four replicates to get a total average error rate.
pc_error_summary <- inner_join(pc_errors, pc_counts, by="seq_name") %>%
	separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
	select(-mock, -dna) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(numparents == 1, mismatches <= 20) %>%
	group_by(rounds, polymerase) %>%
	summarize(pc_error = sum(count*mismatches)/sum(count*total)) %>%
	ungroup()

# This dataframe contains the fraction of sequences that were flagged as being chimeras after
# running pre.cluster. We have pooled the four replicates to get an overall chimera rate.
pc_chimera_summary <- inner_join(pc_errors, pc_counts, by="seq_name") %>%
	separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
	select(-mock, -dna) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	filter(mismatches <= 20) %>%
	group_by(rounds, polymerase) %>%
	summarize(pc_frac_chimera = sum(count * (numparents == 2))/sum(count)) %>%
	ungroup()


chimera_screened_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.count_table"
screened_counts <- read_tsv(chimera_screened_file_name) %>%
						select(-total) %>%
						rename(seq_name = Representative_Sequence) %>%
						gather(sample, count, -seq_name)


# This dataframe contains the sensitivity and specificity for each polymerase and round number for
# the vsearch implementation of uchime, which was run after the pre.cluster step. We have pooled
# the four replicates to get an overall sensitivities and specificities.
sens_spec <- inner_join(pc_errors, pc_counts, by="seq_name") %>%
						left_join(., screened_counts, by=c("seq_name", "sample")) %>%
						separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
						select(-mock, -dna) %>%
						mutate(rounds = str_replace(rounds, "x", "")) %>%
						filter(mismatches <= 20) %>%
						mutate(tp = (numparents == 2) & is.na(count.y),
									fp = (numparents == 1) & is.na(count.y),
									tn = (numparents == 1) & !is.na(count.y),
									fn = (numparents == 2) & !is.na(count.y)) %>%
						group_by(rounds, polymerase) %>%
						summarize(tp = sum(tp), fp = sum(fp), tn = sum(tn), fn = sum(fn),
											vsearch_sens = tp/(tp+fn), vsearch_spec = tn/(fp+tn)) %>%
						ungroup() %>%
						select(-tp, -fp, -fn, -tn)

# Want to see whether all of the sequences made it into the mock community data
n_references <- inner_join(contig_errors, contig_counts, by="seq_name") %>%
						separate(sample, into=c("rounds", "polymerase", "mock", "dna")) %>%
						select(-mock, -dna) %>%
						mutate(rounds = str_replace(rounds, "x", "")) %>%
						filter(numparents == 1, mismatches <= 20) %>%
						group_by(rounds, polymerase, reference) %>%
						summarize(N = sum(count)) %>%
						summarize(n_references = n()) %>% print(n=50)
					# -> they did.

#	Aggregate everything and spit it out...
inner_join(contig_error_summary, contig_chimera_summary, by=c("rounds", "polymerase")) %>%
	inner_join(., pc_error_summary, by=c("rounds", "polymerase")) %>%
	inner_join(., pc_chimera_summary, by=c("rounds", "polymerase")) %>%
	inner_join(., sens_spec, by=c("rounds", "polymerase")) %>%
	inner_join(., n_references, by=c("rounds", "polymerase")) %>%
	write_tsv("data/process/error_chimera_rates.tsv")
