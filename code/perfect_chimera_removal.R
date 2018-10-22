####################################################################################################
#
#	perfect_chimera_removal.R
#
#	This script will take in a fasta, count_table, and error.summary file and remove those sequences
#	that are flagged as having two parents
#
####################################################################################################

library(tidyverse)

# Read in the error.summary file and identify those sequences with max_diffs or fewer mismatches to
# a reference sequence
error_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.error.summary"
good_seqs <- read_tsv(error_file_name) %>%
								filter(numparents == 1) %>%
								select(query)


# Read in the fasta file and select for those sequences that were identified in good_seqs
fasta_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.fasta"
new_fasta_file_name <- str_replace(fasta_file_name, "fasta", "perfect.fasta")

fasta_file <- scan(fasta_file_name, sep="\n", what=character(), quiet=T)
tibble(query = fasta_file[c(TRUE, FALSE)], sequence = fasta_file[c(FALSE, TRUE)]) %>%
	mutate(query = str_replace(query, ">(.*)\\t$", "\\1")) %>%
	inner_join(., good_seqs, by="query") %>%
	mutate(output = paste0(">", query, "\n", sequence)) %>%
	pull(output) %>%
	write(., new_fasta_file_name)


# Read in the count_table file and select for those sequences that were identified in good_seqs
count_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.count_table"
new_count_file_name <- str_replace(count_file_name, "count_table", "perfect.count_table")

read_tsv(count_file_name) %>%
	inner_join(., good_seqs, by=c("Representative_Sequence" = "query")) %>%
	write_tsv(new_count_file_name)
