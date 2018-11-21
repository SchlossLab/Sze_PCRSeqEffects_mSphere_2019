library(tidyverse)


read_dist <- function(dist_file_name){
	linear_data <- scan(dist_file_name, what="character", quiet=TRUE)[-1]

	samples <- str_subset(linear_data, "_")
	n_samples <- length(samples)
	distance_strings <- str_subset(linear_data, "\\.")

	distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
	colnames(distance_matrix) <- samples
	as.tibble(cbind(rows=samples, distance_matrix)) %>%
		gather(columns, distances, -rows) %>%
		filter(rows < columns) %>%
		arrange(columns, rows) %>%
		mutate(distances = as.numeric(distance_strings))
}

dist_file <- "data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist"

read_dist(dist_file) %>%
	mutate(rows = str_replace(rows, "_PMM_FS", ""),
				columns = str_replace(columns, "_PMM_FS", "")) %>%
	write_tsv("data/process/stool_beta_diversity.tsv")


read_tsv("data/process/stool_beta_diversity.tsv") %>%
	filter(str_detect(rows, "25x")) %>%
	separate(rows, into=c("row_rounds", "polymerase", "subject")) %>%
	separate(columns, into=c("rounds", "col_polymerase", "col_subject")) %>%
	filter(polymerase == col_polymerase, subject == col_subject) %>%
	select(polymerase, rounds, subject, distances) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=subject, y=distances, fill=rounds))+
		geom_col(position=position_dodge()) +
		facet_wrap(.~polymerase) +
		ggsave("test1.pdf")


read_tsv("data/process/stool_beta_diversity.tsv") %>%
	separate(rows, into=c("rounds", "polymerase", "row_subject")) %>%
	separate(columns, into=c("col_rounds", "col_polymerase", "col_subject")) %>%
	filter(rounds == col_rounds, polymerase == col_polymerase) %>%
	select(rounds, polymerase, row_subject, col_subject, distances) %>%
	group_by(rounds, polymerase) %>%
	summarize(ave_dist = mean(distances), n = n()) %>%
	ungroup() %>%
	filter(n == 6) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=rounds, y=ave_dist, group=polymerase, color=polymerase)) +
		geom_line() +
		ggsave('test2.pdf')
