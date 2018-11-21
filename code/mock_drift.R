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

perfect_dist_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.braycurtis.0.03.lt.ave.dist"

perfect_dist <- read_dist(perfect_dist_file) %>% mutate(method="perfect")

vsearch_dist_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.dist"

vsearch_dist <- read_dist(vsearch_dist_file) %>% mutate(method="vsearch")

bind_rows(perfect_dist, vsearch_dist) %>%
	mutate(rows = str_replace(rows, "_Mock_DNA.", ""),
				columns = str_replace(columns, "_Mock_DNA.", "")) %>%
	filter(rows == columns) %>%
	group_by(rows, method) %>%
	summarize(ave_dist = mean(distances)) %>%
	ungroup() %>%
	separate(rows, into=c("rounds", "polymerase")) %>%
	mutate(rounds = as.numeric(str_replace(rounds, "x", ""))) %>%
	write_csv("data/process/mock_beta_drift.csv")

read_csv("data/process/mock_beta_drift.csv") %>%
	filter(method == "vsearch") %>%
	ggplot(aes(x=rounds, y=ave_dist, color=polymerase)) +
		geom_line() +
		ggsave("test.pdf")
