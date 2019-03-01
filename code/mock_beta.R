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

perfect_dist_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.perfect.opti_mcc.merge.braycurtis.0.03.lt.ave.dist"

perfect_dist <- read_dist(perfect_dist_file) %>% mutate(method="perfect")

vsearch_dist_file <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.merge.braycurtis.0.03.lt.ave.dist"

vsearch_dist <- read_dist(vsearch_dist_file) %>% mutate(method="vsearch")


bind_rows(perfect_dist, vsearch_dist) %>%
	write_tsv("data/process/mock_beta_diversity.tsv")
