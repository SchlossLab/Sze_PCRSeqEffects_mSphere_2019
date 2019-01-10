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
