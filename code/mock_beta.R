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


read_tsv("data/process/mock_beta_diversity.tsv") %>%
	filter(str_detect(rows, "25x")) %>%
	separate(rows, into=c("row_rounds", "polymerase", "row_mock")) %>%
	separate(columns, into=c("col_rounds", "col_polymerase", "col_mock")) %>%
	filter(polymerase == col_polymerase, method=="vsearch") %>%
	select(polymerase, col_rounds, distances, method) %>%
	arrange(polymerase) %>%
	mutate(col_rounds = str_replace(col_rounds, "x", "")) %>%
	ggplot(aes(x=polymerase, y=distances, fill=col_rounds))+
		geom_col(position=position_dodge()) +
		ggsave("results/figures/mock_beta.pdf")
