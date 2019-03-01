library(tidyverse)
library(vegan)

set.seed(19760620) #pds birthday

####################################################################################################

read_dist <- function(dist_file_name){

	linear_data <- scan(dist_file_name, what="character", sep="\n", quiet=TRUE)

	n_samples <- as.numeric(linear_data[1])
	linear_data <- linear_data[-1]

	samples <- str_replace(linear_data, "\t.*", "")
	linear_data <- str_replace(linear_data, "[^\t]*\t", "")
	linear_data <- linear_data[-1]

	distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)

	for(i in 1:(n_samples-1)){
		row <- as.numeric(unlist(str_split(linear_data[i], "\t")))
		distance_matrix[i+1,1:length(row)] <- row
	}

	distance_matrix <- distance_matrix + t(distance_matrix)
	rownames(distance_matrix) <- samples

	as.dist(distance_matrix)
}

####################################################################################################

mock_file_name <- "data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.dist"

mock_dist <- read_dist(mock_file_name)

mock_samples <- tibble(names = attr(mock_dist, "Labels")) %>% mutate(names=str_replace(names, "_Mock_.*", "")) %>% separate(names, into=c("rounds", "polymerase"))

mock_adonis <- adonis(mock_dist~mock_samples$rounds+mock_samples$polymerase, permutations=9999)

tibble(effects = c("rounds", "polymerase"),
			r_sq = mock_adonis$aov.tab$R2[1:2],
			p = mock_adonis$aov.tab$Pr[1:2]) %>%
	write_tsv("data/process/vegan_mock.tsv")

####################################################################################################

stool_file_name <- "data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.dist"

stool_dist <- read_dist(stool_file_name)

stool_samples <- tibble(names = attr(stool_dist, "Labels")) %>% mutate(names=str_replace(names, "_PMM_FS", "")) %>% separate(names, into=c("rounds", "polymerase", "subject"))

stool_adonis <- adonis(stool_dist~stool_samples$subject+stool_samples$rounds+stool_samples$polymerase, permutations=9999)

tibble(effects = c("subject", "rounds", "polymerase"),
			r_sq = stool_adonis$aov.tab$R2[1:3],
			p = stool_adonis$aov.tab$Pr[1:3]) %>%
	write_tsv("data/process/vegan_stool.tsv")


####################################################################################################
