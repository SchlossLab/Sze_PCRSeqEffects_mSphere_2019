library(tidyverse)
library(viridis)

read_tsv("data/process/mock_bias_salmonella.tsv") %>%
	ggplot(aes(x=rounds, y=ratio, group=polymerase, color=polymerase)) +
		geom_hline(aes(yintercept=6), col="gray") +
		geom_line() +
		coord_cartesian(ylim=c(4,8)) +
		labs(y="Ratio of relative abundances for each operon", x="Number of rounds of PCR") +
		scale_color_manual(name=NULL,
												breaks=c("ACC", "K", "PHU", "PL", "Q5"),
												labels=c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"),
												values=viridis(5)) +
		theme_classic()

ggsave("results/figures/salmonella_bias.pdf", width=6.875, height=6)
