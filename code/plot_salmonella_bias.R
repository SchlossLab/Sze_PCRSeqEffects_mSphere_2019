library(tidyverse)
library(RColorBrewer)

polymerase_colors <- brewer.pal(5, "RdBu")
polymerase_colors[3] <- "darkgray"

read_tsv("data/process/mock_bias_salmonella.tsv") %>%
	ggplot(aes(x=rounds, y=ratio, group=polymerase, color=polymerase)) +
		geom_hline(aes(yintercept=6), col="gray") +
		geom_line(size=1) +
		coord_cartesian(ylim=c(4,8)) +
		labs(y="Ratio of relative abundances\nfor each operon", x="Number of rounds of PCR") +
		scale_color_manual(name=NULL,
												breaks=c("ACC", "K", "PHU", "PL", "Q5"),
												labels=c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"),
												values=polymerase_colors) +
		theme_classic() +
		theme(
			legend.key.height = unit(0.8, "line")
		)

ggsave("results/figures/salmonella_bias.pdf", width=4.4, height=3)
