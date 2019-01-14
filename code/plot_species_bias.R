library(tidyverse)
library(RColorBrewer)

mock <- read_tsv("data/process/mock_bias_species.tsv")

polymerase_colors <- brewer.pal(5, "RdBu")
polymerase_colors[3] <- "darkgray"

species_order <- mock %>% arrange(percent_gc) %>% pull(species) %>% str_replace("_", " ") %>% unique

mock %>%
	mutate(species = factor(str_replace(species, "_", " "), levels=species_order),
				obs_rel_abund = 100 * obs_rel_abund,
				actual_rel_abund = 100 * actual_rel_abund,
				) %>%
	ggplot(aes(x=rounds, y=obs_rel_abund, group=polymerase, color=polymerase)) +
		facet_wrap(~species)+#, labeller = label_parsed) +
		geom_hline(aes(yintercept=actual_rel_abund), col="gray") +
		geom_line(size=1) +
		geom_text(aes(x=20.75, y=0.5, label=format(percent_gc, nsmall=1L)), size=3, color="black") +
		coord_cartesian(ylim=c(0,25), xlim=c(20,35)) +
		labs(y="Relative abundance (%)", x="Number of rounds of PCR") +
		scale_color_manual(name=NULL,
												breaks=c("ACC", "K", "PHU", "PL", "Q5"),
												labels=c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"),
												values=polymerase_colors) +
		theme_classic() +
		theme(strip.text=element_text(face="italic"),
					strip.text.x=element_text(hjust=0, size=9),
					strip.background = element_rect(color = 'white', fill = 'white', size = 10),
					legend.position=c(0.85, 0.1),
					legend.key.height = unit(0.8, "line")
				)

ggsave("results/figures/species_bias.pdf", width=6.875, height=6)
