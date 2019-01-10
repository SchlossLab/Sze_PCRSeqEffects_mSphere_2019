library(tidyverse)
library(viridis)
library(RColorBrewer)

library(cowplot)

facet_names <- c(
									'percentage' = "Total chimeras (%)",
									'sensitivity' = "Sensitivity (%)",
									'specificity' = "Specificity (%)"
								)

mock <- read_tsv("data/process/error_chimera_rates.tsv") %>%
	mutate(percentage = 100 * pc_frac_chimera,
					sensitivity = 100 * vsearch_sens,
					specificity = 100 * vsearch_spec) %>%
	select(rounds, polymerase, percentage, sensitivity, specificity) %>%
	gather(per_sens_spec, value, percentage, sensitivity, specificity) %>%
	ggplot(aes(x=rounds, y=value, color=polymerase)) +
		geom_line(size=1) +
		geom_blank(aes(y = 0)) +
		facet_wrap(~per_sens_spec, scales="free_y", labeller = as_labeller(facet_names), strip.position="left") +
		labs(y=NULL, x="Number of rounds of PCR") +
		scale_color_manual(name=NULL,
												breaks=c("ACC", "K", "PHU", "PL", "Q5"),
												labels=c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"),
												values=viridis(5)) +
		theme_classic() +
		theme(
			panel.spacing=unit(0.5,"cm"),
			strip.background = element_blank(),
			legend.key.height = unit(0.8, "line"),
			strip.placement="outside",
			strip.text=element_text(size=11)
		)



facet_names <- c(
									'ACC' = "Accuprime",
									'K' = "Kappa",
									'PHU' = "Phusion",
									'PL' = "Platinum",
									'Q5' = "Q5"
								)

stool <- read_tsv("data/process/stool_chimera.tsv") %>%
	mutate(per_chimera = 100 * frac_chimera) %>%
	ggplot(aes(x=rounds, y=per_chimera, color=subject)) +
		geom_line(size=1) +
		facet_grid(.~polymerase, labeller = as_labeller(facet_names)) +
		scale_color_manual(name=NULL,
												values=brewer.pal(5, "Set1")[-3]) +
		labs(y="Detected chimeras (%)", x="Number of rounds of PCR") +
		theme_classic() +
		theme(
			panel.spacing=unit(0.5,"cm"),
			strip.background = element_blank(),
			legend.key.height = unit(0.8, "line")
		)

plot_grid(mock, stool, nrow=2, labels=c("A", "B"))

ggsave("results/figures/chimera_plots.pdf", width=6.875, height=4.5, units="in")
