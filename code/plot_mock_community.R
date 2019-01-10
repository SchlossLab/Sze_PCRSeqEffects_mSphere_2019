library(tidyverse)
library(cowplot)
library(RColorBrewer)

polymerase_names <- c(
									'ACC' = "Accuprime",
									'K' = "Kappa",
									'PHU' = "Phusion",
									'PL' = "Platinum",
									'Q5' = "Q5"
								)

alpha_names <- c(
									'sobs' = "OTU Richness",
									'shannon' = "Shannon\nDiversity",
									'invsimpson' = "Inverse Simpson\nDiversity"
								)


alpha <- read_tsv("data/process/mock_alpha_diversity.tsv") %>%
	filter(metric == "shannon" | metric == "sobs") %>%
	mutate(metric = factor(metric, levels = c("sobs", "shannon"))) %>%
	ggplot(aes(x=rounds, y=value, group=method, color=method)) +
	geom_hline(aes(yintercept=ideal)) +
	geom_line(size=1) +
	facet_grid(metric~polymerase,
						scales="free_y",
						switch="y",
						labeller = labeller(.cols = polymerase_names,
																.rows = alpha_names)
						) +
	labs(y=NULL, x="Number of rounds of PCR") +
	scale_color_manual(
		name=NULL,
		breaks=c("no_error", "perfect", "vsearch"),
		labels=c("No errors or\nchimeras", "Errors and\nperfect chimera\nremoval", "Errors and\nchimera removal\nwith VSEARCH"),
		values=brewer.pal(3, "Set2")
	) +
	theme_classic() +
	theme(
		panel.spacing=unit(0.5,"cm"),
		strip.background = element_blank(),
		legend.key.height = unit(2.3, "line"),
		strip.placement="outside",
		strip.text=element_text(size=11)
	)

# ggsave("results/figures/mock_alpha.pdf", width=6.875, height=3)


rel_distance <- read_tsv("data/process/mock_beta_diversity.tsv") %>%
	filter(str_detect(rows, "25x")) %>%
	separate(rows, into=c("row_rounds", "polymerase", "row_mock")) %>%
	separate(columns, into=c("col_rounds", "col_polymerase", "col_mock")) %>%
	filter(polymerase == col_polymerase, method=="vsearch") %>%
	select(polymerase, col_rounds, distances, method) %>%
	arrange(polymerase) %>%
	mutate(col_rounds = str_replace(col_rounds, "x", "")) %>%
	ggplot(aes(x=polymerase, y=distances, fill=col_rounds))+
		geom_col(position=position_dodge()) +
		labs(x=NULL, y="Bray-Curtis Distance") +
		scale_x_discrete(labels=polymerase_names) +
		scale_fill_manual(name=NULL, breaks=c(30, 35), labels=c("30 rounds", "35 rounds"), values=c("blue", "red")) +
		theme_classic() +
		theme(
			legend.position=c(0.22, 0.90),
			axis.text.x=element_text(size=8)
		)


label_positions <- tibble(x=c(-0.08, 0, 0.08, 0.16), y=c(0.12, 0.086, 0.088, 0.04), round = c('20', '25', '30', '35'), label=c('20 rounds', '25 rounds', '30 rounds', '35 rounds'))

pcoa <- read_tsv("data/mothur/mock.trim.contigs.good.unique.good.filter.unique.pick.pick.precluster.vsearch.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
	select(group, axis1, axis2) %>%
	separate(group, into=c("rounds", "polymerase", "mock", "dna")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=axis1, y=axis2, color=rounds, shape=polymerase)) +
		geom_point() +
		stat_ellipse(aes(group=rounds), show.legend=FALSE) +
		geom_text(data = label_positions, aes(x, y, label=label, color=round), inherit.aes=FALSE, show.legend=FALSE, size=3) +
		scale_shape(name=NULL, labels=polymerase_names) +
		scale_color_manual(name=NULL, values=brewer.pal(4, "Set1")) +
		guides(color=FALSE) +
		labs(x="PCo Axis 1", y="PCo Axis 2") +
		theme_classic() +
		theme(
			legend.position = c(0.88, 0.9),
			legend.key.height = unit(0.8, "line"),
			legend.key.size = unit(0.1, "cm")
		)

bottom_row <- plot_grid(rel_distance, pcoa, labels=c("B", "C"), rel_widths=c(1,1.2))
plot_grid(alpha, bottom_row, nrow=2, labels=c("A", ""))

ggsave("results/figures/mock_community.pdf", width=6.875)
