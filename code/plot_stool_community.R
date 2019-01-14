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


alpha <- read_tsv("data/process/stool_alpha_diversity.tsv") %>%
	filter(metric == "shannon" | metric == "sobs") %>%
	mutate(metric = factor(metric, levels = c("sobs", "shannon"))) %>%
	ggplot(aes(x=rounds, y=value, group=subject, color=subject)) +
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
		# breaks=c("no_error", "perfect", "vsearch"),
		# labels=c("No errors or\nchimeras", "Errors and\nperfect chimera\nremoval", "Errors and\nchimera removal\nwith VSEARCH"),
		values=brewer.pal(4, "Set2")
	) +
	theme_classic() +
	theme(
		panel.spacing=unit(0.5,"cm"),
		strip.background = element_blank(),
		legend.key.height = unit(0.8, "line"),
		strip.placement="outside",
		strip.text=element_text(size=11)
	)


rel_distance <- read_tsv("data/process/stool_beta_diversity.tsv") %>%
 	filter(str_detect(rows, "25x")) %>%
 	separate(rows, into=c("row_rounds", "polymerase", "subject")) %>%
 	separate(columns, into=c("col_rounds", "col_polymerase", "col_subject")) %>%
 	filter(polymerase == col_polymerase, subject == col_subject) %>%
 	select(polymerase, col_rounds, distances, subject) %>%
	arrange(polymerase) %>%
	mutate(col_rounds = str_replace(col_rounds, "x", "")) %>%
	ggplot(aes(x=polymerase, y=distances, col=col_rounds, group=paste(subject, polymerase))) +
		geom_line(position=position_dodge(width=0.5), color="black")+
		geom_point(position=position_dodge(width=0.5)) +
		labs(x=NULL, y="Bray-Curtis distance relative\nto 25 rounds of PCR") +
		scale_x_discrete(labels=polymerase_names) +
		scale_color_manual(name=NULL, breaks=c(30, 35), labels=c("30 rounds", "35 rounds"), values=c("blue", "red")) +
		coord_cartesian(ylim=c(0,0.35)) +
		theme_classic() +
		theme(
			legend.position=c(0.22, 0.20),
			axis.text.x=element_text(size=8),
			legend.key.height = unit(0.8, "line"),
			legend.key.size = unit(0.1, "cm")
		)


label_positions <- tibble(x=c(-0.3, 0.095, 0.03, 0.12), y=c(0.05, -0.2, 0.2, -0.04), round = c('DA10016', 'DA10029', 'DA10001', 'DA10040'), label=c('DA10016', 'DA10029', 'DA10001', 'DA10040'))

pcoa <- read_tsv("data/mothur/stool.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
	select(group, axis1, axis2) %>%
	separate(group, into=c("rounds", "polymerase", "pmm", "FS", "subject")) %>%
	mutate(rounds = str_replace(rounds, "x", "")) %>%
	ggplot(aes(x=axis1, y=axis2, color=rounds, shape=polymerase)) +
		stat_ellipse(aes(group=subject), type="norm", geom="polygon", fill="#DDDDDD", color="#DDDDDD", show.legend=FALSE) +
		geom_point() +
		geom_text(data = label_positions, aes(x, y, label=label), inherit.aes=FALSE, show.legend=FALSE, size=3) +
		scale_shape(name=NULL, labels=polymerase_names) +
		scale_color_manual(name=NULL, labels=c("20 rounds", "25 rounds", "30 rounds", "35 rounds"), values=brewer.pal(4, "Set1")) +
		guides(fill=FALSE) +
		labs(x="PCo Axis 1", y="PCo Axis 2") +
		theme_classic() +
		theme(
			legend.box="horizontal",
			legend.position = c(0.28, 0.18),
			legend.key.height = unit(0.8, "line"),
			legend.key.size = unit(0.06, "cm"),
			legend.margin = margin(0,0,0,0)
		)

bottom_row <- plot_grid(rel_distance, pcoa, labels=c("B", "C"), rel_widths=c(1,1.2))
plot_grid(alpha, bottom_row, nrow=2, labels=c("A", ""))

ggsave("results/figures/stool_community.pdf", width=6.875)
