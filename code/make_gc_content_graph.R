### Plot GC content graph
# Marc Sze


## Load in needed functions 
source("code/functions.R")

## Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "scales", "viridis", "RColorBrewer", "ggrepel"))

# Load in needed data
summary_data <- read_csv("data/process/tables/gc_summary_table.csv")

# summary_data <- read_csv("data/process/tables/gc_content_amp_summary.csv") %>% 
#   filter(cycles != "fifteen")

gc_graph <- summary_data %>% 
  mutate(cycles = factor(cycles, 
                         levels = c("twenty", "twenty_five", "thirty", "thirty_five"), 
                         labels = c("20x", "25x", "30x", "35x")), 
         taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")), 
         reference = factor(reference, 
                            levels = c("Staphylococcus_aureus_16S", "Listeria_monocytogenes_16S", 
                                       "Enterococcus_faecalis_16S", "Bacillus_subtilis_16S", 
                                       "Lactobacillus_fermentum_16S", "Escherichia_coli_16S_1", 
                                       "Salmonella_enterica_16S_1", "Pseudomonas_aeruginosa_16S"), 
                            labels = c("Staphylococcus aureus\n(GC = 32.7%)", "Listeria monocytogenes\n(GC = 38.0%)", 
                                       "Enterococcus faecalis\n(GC = 38.0%)", "Bacillus subtilis\n(GC = 43.6%)", 
                                       "Lactobacillus fermentum\n(GC = 49.7%)", "Escherichia coli\n(GC = 50.6%)", 
                                       "Salmonella enterica\n(GC = 52.1%)", "Pseudomonas aeruginosa\n(GC = 66.6%)"))) %>% 
  ggplot(aes(cycles, median_abund, color = reference, group = reference)) + 
  geom_hline(yintercept = 0.12, linetype = "dashed") + 
  geom_pointrange(aes(ymin = min_abund, ymax = max_abund)) + 
  geom_line() + theme_bw() + coord_cartesian(ylim = c(0, 0.3)) + 
  labs(x = "", y = "Relative Abundance") + 
  facet_wrap(~taq) + 
  scale_color_brewer(name = "Bacterial Species", type = "div", palette = "RdBu") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        legend.position = c(0.85, 0.3), 
        legend.text = element_text(size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("results/figures/Figure1.pdf", gc_graph, width = 6, height = 8, dpi = 300)






