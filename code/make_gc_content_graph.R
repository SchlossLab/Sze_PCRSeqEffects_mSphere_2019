### Plot RF results MDA for classifying samples
# Marc Sze


## Load in needed functions 
source("code/functions.R")

## Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "scales", "viridis", "RColorBrewer", "ggrepel"))

# Load in needed data
summary_data <- read_csv("data/process/tables/gc_content_amp_summary.csv")


gc_graph <- summary_data %>% 
  mutate(cycles = factor(cycles, 
                         levels = c("fifteen", "twenty", "twenty_five", "thirty", "thirty_five"), 
                         labels = c("15x", "20x", "25x", "30x", "35x")), 
         higher_gc = factor(higher_gc, 
                            levels = c("yes", "no"), 
                            labels = c("Yes", "No")), 
         taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log10(median_of_median_abund+1), color = higher_gc)) + 
  geom_pointrange(aes(ymin = log10(min_abund+1), ymax = log10(max_abund+1)), position = position_dodge(width = 0.7)) + 
  facet_wrap(~taq) + theme_bw() + 
  labs(x = "", y = expression(Log["10"]~median~counts)) + 
  scale_color_manual(name = "Higher 16S V4\nGC Content", values = c('#FF8000', '#8B5A2B')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        legend.position = c(0.85, 0.3), 
        legend.text = element_text(size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("results/figures/gc_content_graph.pdf", gc_graph, width = 6, height = 4.5, dpi = 300)






