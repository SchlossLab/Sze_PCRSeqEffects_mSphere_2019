### Plot RF results MDA for classifying samples
# Marc Sze


## Load in needed functions 
source("code/functions.R")

## Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "scales", "viridis", "RColorBrewer", "ggrepel"))

# Load in the needed data
rf_data <- read_csv("data/process/tables/sample_group_rf_imp_vars_summary.csv")

pg_data <- read_csv("data/process/tables/polymerase_group_top10_rf_imp_vars_summary.csv")
cg_data <- read_csv("data/process/tables/cycle_group_top10_rf_imp_vars_summary.csv")


rf_data %>% 
  mutate(otu = str_replace(otu, "Otu", ""), 
         otu = as.numeric(otu), 
         otu = paste("OTU", otu, sep = ""), 
         otu = factor(otu, 
                      levels = otu, 
                      labels = otu), 
         genus = factor(genus, 
                        levels = unique(genus), 
                        labels = unique(genus)), 
         other_model = factor(other_model, 
                              levels = c("neither", "cycles", "polymerase", "both"), 
                              labels = c("Not in Top 10", "Cycles", "Polymerase", "Cycles & Polymerase"))) %>% 
  ggplot(aes(otu, median_mda, color = other_model)) + 
  # geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 2) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c('#8B8386', '#B0171F', '#0000FF', '#308014')) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



  geom_text_repel(aes(label = otu), point.padding = 0.5, size = 3.0) + 
  scale_color_manual(name = "", values = getPalette(total), guide = guide_legend(nrow = 5)) + 
  labs(y = "Difference in MDA") + 
  annotate("text", label = 'atop(bold("More Important in Best Model"))', x = 18, y = 1.5, size = 3.5, parse = T) +
  annotate("text", label = 'atop(bold("More Important in Worse Model"))', x = 6, y = -0.75, size = 3.5, parse = T) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic", size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.)



