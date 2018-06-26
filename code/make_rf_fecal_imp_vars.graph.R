### Plot RF results MDA for classifying samples
# Marc Sze


## Load in needed functions 
source("code/functions.R")

## Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "scales", "viridis", "RColorBrewer", "ggrepel"))

# Load in the needed data
rf_data <- read_csv("data/process/tables/sample_group_rf_imp_vars_summary.csv") %>% 
  mutate(fill_data = ifelse(other_model == "neither", NA, other_model))

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
                              labels = c("Not in Top 10", "Cycles", "Polymerase", "Cycles & Polymerase")), 
         fill_data = factor(fill_data, 
                            levels = c("cycles", "polymerase", "both"), 
                            labels = c("Cycles", "Polymerase", "Cycles & Polymerase"))) %>% 
  ggplot(aes(otu, median_mda, color = other_model, shape = other_model, fill = other_model)) + 
  # geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3) + 
  theme_bw() + labs(x = "", y = "Median MDA") + 
  scale_shape_manual(name = "", values = c(1, 21, 21, 21)) + 
  scale_color_manual(name = "", values = c('#8B8386', '#B0171F', '#0000FF', '#308014')) + 
  scale_fill_manual(name = "", values = c('#8B8386', '#B0171F', '#0000FF', '#308014')) + 
  annotate("text", label = paste("Roseburia, ", "Cycles Model Median MDA = ", 
                                 as.character(format(filter(cg_data, otu == "Otu0006") %>% pull(median_mda), 
                                                                  digits = 2, nsmall = 2)), sep = ""), 
           x = 67, y = 4.27, size = 3.5, parse = F, color = '#B0171F') + 
  annotate("text", label = paste("Bacteroides, ", "Cycles Model Median MDA = ", 
                                 as.character(format(filter(cg_data, otu == "Otu0041") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 73, y = 3.9, size = 3.5, parse = F, color = '#B0171F') + 
  annotate("text", label = paste("Ruminococcus, ", "Polymerase Model Median MDA = ", 
                                 as.character(format(filter(pg_data, otu == "Otu0018") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 85, y = 3.31, size = 3.5, parse = F, color = '#0000FF') + 
  annotate("text", label = paste("Slackia, ", "Polymerase Model Median MDA = ", 
                                 as.character(format(filter(pg_data, otu == "Otu0047") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 88, y = 2.62, size = 3.5, parse = F, color = '#0000FF') + 
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



