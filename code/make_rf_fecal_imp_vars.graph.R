### Plot RF results MDA for classifying samples
# Marc Sze


## Load in needed functions 
source("code/functions.R")

## Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "scales", "viridis", "RColorBrewer", "ggrepel"))

# Load in the needed data
rf_data <- read_csv("data/process/tables/sample_group_rf_imp_vars_summary.csv") %>% 
  mutate(fill_data = ifelse(other_model == "neither", NA, other_model)) %>% 
  bind_rows(c(otu = "blank", genus = "blank")) %>% 
  mutate(median_mda = case_when(otu == "blank" ~ 100, 
                                TRUE ~ median_mda), 
         other_model = case_when(otu == "blank" ~ "neither", 
                   TRUE ~ other_model)) %>% 
  arrange(desc(median_mda))

pg_data <- read_csv("data/process/tables/polymerase_group_top10_rf_imp_vars_summary.csv")
cg_data <- read_csv("data/process/tables/cycle_group_top10_rf_imp_vars_summary.csv")

model_data <- read_csv("data/process/tables/sample_group_rf_model_summary.csv") %>% 
  mutate(model_type = "fecal_sample_group") %>% 
  bind_rows(read_csv("data/process/tables/polymerase_group_rf_model_summary.csv") %>% 
              mutate(model_type = "fecal_polymerase"), 
            read_csv("data/process/tables/cycles_group_rf_model_summary.csv") %>% 
              mutate(model_type = "fecal_cycles"), 
            read_csv("data/process/tables/polymerase_group_mock_rf_model_summary.csv") %>% 
              mutate(model_type = "mock_polymerase"), 
            read_csv("data/process/tables/cycles_group_mock_rf_model_summary.csv") %>% 
              mutate(model_type = "mock_cycles"))


graph_data <- model_data %>% 
  mutate(probs = exp(-logLoss)) %>% 
  select(logLoss, probs, model_type)

# Model comparisons
model_graph <- graph_data %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("fecal_sample_group", "fecal_polymerase", "fecal_cycles", 
                                        "mock_polymerase", "mock_cycles"), 
                             labels = c("Fecal Individual Model", "Fecal Polymerase Model", 
                                        "Fecal Cycles Model", "Mock Polymerase Model", "Mock Cycles Model"))) %>% 
  ggplot(aes(model_type, probs)) + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  geom_point(color = "darkgray") + theme_bw() + 
  labs(x = "", y = "Probablity of Correct Classification") + 
  coord_flip(ylim = c(0, 1)) + ggtitle("A") + 
  theme(plot.title = element_text(face="bold", hjust = -0.55, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




# Important Variables 
imp_taxa_graph <- rf_data %>% 
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
  geom_point(size = 3) + ggtitle("B") + coord_cartesian(ylim = c(0, 7)) + 
  theme_bw() + labs(x = "", y = "Median MDA") + 
  scale_shape_manual(name = "", values = c(1, 21, 21, 21)) + 
  scale_color_manual(name = "", values = c('#8B8386', '#B0171F', '#0000FF', '#308014')) + 
  scale_fill_manual(name = "", values = c('#8B8386', '#B0171F', '#0000FF', '#308014')) + 
  annotate("text", label = paste("Roseburia, ", "Cycles Model Median MDA = ", 
                                 as.character(format(filter(cg_data, otu == "Otu0006") %>% pull(median_mda), 
                                                                  digits = 2, nsmall = 2)), sep = ""), 
           x = 64, y = 4.27, size = 3.5, parse = F, color = '#B0171F') + 
  annotate("text", label = paste("Bacteroides, ", "Cycles Model Median MDA = ", 
                                 as.character(format(filter(cg_data, otu == "Otu0041") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 70, y = 3.9, size = 3.5, parse = F, color = '#B0171F') + 
  annotate("text", label = paste("Ruminococcus, ", "Polymerase Model Median MDA = ", 
                                 as.character(format(filter(pg_data, otu == "Otu0018") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 82, y = 3.31, size = 3.5, parse = F, color = '#0000FF') + 
  annotate("text", label = paste("Slackia, ", "Polymerase Model Median MDA = ", 
                                 as.character(format(filter(pg_data, otu == "Otu0047") %>% pull(median_mda), 
                                                     digits = 2, nsmall = 2)), sep = ""), 
           x = 86, y = 2.62, size = 3.5, parse = F, color = '#0000FF') + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



test <- grid.arrange(model_graph, imp_taxa_graph, layout_matrix = rbind(c(1, 1, 2, 2, 2)))

ggsave("results/figures/Figure5.pdf", test, width = 11, height = 5, dpi = 300)

             