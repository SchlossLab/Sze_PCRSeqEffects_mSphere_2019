### Create graphs of the distances compared in 5 cycle increments (Bray)
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra"))

### Had to remove Kappa since it was missing 30x so could not do the same 5 increment comparison ###


###########################################################################################################################
############################### List of functions to make analysis work ###################################################
###########################################################################################################################

# Function to read in the data
read_data <- function(pathing, start_name, end_name, differentiator){
  # pathing is the variable that controls the path
  # start_name is the beginning of the file
  # end_name is how the file ends (e.g. .csv)
  # differentiator is the unique component of the file name
  
  # reads in the respective file 
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  # returns the file to the global environment
  return(tempData)
  
}


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000")

dist_data <- sapply(sub_sample_level, 
                    function(x) read_csv(paste("data/process/tables/bray_5_cycle_dist_", 
                                               x, "_data.csv", sep = "")), simplify = F)

mock_dist_data <- sapply(sub_sample_level, 
                    function(x) read_csv(paste("data/process/tables/mock_bray_5_cycle_dist_", 
                                               x, "_data.csv", sep = "")), simplify = F)

ann_text <- data.frame(cycle_compare = 2.2, distance = 0.95, lab = "Sub-Sampled", 
                       taq = factor("ACC", levels = c("ACC", "PHU", "PL", "Q5"), 
                                    labels = c("Accuprime", "Phusion", "Platinum", "Q5")))

fs_1000 <- dist_data[["1000"]] %>% 
  filter(taq != "K") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("A") +  
  geom_text(data = ann_text, color = "black", label = "Sub-Sampled to\n1000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))
  

fs_5000 <- dist_data[["5000"]] %>% 
  filter(taq != "K", !is.na(taq), cycle_compare != "15to25") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("B") +  
  geom_text(data = ann_text, color = "black", label = "Sub-Sampled to\n5000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


ann_text_1000 <- data.frame(cycle_compare = 1.8, distance = 0.95, lab = "Sub-Sampled", 
                       taq = factor("ACC", levels = c("ACC", "PHU", "PL", "Q5"), 
                                    labels = c("Accuprime", "Phusion", "Platinum", "Q5")))

fs_10000 <- dist_data[["10000"]] %>% 
  filter(taq != "K", !is.na(taq), cycle_compare != "15to25") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("C") +  
  geom_text(data = ann_text_1000, color = "black", label = "Sub-Sampled to\n10000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


ann_text <- data.frame(cycle_compare = 1.8, distance = 0.95, lab = "Sub-Sampled", 
                       taq = factor("ACC", levels = c("ACC", "PHU", "PL", "Q5"), 
                                    labels = c("Accuprime", "Phusion", "Platinum", "Q5")))

m_1000 <- mock_dist_data[["1000"]] %>% 
  filter(taq != "K") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("D") +  
  geom_text(data = ann_text, color = "black", label = "Sub-Sampled to\n1000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


m_5000 <- mock_dist_data[["5000"]] %>% 
  filter(taq != "K", !is.na(taq), cycle_compare != "15to25") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("E") +  
  geom_text(data = ann_text, color = "black", label = "Sub-Sampled to\n5000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


m_10000 <- mock_dist_data[["10000"]] %>% 
  filter(taq != "K", !is.na(taq), cycle_compare != "15to25", cycle_compare != "25to35") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x\nvs\n20x", "20x\nvs\n25x", "25x\nvs\n30x", "30x\nvs\n35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_point(size = 2) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  facet_grid(~taq) + 
  scale_color_manual(name = "Samples", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  labs(x = "", y = "Bray Curtis Index") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("F") +  
  geom_text(data = ann_text_1000, color = "black", label = "Sub-Sampled to\n10000 Sequences", size = 2) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))

combined_graph <- grid.arrange(fs_1000, fs_5000, fs_10000, m_1000, m_5000, m_10000, ncol = 3, nrow = 2)

ggsave("results/figures/Figure3.pdf", combined_graph, width = 15, height = 7, dpi = 300)
