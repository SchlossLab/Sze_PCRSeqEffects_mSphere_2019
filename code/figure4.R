### Create fecal specific community graph
### Fecal Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra"))

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
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = "")) %>% 
    filter(sample_name %in% c("DA10001", "DA10016", "DA10029", "DA10040")) %>% 
    mutate(cycle_num = as.numeric(str_replace(cycles, "x", "")))
  
  # returns the file to the global environment
  return(tempData)
  
}

###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000", "15000", "20000")

# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                      function(x) read_data("data/process/tables/", "fecal_sub_sample_", "_count_table.csv", x) %>% 
                        mutate(depth_level = as.numeric(x)), 
                      simplify = F) %>% bind_rows() %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")), 
         cycle_num = factor(cycle_num, 
                            levels = c(15, 20, 25, 30, 35), 
                            labels = c("15x", "20x", "25x", "30x", "35x")))


summary_data <- numOTU_data %>% group_by(depth_level, taq, cycle_num) %>% 
  summarise(taq_mean = median(numOTUs, na.rm = T), 
            taq_min = min(numOTUs, na.rm = T), 
            taq_max = max(numOTUs, na.rm = T))

overall_taq_mean_data <- numOTU_data %>% 
  filter(depth_level == 1000) %>% 
  group_by(taq) %>% 
  summarise(overall_mean = median(numOTUs, na.rm = T))

# Vector of sub samples used
dist_data <- read_csv("data/process/tables/bray_5_cycle_dist_1000_data.csv") %>% 
  filter(cycle_compare != "25to35") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")),
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x vs 20x", "20x vs 25x", 
                                           "25x vs 30x", "30x vs 35x"))) %>% 
  group_by(taq, cycle_compare) %>% 
  summarise(mean_dist = median(distance, na.rm = T), 
            max_dist = max(distance, na.rm = T), 
            min_dist = min(distance, na.rm = T))

fecal_samples <- summary_data %>% 
  filter(depth_level == 1000) %>% 
  ggplot(aes(cycle_num, taq_mean, color = taq, group = taq)) + 
  #  geom_hline(aes(yintercept = overall_mean, color = taq), overall_taq_mean_data, size = 1, 
  #             alpha = 0.5, linetype = "dashed", show.legend = F) + 
  geom_pointrange(aes(ymin = taq_min, ymax = taq_max), size = 0.4, alpha = 0.9, show.legend = T, position = position_dodge(width = 0.8)) + 
  geom_vline(xintercept=seq(1.5, length(unique(summary_data$depth_level))-0.5, 1), 
             lwd=1, colour="gray") + 
  theme_bw() + coord_cartesian(ylim = c(0, 150)) + ggtitle("A") + 
  scale_color_manual(name = "HiFi DNA Polymerase", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Number of Cycles", y = "Number of OTUs") + 
  annotate("text", label = paste("Fecal Data"), x = 5.1, y = 150, size = 3.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.position = "bottom", 
        legend.background = element_rect(color = "black"))


fecal_graph <- dist_data %>% 
  ggplot(aes(cycle_compare, mean_dist, color = taq, group = taq)) + 
  geom_pointrange(aes(ymin = min_dist, ymax = max_dist), size = 0.4, show.legend = T, position = position_dodge(width = 0.8)) + 
  geom_vline(xintercept=seq(1.5, length(unique(dist_data$taq))-1.5, 1), 
             lwd=1, colour="gray") + 
  scale_color_manual(name = "HiFi DNA Polymerase", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  theme_bw() +  coord_cartesian(ylim = c(0, 1)) + ggtitle("B") + 
  labs(x = "Cycle Comparison", y = "Bray-Curtis Index") + 
  annotate("text", label = paste("Fecal Data"), x = 4.2, y = 1.0, size = 3.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))

combined_graph <- grid.arrange(fecal_samples, fecal_graph, layout_matrix = cbind(c(1), c(2)))

ggsave("results/figures/Figure4.pdf", combined_graph, width = 11, height = 5, dpi = 300)
