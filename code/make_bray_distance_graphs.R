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
                                               x, "_data.csv", sep = "")) %>% 
                      mutate(depth_level = as.numeric(x)), simplify = F) %>% bind_rows() %>% 
  filter(taq != "K") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x vs 20x", "20x vs 25x", 
                                           "25x vs 30x", "30x vs 35x")))
  

mock_dist_data <- sapply(sub_sample_level, 
                    function(x) read_csv(paste("data/process/tables/mock_bray_5_cycle_dist_", 
                                               x, "_data.csv", sep = "")) %>% 
                      mutate(depth_level = as.numeric(x)), simplify = F) %>% bind_rows() %>% 
  filter(taq != "K") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x vs 20x", "20x vs 25x", 
                                           "25x vs 30x", "30x vs 35x")))

ann_text <- data.frame(cycle_compare = 2.2, distance = 0.95, lab = "Sub-Sampled", 
                       taq = factor("ACC", levels = c("ACC", "PHU", "PL", "Q5"), 
                                    labels = c("Accuprime", "Phusion", "Platinum", "Q5")))



fecal_graph <- dist_data %>% filter(!is.na(cycle_compare)) %>% 
  ggplot(aes(depth_level, distance, color = cycle_compare, group = taq)) + 
  geom_point(size = 2) + 
  facet_grid(~taq) + 
  scale_color_manual(name = "Cycle Comparison", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00")) + 
  scale_x_continuous(breaks = c(1000, 5000, 10000), 
                     labels = c("1000", "5000", "10000")) + 
  theme_bw() +  coord_cartesian(ylim = c(0, 1)) + ggtitle("A") + 
  labs(x = "Sub-sampling Depth", y = "Bray Curtis Index") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 70, hjust = 1), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


mock_graph <- mock_dist_data %>% filter(!is.na(cycle_compare)) %>% 
  ggplot(aes(depth_level, distance, color = cycle_compare, group = taq)) + 
  geom_point(size = 2) + 
  facet_grid(~taq) + 
  scale_color_manual(name = "Cycle Comparison", 
                     values = c("#FF00FF", "#0000FF", "#8B1A1A", "#EE9A00"), 
                     limits = c("15x vs 20x", "20x vs 25x", "25x vs 30x", "30x vs 35x")) + 
  scale_x_continuous(breaks = c(1000, 5000, 10000), 
                     labels = c("1000", "5000", "10000")) +
  theme_bw() +  coord_cartesian(ylim = c(0, 1)) + ggtitle("B") + 
  labs(x = "Sub-sampling Depth", y = "Bray Curtis Index") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 70, hjust = 1), 
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 7), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


combined_graph <- grid.arrange(fecal_graph, mock_graph, ncol = 2, nrow = 1)

ggsave("results/figures/Figure3.pdf", combined_graph, width = 8, height = 4, dpi = 300)
