### Create graphs of the distances compared in 5 cycle increments (Bray)
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra"))

### Had to remove Kappa since it was missing 30x so could not do the same 5 increment comparison ###

#Figure 3. the x-axis should be the number of cycles and the y-axis should be the BC distance between the current and previous number of cycles. 
# Again, present the median and IQR or range

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


mock_dist_data <- read_csv("data/process/tables/mock_bray_5_cycle_dist_1000_data.csv") %>%  
  filter(!is.na(distance), cycle_compare != "25to35") %>% 
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

fecal_graph <- dist_data %>% 
  ggplot(aes(cycle_compare, mean_dist, color = taq, group = taq)) + 
  geom_pointrange(aes(ymin = min_dist, ymax = max_dist), size = 0.4, show.legend = T, position = position_dodge(width = 0.8)) + 
  geom_vline(xintercept=seq(1.5, length(unique(dist_data$taq))-1.5, 1), 
             lwd=1, colour="gray") + 
  scale_color_manual(name = "HiFi DNA Polymerase", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  theme_bw() +  coord_cartesian(ylim = c(0, 1)) + ggtitle("A") + 
  labs(x = "Cycle Comparison", y = "Bray-Curtis Index") + 
  annotate("text", label = paste("Fecal Data"), x = 4.3, y = 1.0, size = 3.5) + 
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


mock_graph <- mock_dist_data %>% 
  ggplot(aes(cycle_compare, mean_dist, color = taq, group = taq)) + 
  geom_pointrange(aes(ymin = min_dist, ymax = max_dist), size = 0.4, show.legend = T, position = position_dodge(width = 0.8)) + 
  geom_vline(xintercept=seq(1.5, length(unique(mock_dist_data$taq))-2.5, 1), 
             lwd=1, colour="gray") + 
  scale_color_manual(name = "HiFi DNA Polymerase", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  theme_bw() +  coord_cartesian(ylim = c(0, 1)) + ggtitle("B") + 
  labs(x = "Cycle Comparison", y = "Bray-Curtis Index") + 
  annotate("text", label = paste("Mock Data"), x = 3.35, y = 1.0, size = 3.5) + 
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

combined_graph <- grid.arrange(fecal_graph, mock_graph, ncol = 1, nrow = 2)

ggsave("results/figures/Figure4.pdf", combined_graph, width = 6, height = 6.5, dpi = 300)
