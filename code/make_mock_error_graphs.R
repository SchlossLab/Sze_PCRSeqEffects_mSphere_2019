### Create graphs of the error based on cycle number and sub sampling
### Mock Community only
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
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  # returns the file to the global environment
  return(tempData)
  
}




###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

error_files <- c("mock_precluster_error", "mock_chimera_error", "mock_error") 

# # Read in the count data
error_data <- sapply(error_files, function(x) read_csv(paste("data/process/tables/", x, "_summary.csv", sep = "")) %>% 
                       filter(total_seqs >= 1000), simplify = F)


# Generate graph of Mock DNA samples 
full_pipeline <- error_data[["mock_error"]] %>% 
  filter(total_seqs >= 1000) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  group_by(taq, cycles) %>% 
  summarise(group_median = median(mean_error), group_iqr25 = quantile(mean_error)["25%"], 
            group_iqr75 = quantile(mean_error)["75%"]) %>% 
  ggplot(aes(cycles, group_median, color = taq, group = taq)) + 
  geom_line(show.legend = T) + 
  geom_errorbar(aes(ymin=group_iqr25, 
                    ymax=group_iqr75), show.legend = T, 
                width = 0.1, size = 0.5, alpha = 0.4) + 
  geom_point(size = 2, alpha = 0.7, show.legend = T) + 
  theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Number of Cycles", y = "Median Error Rate") + 
  coord_cartesian(ylim = c(0, 0.015)) + 
  annotate("text", label = paste("Mock Data"), x = 2.5, y = 0.015, size = 3.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.8, 0.8), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


 
 
#combined_graph <- grid.arrange(before_precluster, before_chimera_vsearch, full_pipeline, ncol = 3)

ggsave("results/figures/Figure4.pdf", full_pipeline, width = 4, height = 6, dpi = 300)
