### Create graphs of the different numOTUs based on cycle number and sub sampling
### Fecal Sample Community only
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
                      function(x) read_data("data/process/tables/", "fecal_zscore_sub_sample_", "_count_table.csv", x) %>% 
                        mutate(depth_level = as.numeric(x)), 
                      simplify = F) %>% bind_rows() %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")), 
         cycle_num = factor(cycle_num, 
                            levels = c(15, 20, 25, 30, 35), 
                            labels = c("15x", "20x", "25x", "30x", "35x")))
  

summary_data <- numOTU_data %>% group_by(taq) %>% 
  summarise(taq_median = median(scaled_numOTU, na.rm = T))

fecal_samples <- numOTU_data %>% 
  ggplot(aes(depth_level, scaled_numOTU, color = cycle_num, group = taq)) + 
  geom_point(size = 2, alpha = 0.9, show.legend = T) + theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) + 
  geom_hline(aes(yintercept = taq_median), summary_data, color = "red", size = 1) + 
  facet_grid(. ~ taq) + 
  labs(x = "Sub-sampling Depth", y = "Z-Score Normalized Number of OTUs") + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#006400", "#00FF7F", "#63B8FF", "#104E8B")) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))
  

ggsave("results/figures/Figure1.pdf", fecal_samples, width = 10, height = 5, dpi = 300)
