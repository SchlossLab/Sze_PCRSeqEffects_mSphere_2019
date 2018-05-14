### Create graphs of the different numOTUs based on cycle number and sub sampling
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

# Vector of sub samples used
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                          function(x) 
                            read_data("data/process/tables/", 
                                      "mock_sub_sample_", "_count_table.csv", x) %>% 
                        mutate(depth_level = as.numeric(x)), 
                          simplify = F) %>% bind_rows() %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))


summary_data <- numOTU_data %>% 
  group_by(depth_level, taq, cycles) %>% 
  summarise(taq_mean = median(numOTUs, na.rm = T), 
            taq_min = min(numOTUs, na.rm = T), 
            taq_max = max(numOTUs, na.rm = T))
  

ref_otu_num <- read_tsv("data/process/zymo_mock_only.shared") %>% 
  select(numOtus)

# Generate graph of Mock DNA samples 
thousand_graph <- summary_data %>% filter(depth_level == 1000) %>% 
  ggplot(aes(cycles, taq_mean, color = taq, group = taq)) + 
  geom_pointrange(aes(ymin = taq_min, ymax = taq_max), size = 0.4, alpha = 0.9, show.legend = T, position = position_dodge(width = 0.8)) +
  geom_vline(xintercept=seq(1.5, length(unique(summary_data$depth_level))-0.5, 1), 
             lwd=1, colour="gray") + 
  geom_hline(yintercept = log2(ref_otu_num$numOtus), linetype = "dashed", size = 1) + 
  theme_bw() + 
  scale_color_manual(name = "HiFi DNA Polymerase", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Number of Cycles", y = "Number of OTUs") + 
  coord_cartesian(ylim = c(0, 70)) + 
  annotate("text", label = paste("Mock Data"), x = 0.8, y = 70, size = 3.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


ggsave("results/figures/Figure2.pdf", thousand_graph, width = 6, height = 5, dpi = 300)



