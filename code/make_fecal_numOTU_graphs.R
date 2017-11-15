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
    filter(sample_name %in% c("DA10001", "DA10016", "DA10029", "DA10040"))
  # returns the file to the global environment
  return(tempData)
  
}

###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000")

# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                      function(x) read_data("data/process/tables/", "fecal_sub_sample_", "_count_table.csv", x), 
                      simplify = F)


thousand_graph <- numOTU_data[["1000"]] %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_point(size = 2, alpha = 0.7, show.legend = F) + theme_bw() + 
  facet_grid(. ~ sample_name) + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("A") + coord_cartesian(ylim = c(0, 9)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.08, 0.3), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


five_thousand_graph <- numOTU_data[["5000"]] %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_point(size = 2, alpha = 0.7, show.legend = F) + theme_bw() + 
  facet_grid(. ~ sample_name) + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("B") + coord_cartesian(ylim = c(0, 9)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.08, 0.3), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


ten_thousand_graph <- numOTU_data[["10000"]] %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(. ~ sample_name) + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("C") + coord_cartesian(ylim = c(0, 9)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.08, 0.3), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


combined_graph <- grid.arrange(thousand_graph, five_thousand_graph, ten_thousand_graph, nrow = 3)

ggsave("results/figures/fecal_numOTU_graph.pdf", combined_graph, width = 8, height = 11, dpi = 300)
