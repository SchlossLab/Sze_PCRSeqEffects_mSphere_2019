### Create graphs of the different numOTUs based on cycle number and sub sampling
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis"))


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



# Generate graph of Mock DNA samples (not subsampled)
numOTU_data[["1000"]] %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_smooth(size = 1, method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("A") + coord_cartesian(ylim = c(0, 8)) + 
  annotate("text", label = paste("Sub-sampled to 1000 Sequences"), x = 1, y = 8.2, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.10, 0.18), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))








###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                          function(x) read_data("data/process/tables/", "mock_sub_sample_", "_count_table.csv", x), 
                          simplify = F)







