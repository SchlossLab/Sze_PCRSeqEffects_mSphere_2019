### Correlations between error and chimera with number of OTUs
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


# FUnction to create a new common column
make_common_column <- function(i, dataList){
  
  tempData <- dataList[[i]] %>% 
    mutate(common_column = paste(cycles, "_", taq, "_", sample_name, sep = ""))
  
  return(tempData)
  
}

# Function to join data together
combine_data <- function(i, dataList, to_match_list){
  
  tempData <- dataList[[i]] %>% 
    left_join(to_match_list[[i]], by = "common_column")
  
  return(tempData)
    
}


# Read in subsample.shared files
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in the count data
error_data <- sapply(sub_sample_level, 
                     function(x) read_data("data/process/tables/error_", "", "_summary.csv", x), 
                     simplify = F)


# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                      function(x) read_data("data/process/tables/", "mock_sub_sample_", "_count_table.csv", x), 
                      simplify = F)


# Create new common column
up_error_data <- sapply(sub_sample_level, 
                        function(x) make_common_column(x, error_data), simplify = F)

up_numOTU_data <- sapply(sub_sample_level, 
                         function(x) make_common_column(x, numOTU_data), simplify = F)

# Combine the data together
combined_list <- sapply(sub_sample_level, 
                        function(x) combine_data(x, up_error_data, up_numOTU_data), simplify = F)


combined_list[["1000"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_smooth(size = 1, method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("A") + coord_cartesian(ylim = c(0, 8)) + 
  annotate("text", label = paste("Sub-sampled to 1000 Sequences"), x = 1.5, y = 8.2, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.20, 0.18), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))



###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################





