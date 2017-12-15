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
sub_sample_level <- c("1000", "5000", "10000", "15000", "20000")

dist_data <- sapply(sub_sample_level, 
                    function(x) read_csv(paste("data/process/tables/bray_5_cycle_dist_", 
                                               x, "_data.csv", sep = "")), simplify = F)



dist_data[["1000"]] %>% 
  filter(taq != "K") %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5")), 
         cycle_compare = factor(cycle_compare, 
                                levels = c("15to20", "20to25", "25to30", "30to35"), 
                                labels = c("15x to 20x", "20x to 25x", "25x to 30x", "30x to 35x"))) %>% 
  ggplot(aes(cycle_compare, distance, color = sample_name, group = taq)) + 
  geom_jitter(width = 0.2) + 
  facet_grid(~taq)
  








