### Analyzing Percent Chimeric and % chimeric removed
### Analysis of Mock communities and fecal communities
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis"))


# Read in Meta data file
metadata <- read_csv("data/process/tables/meta_data.csv")



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
  tempData <- read_tsv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  # returns the file to the global environment
  return(tempData)
  
}


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in subsample.shared files
sub_sample_level <- c("1000", "5000", "10000")

sub_shared_data <- sapply(sub_sample_level, 
                          function(x) read_data("data/process/", "all_amp.0.03.subsample.", ".shared", x), 
                          simplify = F)
