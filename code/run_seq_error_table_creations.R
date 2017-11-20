### Code to read in and process seq error files
### transform to make it easier for R analysis
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr"))


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
  tempData <- read_tsv(paste(pathing, start_name, differentiator, end_name, sep = ""), na = character())
  # returns the file to the global environment
  return(tempData)
  
}







###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Create vector with different subsample levels
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in needed count files
count_tables <- sapply(sub_sample_level, 
                          function(x) read_data("data/process/", "mock_error.", ".count_table", x), 
                          simplify = F)

test <- read_tsv("data/process/mock_error.1000.summary", na = character())

# Read in seq data


