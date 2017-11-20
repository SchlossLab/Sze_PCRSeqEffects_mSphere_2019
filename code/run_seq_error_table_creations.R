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

error_summary <- sapply(sub_sample_level, 
               function(x) read_data("data/process/", "mock_error.", ".summary", x), 
               simplify = F)

# Read in master meta data file
metadata <- read_csv("data/process/tables/meta_data.csv")

test_count <- count_tables[[4]]
test_error <- error_summary[[4]]

test <- test_error %>% 
  left_join(test_count, by = c("query" = "Representative_Sequence")) %>% 
  gather("sample_name", "total_seqs", colnames(test_count)[-c(1:2)]) %>% 
  left_join(metadata, by = c("sample_name" = "full_name"))








