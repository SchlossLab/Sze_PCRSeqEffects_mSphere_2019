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

# Function to create a full version and a slim version of the summary data files
combine_data <- function(countList, errorList, meta, simple = F){
  
  tempCount <- countList
  tempError <- errorList
  
  if(simple == F){
    
    tempData <- tempError %>% 
      left_join(tempCount, by = c("query" = "Representative_Sequence")) %>% 
      gather("sample_name", "total_seqs", colnames(tempCount)[-c(1:2)]) %>% 
      left_join(meta, by = c("sample_name" = "full_name")) %>% 
      filter(total_seqs != 0)
    
  } else{
    
    tempData <- tempError %>% 
      left_join(tempCount, by = c("query" = "Representative_Sequence")) %>% 
      gather("sample_name", "total_seqs", colnames(tempCount)[-c(1:2)]) %>% 
      left_join(meta, by = c("sample_name" = "full_name")) %>% 
      select(query, sample_name, total_seqs, cycles, taq, sample_type, sample_name.y, error, numparents) %>% 
      rename(full_name = sample_name, sample_name = sample_name.y) %>% 
      filter(total_seqs != 0)
  }
  
  return(tempData)
}





###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Create vector with different subsample levels
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in needed count files
count_table <- read_data("data/process/", "mock_error", ".count_table", "")

error_summary <- read_data("data/process/", "mock_error", ".summary", "")


# Read in master meta data file
metadata <- read_csv("data/process/tables/meta_data.csv")

# Combine and make full data tables
full_combined <- combine_data(count_table, error_summary, metadata)

# Combine and make slim data table
slim_combined <- combine_data(count_table, error_summary, metadata, simple = T)
  

##### I need to incorporate the number of times a sequence occurs





