### Analyzing OTUS and diversity metrics from cycle number
### Analysis of Actual Samples communities only
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


# Function to seperate data tables between DNA ctrls and FS
pare_down_table <- function(i, dataList, meta_table, remove_samples1, remove_samples2){
  # i is the respective sub sample level (e.g. 50, 100)
  # dataList is the list where all the data is stored
  # meta_table is the master list of what each sample is
  # removes_samples is the variable to filter the data table by
  
  # Creates a subsetted meta data file based on samples that are to be removed
  tempMeta <- meta_table %>% filter(sample_type != remove_samples1, sample_type != remove_samples2)
  # matches the data table with the sub setted meta data file
  tempData <- dataList[[i]] %>% 
    slice(match(tempMeta$full_name, Group))
  # Remove label and numOtus variable and column sums are not 0  
  selected_tempData <- tempData %>% 
    select(-label, -Group, -numOtus) %>% 
    select(which(colSums(.) > 0)) %>% 
    mutate(Group = tempData$Group) %>% 
    select(Group, everything())
  # return the file to the global work environment
  return(selected_tempData)
}


# Function to generate OTU counts from subsample shared file
make_count_table <- function(i, dataList){
  # i is the sub sampled level of interest (e.g. 50, 100)
  # dataList is the pared down data table with only non-zero counts
  
  # select the specific data table
  tempData <- dataList[[i]]
  
  # generate a present/abscence data table
  count_tempData <- tempData %>% 
    select(-Group) %>% 
    mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1))) %>% 
    mutate(Group = tempData$Group) %>% 
    select(Group, everything())
  # Track the progress
  print(paste("Completed table ", i, " presence/absence analysis", sep = ""))
  # Return the table to the global work environment
  return(count_tempData)
  
}


# Function to sum counts and save them per sample
sum_OTU_counts <- function(i, dataList){
  # i is the sub sampled level of interest
  # dataList is the list with all the precense/absence counts
  
  # Pull the specific data table to interested
  tempData <- dataList[[i]]
  # Sums all the OTUs for a given sample
  sum_tempData <- tempData %>% 
    select(-Group) %>% 
    mutate(numOTUs = rowSums(.), 
           Group = tempData$Group) %>% 
    select(Group, numOTUs)
  # returns the final table to the global work environment  
  return(sum_tempData)
}


# Function to add final counts with meta data information
combine_tables <- function(i, dataList, meta_table){
  # i is the sub sampled level of interest
  # dataList is the list of OTU count data tables
  # meta_table is the selected meta data
  
  # combine the OTU count data with the meta data table
  tempData <- dataList[[i]] %>% 
    left_join(meta_table, by = c("Group" = "full_name"))
  # return the result to the global work environment
  return(tempData)
  
}



###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in subsample.shared files
sub_sample_level <- c("1000", "5000", "10000", "15000", "20000")

sub_shared_data <- sapply(sub_sample_level, 
                          function(x) read_data("data/process/", "all_amp.0.03.subsample.", ".shared", x), 
                          simplify = F)


# Select only DNA MOCK samples and keep only OTUs that are not 0 for every sample
mock_samples_data <- sapply(sub_sample_level, 
                            function(x) pare_down_table(x, sub_shared_data, metadata, "Mock", "Water"), simplify = F)


# Generate OTU presence/absence table
mock_presence_table <- sapply(sub_sample_level, 
                              function(x) make_count_table(x, mock_samples_data), simplify = F)

# Generate OTU counts for each sample 
mock_OTU_counts <- sapply(sub_sample_level, 
                          function(x) sum_OTU_counts(x, mock_presence_table), simplify = F)
# Combine counts with metadata
mock_OTU_combined_table <- sapply(sub_sample_level, 
                                  function(x) combine_tables(x, mock_OTU_counts, metadata), simplify = F)


# Write meta data with count tables for later use
sapply(sub_sample_level, 
       function(x) write_csv(mock_OTU_combined_table[[x]], 
                             paste("data/process/tables/fecal_sub_sample_", x, "_count_table.csv", sep = "")))







