### Analyzing OTUS and diversity metrics from cycle number
### Analyze both actual samples and Mock communities

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr"))


# Read in Meta data file
metadata <- read_csv("data/process/tables/meta_data.csv")



###########################################################################################################################
############################### List of functions to make analysis work ###################################################
###########################################################################################################################

# Function to read in the data
read_data <- function(pathing, start_name, end_name, differentiator){
  
  tempData <- read_tsv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  
  return(tempData)
  
}


# Function to seperate data tables between DNA ctrls and FS
pare_down_table <- function(i, dataList, meta_table, remove_samples){
  
  tempMeta <- meta_table %>% filter(sample_type != remove_samples)
  
  tempData <- dataList[[i]] %>% 
    slice(match(Group, tempMeta$full_name))
    
  selected_tempData <- tempData %>% 
    select(-label, -Group, -numOtus) %>% 
    select(which(colSums(.) > 0)) %>% 
    mutate(Group = tempData$Group) %>% 
    select(Group, everything())
  
  return(selected_tempData)
}


# Function to generate OTU counts from subsample shared file
make_count_table <- function(i, dataList){
  
  # generate a present/abscence data table
  tempData <- dataList[[i]] %>% 
    select(-Group) %>% 
    mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1)))
  
  return(tempData)
  
}



# generate a present/abscence data table
dna_pres_absen_test <- dna_test %>% 
  select(-Group) %>% 
  mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1)))

sub_dna_pres_absen_test <- sub_dna_test %>% 
  select(-Group) %>% 
  mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1)))

stool_pres_absen_test <- stool_test %>% 
  select(-Group) %>% 
  mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1)))

sub_stool_pres_absen_test <- sub_stool_test %>% 
  select(-Group) %>% 
  mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1)))








###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################


# Read in subsample.shared files
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

sub_shared_data <- sapply(sub_sample_level, 
                          function(x) read_data("data/process/", "all_amp.0.03.subsample.", ".shared", x), 
                          simplify = F)

# Select only DNA MOCK samples and keep only OTUs that are not 0 for every sample
mock_samples_data <- sapply(sub_sample_level, 
               function(x) pare_down_table(x, sub_shared_data, metadata, "FS"), simplify = F)

# Generate OTU presence/absence table
mock_presence_table <- sapply(sub_sample_level, 
                              function(x) make_count_table(x, mock_samples_data), simplify = F)





