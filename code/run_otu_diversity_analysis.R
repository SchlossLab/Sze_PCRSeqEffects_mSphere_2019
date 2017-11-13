### Analyzing OTUS and diversity metrics from cycle number
### Analyze both actual samples and Mock communities

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
  
  tempData <- read_tsv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  
  return(tempData)
  
}


# Function to seperate data tables between DNA ctrls and FS
pare_down_table <- function(i, dataList, meta_table, remove_samples){
  
  tempMeta <- meta_table %>% filter(sample_type != remove_samples)
  
  tempData <- dataList[[i]] %>% 
    slice(match(tempMeta$full_name, Group))
    
  selected_tempData <- tempData %>% 
    select(-label, -Group, -numOtus) %>% 
    select(which(colSums(.) > 0)) %>% 
    mutate(Group = tempData$Group) %>% 
    select(Group, everything())
  
  return(selected_tempData)
}


# Function to generate OTU counts from subsample shared file
make_count_table <- function(i, dataList){
  
  tempData <- dataList[[i]]
  
  # generate a present/abscence data table
  count_tempData <- tempData %>% 
    select(-Group) %>% 
    mutate_all(function(x) ifelse(x == 0, invisible(0), invisible(1))) %>% 
    mutate(Group = tempData$Group) %>% 
    select(Group, everything())
  
  print(paste("Completed table ", i, " presence/absence analysis", sep = ""))
  
  return(count_tempData)
  
}


# Function to sum counts and save them per sample
sum_OTU_counts <- function(i, dataList){
  
  tempData <- dataList[[i]]
  
  sum_tempData <- tempData %>% 
    select(-Group) %>% 
    mutate(numOTUs = rowSums(.), 
           Group = tempData$Group) %>% 
    select(Group, numOTUs)
    
  return(sum_tempData)
}


# Function to add final counts with meta data information
combine_tables <- function(i, dataList, meta_table){
  
  tempData <- dataList[[i]] %>% 
    left_join(meta_table, by = c("Group" = "full_name"))
  
  return(tempData)
  
}




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

# Generate OTU counts for each sample 
mock_OTU_counts <- sapply(sub_sample_level, 
                          function(x) sum_OTU_counts(x, mock_presence_table), simplify = F)

# Combine counts with metadata
mock_OTU_combined_table <- sapply(sub_sample_level, 
                                  function(x) combine_tables(x, mock_OTU_counts, metadata), simplify = F)


# Generate graph of Mock DNA samples (not subsampled)
mock_OTU_combined_table[[4]] %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(cycles, log2(numOTUs), color = taq, group = taq)) + 
  geom_smooth(size = 1, method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = expression(Log["2"]~Number~of~OTUs)) + 
  ggtitle("Mock DNA") + coord_cartesian(ylim = c(0, 8))










