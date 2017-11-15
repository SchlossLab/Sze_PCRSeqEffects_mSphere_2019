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


# Create ANOVA function to test differences between groups at each amplification cycle
run_comparison <- function(i, dataList){
  # i is the sub sampled level of interest
  # dataList is the combined meta data with OTU count data table list
  
  # Add a vector with all possible cycles analyzed
  amp_cycles <- c("15x", "20x", "25x", "30x", "35x")
  # Pull out data table of interest
  tempData <- dataList[[i]]
  # Runs the ANOVA across all the different amp cycles
  tempTests <- sapply(amp_cycles, 
                      function(x) run_anova(x, i, tempData), simplify = F)
  # removes all the cycles that do not have any values
  tempTests <- tempTests[!sapply(tempTests, is.null)]
  # converts the list to a data frame, run the BH correction, and reorganizes the table
  tempTests <- tempTests %>% bind_rows() %>% 
    mutate(bh = p.adjust(pvalue, method = "BH")) %>% 
    select(Df, Sum.Sq, Mean.Sq, F.value, pvalue, bh, cycle, sub_sample_level)
  # Return the results to the global work environment
  return(tempTests)
  
}


# Function that runs the actual ANOVA with a select data set
run_anova <- function(ac, subsample, dataTable){
  # ac is the amplification cycle of interest
  # subsample is the sub sampled level of interest
  # dataTable is the data frame for a specific sub sampling level
  
  # create a new temp data frame filtered on the amplification cycle
  tempData <- dataTable %>% filter(cycles == ac)
  # run the ANOVA for that specific amplifcation cycle looking a Taq differences
  tempComparison <- as.data.frame.list(try(
    summary(aov(numOTUs ~ taq, data = tempData)), silent = T))["taq", ] %>% tbl_df()
  # Check to see if the length of the table is 1 due to error out in the try function
  if(length(colnames(tempComparison)) == 1){
    # assign a null place holder
    tempComparison <- c()
  } else{
    # rename the pvalue and add a cycle and sub-sample level to the data frame
    tempComparison <- tempComparison %>% 
      rename(pvalue = Pr..F.) %>%  
      mutate(cycle = ac, 
             sub_sample_level = subsample)
  }
  # Return the data table to the work environment of the run_comparison function
  return(tempComparison)
}


# Function to run Tukey post hoc test on only those that are significant after bh correction
run_tukey <- function(i, dataList, rawData){
  # i is the sub sample level of interest
  # dataList is the ANOVA list results
  # rawData is is the combined count and meta data list
  
  # Remove all samples that have a BH corrected value of over 0.05
  tempData <- dataList[[i]] %>% 
    filter(bh < 0.05)
  # Pull out the specific cycle numbers that this occured for
  tempVector <- as.data.frame(tempData)[, "cycle"]
  # Check if the vector is empty or not
  if(length(rownames(tempData)) == 0){
    # Assigns the value a NULL place holder
    tempResults <- c()
  } else{
    # Runs the Tukey post hoc test for each amplification cycle in the tempVector
    tempResults <- sapply(tempVector, 
                          function(x) get_tukey_test(i, x, rawData), simplify = F)
    # Converts the results from a list to a data frame
    tempResults <- tempResults %>% bind_rows()
    
  }
  # Returns the final table to the global work environment
  return(tempResults)
}

# Function that runs that actual Tukey post-hoc test
get_tukey_test <- function(subsample, cycle_num, dataTable){
  # subsample is the sub sampled level of interest
  # cycle_num is the amplification cycle number of interest
  # dataTable is the combined count and meta data list
  
  # filter the data based on subsample level and cycle number of interest
  tempData <- dataTable[[subsample]] %>% filter(cycles == cycle_num)
  # Run the Tukey test and convert to a data frame with cycle number and sub sampling level
  post_hoc_outcome <- TukeyHSD(aov(lm(numOTUs ~ taq, data = tempData)))[["taq"]] %>% 
    as.data.frame() %>% 
    mutate(comparison = rownames(.), cycle = cycle_num, sub_sample_level = subsample)
  # return results to the run_tukey work environment
  return(post_hoc_outcome)
}
  
  
###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in subsample.shared files
sub_sample_level <- c("1000", "5000", "10000")

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


# Run the ANOVA comparisons between amp cycle across subsamplings
anova_tests <- sapply(sub_sample_level, 
                      function(x) run_comparison(x, mock_OTU_combined_table), simplify = F)

combined_anova_table <- anova_tests %>% bind_rows()

# Run the Tukey post-hoc test comparisons on only the ANOVAs that were significant after BH correction
tukey_tests <- sapply(sub_sample_level, 
                      function(x) run_tukey(x, anova_tests, mock_OTU_combined_table), simplify = F)

combined_tukey_table <- tukey_tests %>% bind_rows()


# Add data table write out
write_csv(combined_anova_table, "data/process/tables/fecal_overall_anova_results.csv")
write_csv(combined_tukey_table, "data/process/tables/fecal_overall_tukey_results.csv")

# Write meta data with count tables for later use
sapply(sub_sample_level, 
       function(x) write_csv(mock_OTU_combined_table[[x]], 
                             paste("data/process/tables/fecal_sub_sample_", x, "_count_table.csv", sep = "")))







