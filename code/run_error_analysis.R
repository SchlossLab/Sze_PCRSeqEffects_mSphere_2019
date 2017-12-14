### Create analysis of the mock error and chimera
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "dunn.test"))

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


# Create ANOVA function to test differences between groups at each amplification cycle
run_comparison <- function(i, col_of_int, taq_name, dataList){
  # i is the sub sampled level of interest
  # dataList is the combined meta data with OTU count data table list
  
  # Add a vector with all possible cycles analyzed
  amp_cycles <- c("15x", "20x", "25x", "30x", "35x")
  # Pull out data table of interest
  tempData <- dataList
  #tempData <- dataList[[i]]
  # Runs the ANOVA across all the different amp cycles
  tempTests <- sapply(amp_cycles, 
                      function(x) run_kruskal(x, i, col_of_int, taq_name, tempData), simplify = F)
  # removes all the cycles that do not have any values
  tempTests <- tempTests[!sapply(tempTests, is.null)]
  # converts the list to a data frame, run the BH correction, and reorganizes the table
  tempTests <- tempTests %>% bind_rows() %>% 
    mutate(bh = p.adjust(pvalue, method = "BH")) %>% 
    select(statistic, parameter, pvalue, bh, cycle, sub_sample_level)
  # Return the results to the global work environment
  return(tempTests)
  
}


# Function that runs the actual ANOVA with a select data set
run_kruskal <- function(ac, subsample, comparator_var, taq_var, dataTable){
  # ac is the amplification cycle of interest
  # subsample is the sub sampled level of interest
  # dataTable is the data frame for a specific sub sampling level
  
  # create a new temp data frame filtered on the amplification cycle
  tempData <- dataTable %>% filter(cycles == ac)
  # run the ANOVA for that specific amplifcation cycle looking a Taq differences
  tempComparison <- as.data.frame.list(try(
    kruskal.test(formula(paste(comparator_var, " ~ ", "factor(", taq_var, ")", sep = "")), 
                data = tempData), silent = T))
  
  # Check to see if the length of the table is 1 due to error out in the try function
  if(length(colnames(tempComparison)) == 1){
    # assign a null place holder
    tempComparison <- c()
  } else{
    # rename the pvalue and add a cycle and sub-sample level to the data frame
    tempComparison <- tempComparison %>% 
      rename(pvalue = p.value) %>%  
      mutate(cycle = ac, 
             sub_sample_level = subsample)
  }
  # Return the data table to the work environment of the run_comparison function
  return(tempComparison)
}


# Function to run Tukey post hoc test on only those that are significant after bh correction
run_dunn <- function(i, col_of_int, taq_name, dataList, rawData){
  # i is the sub sample level of interest
  # dataList is the ANOVA list results
  # rawData is is the combined count and meta data list
  
  # Remove all samples that have a BH corrected value of over 0.05
  tempData <- dataList %>% filter(bh < 0.05)
  # tempData <- dataList[[i]] %>% 
  #   filter(bh < 0.05)
  # Pull out the specific cycle numbers that this occured for
  tempVector <- as.data.frame(tempData)[, "cycle"]
  # Check if the vector is empty or not
  if(length(rownames(tempData)) == 0){
    # Assigns the value a NULL place holder
    tempResults <- c()
  } else{
    # Runs the Tukey post hoc test for each amplification cycle in the tempVector
    tempResults <- sapply(tempVector, 
                          function(x) get_dunn_test(i, x, col_of_int, taq_name, rawData), simplify = F)
    # Converts the results from a list to a data frame
    tempResults <- tempResults %>% bind_rows()
    
  }
  # Returns the final table to the global work environment
  return(tempResults)
}

# Function that runs that actual Tukey post-hoc test
get_dunn_test <- function(subsample, cycle_num, comparator_var, taq_var, dataTable){
  # subsample is the sub sampled level of interest
  # cycle_num is the amplification cycle number of interest
  # dataTable is the combined count and meta data list
  
  # filter the data based on subsample level and cycle number of interest
  tempData <- as.data.frame(dataTable %>% filter(cycles == cycle_num))
  #tempData <- as.data.frame(dataTable[[subsample]] %>% filter(cycles == cycle_num))
  
  tempValues <- tempData[, comparator_var]
  tempLabels <- as.factor(tempData[, taq_var])
  
  # Run the Dunn test and convert to a data frame with cycle number and sub sampling level
  post_hoc_outcome <- as.data.frame.list(dunn.test(tempValues, tempLabels, method = "bh")) %>% 
    mutate(cycle = cycle_num, sub_sample_level = subsample) %>% 
    rename(pvalue = P, bh = P.adjusted)
  
  # Run the Tukey test and convert to a data frame with cycle number and sub sampling level
  #post_hoc_outcome <- TukeyHSD(aov(lm(
  #  formula(paste(comparator_var, " ~ ", taq_var, sep = "")), data = tempData)))[["taq"]] %>% 
  #  as.data.frame() %>% 
  #  mutate(comparison = rownames(.), cycle = cycle_num, sub_sample_level = subsample)
  # return results to the run_tukey work environment
  return(post_hoc_outcome)
  
}


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in subsample.shared files
# sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")
# 
# # Read in the count data
error_data <- read_csv("data/process/tables/full_error_summary.csv") %>% 
  filter(total_seqs >= 1000)

# error_data <- sapply(sub_sample_level, 
#                      function(x) read_data("data/process/tables/error_", "", "_summary.csv", x), 
#                      simplify = F)


# Run the kruskal comparisons between amp cycle across subsamplings error per base
krsukal_tests_error <- run_comparison("full_data", "mean_error", "taq", error_data)

#  kruskal_tests_error <- sapply(sub_sample_level, 
#                        function(x) run_comparison(x, "mean_error", "taq", error_data), simplify = F)
# 
# combined_kruskal_table_error <- kruskal_tests_error %>% bind_rows()

# Run the kruskal comparisons between amp cycle across subsamplings number of seqs with error
kruskal_tests_seq_error_count <- run_comparison("full_data", "seq_error_prevalence", "taq", error_data)

# kruskal_tests_seq_error_count <- sapply(sub_sample_level, 
#                             function(x) run_comparison(x, "seq_error_prevalence", "taq", error_data), simplify = F)
# 
# combined_kruskal_table_seq_error_count <- kruskal_tests_seq_error_count %>% bind_rows()


# Run the kruskal comparisons between amp cycle across subsamplings number of chimeras
kruskal_tests_chimera <- run_comparison("full_data", "chimera_prevalence", "taq", error_data)
# kruskal_tests_chimera <- sapply(sub_sample_level, 
#                                       function(x) run_comparison(x, "chimera_prevalence", "taq", error_data), simplify = F)
# 
# combined_kruskal_table_chimera <- kruskal_tests_chimera %>% bind_rows()


# Run the Dunn post-hoc test comparisons on only the Kruskal that were significant after BH correction
dunn_tests_error <- run_dunn("full_data", "mean_error", "taq", krsukal_tests_error, error_data)

# dunn_tests_error <- sapply(sub_sample_level, 
#                       function(x) run_dunn(x, "mean_error", "taq", 
#                                             kruskal_tests_error, error_data), simplify = F)
# 
# combined_dunn_error <- dunn_tests_error %>% bind_rows()

dunn_tests_seq_error_count <- run_dunn("full_data", "seq_error_prevalence", "taq", kruskal_tests_seq_error_count, error_data)
# dunn_tests_seq_error_count <- sapply(sub_sample_level, 
#                            function(x) run_dunn(x, "seq_error_prevalence", "taq", 
#                                                 kruskal_tests_seq_error_count, error_data), simplify = F)
# 
# combined_dunn_seq_error_count <- dunn_tests_seq_error_count %>% bind_rows()


dunn_tests_chimera <- run_dunn("full_data", "chimera_prevalence", "taq", kruskal_tests_chimera, error_data)
# dunn_tests_chimera <- sapply(sub_sample_level, 
#                                      function(x) run_dunn(x, "chimera_prevalence", "taq", 
#                                                           kruskal_tests_chimera, error_data), simplify = F)
# 
# combined_dunn_chimera<- dunn_tests_chimera %>% bind_rows()

# Add data table write out
write_csv(krsukal_tests_error, "data/process/tables/mock_error_overall_kruskal_results.csv")
write_csv(dunn_tests_error, "data/process/tables/mock_error_overall_dunn_results.csv")

write_csv(kruskal_tests_seq_error_count, "data/process/tables/mock_error_count_overall_kruskal_results.csv")
write_csv(dunn_tests_seq_error_count, "data/process/tables/mock_error_count_overall_dunn_results.csv")

write_csv(kruskal_tests_chimera, "data/process/tables/mock_chimera_overall_kruskal_results.csv")
write_csv(dunn_tests_chimera, "data/process/tables/mock_chimera_overall_dunn_results.csv")






