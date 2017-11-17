### Create graphs of the different numOTUs based on cycle number and sub sampling
### Fecal Sample Community only
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
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = "")) %>% 
    filter(sample_name %in% c("DA10001", "DA10016", "DA10029", "DA10040"))
  # returns the file to the global environment
  return(tempData)
  
}


# Function to scale (z-score normalize with centering) the count data
normalize_data <- function(i, dataList){
  
  
  tempData <- dataList[[i]] %>% 
    group_by(sample_name, cycles) %>% 
    mutate(scaled_numOTU = scale(numOTUs)) %>% 
    filter(cycles != "15x") %>% 
    ungroup()
  
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
    mutate(bh = p.adjust(Pr..F., method = "BH")) %>% 
    select(Df, Sum.Sq, Mean.Sq, F.value, Pr..F., bh, cycle, sub_sample_level) %>% 
    rename(pvalue = Pr..F.)
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
    summary(aov(scaled_numOTU ~ taq, data = tempData)), silent = T))["taq", ] %>% tbl_df()
  # Check to see if the length of the table is 1 due to error out in the try function
  if(length(colnames(tempComparison)) == 1){
    # assign a null place holder
    tempComparison <- c()
  } else{
    # rename the pvalue and add a cycle and sub-sample level to the data frame
    tempComparison <- tempComparison %>% 
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

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000", "15000", "20000")

# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                      function(x) read_data("data/process/tables/", "fecal_sub_sample_", "_count_table.csv", x), 
                      simplify = F)

# Normalize the data
zscore_data <- sapply(sub_sample_level, 
                      function(x) normalize_data(x, numOTU_data), simplify = F)

# Run the ANOVA
anova_tests <- sapply(sub_sample_level, 
                      function(x) run_comparison(x, zscore_data), simplify = F)

combined_anova_table <- anova_tests %>% bind_rows()

# Run the Tukey post-hoc test comparisons on only the ANOVAs that were significant after BH correction
tukey_tests <- sapply(sub_sample_level, 
                      function(x) run_tukey(x, anova_tests, zscore_data), simplify = F)

combined_tukey_table <- tukey_tests %>% bind_rows() %>% arrange(`p adj`)


