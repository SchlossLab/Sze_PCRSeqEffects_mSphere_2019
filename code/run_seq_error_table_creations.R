### Code to read in and process seq error files
### transform to make it easier for R analysis
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr"))


### Recognized an error in the error calculation and will need to go back and run through mothur to get correct errors ###


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


# Function to generate sampling vectors
get_rand_vectors <- function(i, dataList){
  
  test_set <- dataList[[i]]$total_seqs
  
  total_unique <- length(test_set)
  
  test_vector <- rep(1:total_unique, test_set)
  
  return(test_vector)
  
}

# Function to check if total sequence is equal to or greater than depth and toss those that arent
get_depth_check <- function(i, depth_check, dataList){
  
  tempData <- dataList[[i]]
  
  ifelse(length(tempData) >= depth_check, 
         return(tempData), return(NULL))
  
}

# Function that randomly samples a set amount
run_random_sampling <- function(depth_used, tempData){
  
  set.seed(12345)
  temp_sampling <- sample(tempData, size = depth_used, replace = F)
  
  return(temp_sampling)
}

# Function to pull full seq information based on the random sampling
run_recombine_data <- function(i, RsamplingList, fullDataList){
  
  # Pull specific samples from full data list to make a new data frame
  tempData <- fullDataList[[i]][RsamplingList[[i]], ] %>% select(-total_seqs) %>% 
    mutate(sequence_number = rownames(.))
  
  return(tempData)
}



# Function to randomly sample each sample to a specified depth
get_random_sample <- function(depth, dataTable, sample_vector){

  tempList <- sapply(sample_vector, 
                     function(x) filter(dataTable, full_name == x), simplify = F)
  
  tempSamplingList <- sapply(sample_vector, 
                             function(x) get_rand_vectors(x, tempList))
  
  above_depth_List <- sapply(sample_vector, 
                             function(x) get_depth_check(x, depth, tempSamplingList), simplify = F)
  
  # Remove null elements
  above_depth_List <- above_depth_List[!sapply(above_depth_List, is.null)] 
  
  # Runs the random sampling across samples
  tempRSamplingList <- lapply(above_depth_List, 
                              function(x) run_random_sampling(depth, x))
  
  # Replace the current sample vectors with only those included in the random sampling
  sample_vector <- names(tempRSamplingList)
  
  # Recombine with full data set
  tempRecombinedList <- sapply(sample_vector, 
                               function(x) run_recombine_data(x, tempRSamplingList, tempList), simplify = F)
  
  return(tempRecombinedList)
}


# Function that checks for whether sequence is chimeric or not based on number of parents and sequence errors
run_checks <- function(dataTable){
  
  tempData <- dataTable %>% 
    mutate(chimera = ifelse(numparents == 1, invisible(0), invisible(1)), 
           seq_error = ifelse(error > 0, invisible(1), invisible(0)))
  
  return(tempData)
}


# Function to create summary data
get_summary_data <- function(dataList, depth){
  
  dataList <- dataList[[as.character(depth)]]
  
  tempData <- dataList %>% 
    bind_rows() %>% 
    group_by(taq, cycles, sample_name) %>% 
    summarise(mean_error = mean(error, na.rm = T), sd_error = sd(error, na.rm = T), 
              chimera_prevalence = sum(chimera)/depth, 
              seq_error_prevalence = sum(seq_error)/depth)
  
  return(tempData)
}


# Function to count the number of sequences type of error
run_nucleotide_error_check <- function(dataTable){
  
  tempData <- dataTable %>% 
    mutate(at = ifelse(AT == 0, invisible(0), invisible(1)), 
           ag = ifelse(AG == 0, invisible(0), invisible(1)), 
           ac = ifelse(AC == 0, invisible(0), invisible(1)), 
           ta = ifelse(TA == 0, invisible(0), invisible(1)), 
           tg = ifelse(TG == 0, invisible(0), invisible(1)), 
           tc = ifelse(TC == 0, invisible(0), invisible(1)), 
           ga = ifelse(GA == 0, invisible(0), invisible(1)), 
           gt = ifelse(GT == 0, invisible(0), invisible(1)), 
           gc = ifelse(GC == 0, invisible(0), invisible(1)), 
           ca = ifelse(CA == 0, invisible(0), invisible(1)), 
           ct = ifelse(CT == 0, invisible(0), invisible(1)), 
           cg = ifelse(CG == 0, invisible(0), invisible(1)))
  
  return(tempData)
}


# Function to create individual nucleotide summary data
get_nucleotide_summary_data <- function(dataList, depth){
  
  dataList <- dataList[[as.character(depth)]]
  
  tempData <- dataList %>% 
    bind_rows() %>% 
    group_by(taq, cycles, sample_name.y) %>% 
    summarise(at_rate = sum(AT)/sum(total.x), ag_rate = sum(AG)/sum(total.x), 
              ac_rate = sum(AC)/sum(total.x), ta_rate = sum(TA)/sum(total.x), 
              tg_rate = sum(TG)/sum(total.x), tc_rate = sum(TC)/sum(total.x), 
              ga_rate = sum(GA)/sum(total.x), gt_rate = sum(GT)/sum(total.x), 
              gc_rate = sum(GC)/sum(total.x), ca_rate = sum(CA)/sum(total.x), 
              ct_rate = sum(CT)/sum(total.x), cg_rate = sum(CG)/sum(total.x))
      
      
      
      #mean_error = mean(error, na.rm = T), sd_error = sd(error, na.rm = T), 
      #        chimera_prevalence = sum(chimera)/depth, 
       #       seq_error_prevalence = sum(seq_error)/depth)
  
  return(tempData)
}




###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Create vector with different subsample levels
sub_sample_level <- c(50, 100, 500, 1000, 5000, 10000)
#taqs_used <- c("acc", "k", "phu", "pl", "q5")


# Read in needed count files
count_table <- read_data("data/process/", "mock_error", ".count_table", "")

#good_counts <- sapply(taqs_used, 
#                      function(x) read_data("data/process/", x, ".count_table", "_data_only"), simplify = F)

error_summary <- read_data("data/process/", "mock_error", ".summary", "")

#good_error <- sapply(taqs_used, 
#                      function(x) read_data("data/process/", x, ".summary", "_data_only"), simplify = F)

# Read in master meta data file
metadata <- read_csv("data/process/tables/meta_data.csv")

# Combine and make full data tables
full_combined <- combine_data(count_table, error_summary, metadata) %>% 
  rename(full_name = sample_name)

# Combine and make slim data table
slim_combined <- combine_data(count_table, error_summary, metadata, simple = T)
  
# Generate the random sampling 
testRecombine <- sapply(sub_sample_level, 
                        function(x) get_random_sample(x, slim_combined, unique(slim_combined$full_name)), simplify = F)

names(testRecombine) <- sub_sample_level  


# Add Chimera checked column
testRecombine <- lapply(testRecombine, 
                        function(x) lapply(x, function(y) run_checks(y)))

# Generate Summary Data
good_summary_data <- sapply(sub_sample_level, 
                            function(x) get_summary_data(testRecombine, x), simplify = F)

# Write out summarized data tables for graphing
sapply(c(1:length(good_summary_data)), 
       function(x) write_csv(good_summary_data[[x]], 
                             paste("data/process/tables/error_", sub_sample_level[x], "_summary.csv", sep = "")))


###########################################################################################################################
############################### Run actual analysis for actual seq error ##################################################
###########################################################################################################################

# Generate the random sampling 
full_Recombine <- sapply(sub_sample_level, 
                        function(x) get_random_sample(x, full_combined, unique(full_combined$full_name)), simplify = F)

names(full_Recombine) <- sub_sample_level  

# Add nucleotide changes checked column
full_Recombine <- lapply(full_Recombine, 
                        function(x) lapply(x, function(y) run_nucleotide_error_check(y)))

# Generate Nucleotide Summary Data
nucleotide_summary_data <- sapply(sub_sample_level, 
                            function(x) get_nucleotide_summary_data(full_Recombine, x), simplify = F)







