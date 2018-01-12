### Analyzing Percent Chimeric and % chimeric removed
### Analysis of Mock communities and fecal communities
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


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in subsample.shared files
sub_sample_level <- c("1000", "5000", "10000")
error_files <- c("mock_precluster_error", "mock_chimera_error", "mock_error") 

mock_count_tables <- sapply(error_files, 
                          function(x) read_data("data/process/tables/", "mock_sub_sample_", "_count_table.csv", x), 
                          simplify = F)

fecal_count_tables <- sapply(sub_sample_level, 
                             function(x) read_data("data/process/tables/", "fecal_sub_sample_", "_count_table.csv", x), 
                             simplify = F)
  
# Read in chimera data tables and respective count table

combined_table <- read_tsv("data/process/before_chimera_vsearch.count_table") %>% 
  left_join(
    read_tsv("data/process/all_amp.chimeras", col_names = F) %>% 
      select(X1, X2, X18) %>% rename(sequence_name = X2, chimera = X18) %>% 
      separate(sequence_name, c("sequence_name", "extra"), sep = ";") %>% 
      select(X1, sequence_name, chimera), 
    by = c("Representative_Sequence" = "sequence_name"))

# Rearrange the table and generate summary stats on chimera number versus totals
rearranged_combined_table <- combined_table %>% gather(sample_name, seq_count, 3:281) %>% 
  mutate(chimera_positive = ifelse(chimera == "N", invisible(0), invisible(1))) %>% 
  group_by(sample_name, chimera) %>% 
  summarise(total_seqs = sum(seq_count))


