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


# Function to generate OTU counts from subsample shared file
make_count_table <- function(){
  
  
}




# Generate absence/presence groups
dna_test <- dna_shared %>% 
  select(-label, -Group, -numOtus) %>% 
  select(which(colSums(.) > 0)) %>% 
  mutate(Group = dna_shared$Group) %>% 
  select(Group, everything())

sub_dna_test <- sub_dna_shared %>% 
  select(-label, -Group, -numOtus) %>% 
  select(which(colSums(.) > 0)) %>% 
  mutate(Group = sub_dna_shared$Group) %>% 
  select(Group, everything())

stool_test <- stool_shared %>% 
  select(-label, -Group, -numOtus) %>% 
  select(which(colSums(.) > 0)) %>% 
  mutate(Group = stool_shared$Group) %>% 
  select(Group, everything())

sub_stool_test <- sub_stool_shared %>% 
  select(-label, -Group, -numOtus) %>% 
  select(which(colSums(.) > 0)) %>% 
  mutate(Group = sub_stool_shared$Group) %>% 
  select(Group, everything())


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










