### Create custom files file for error analysis within mothur
## Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse"))

file_data <- read_tsv("data/process/stability_v2.files", col_names = F)

temp_data <- file_data %>% 
  filter(
    grepl("_ACC_", X1) == T | 
      grepl("_K_", X1) == T | 
      grepl("_PHU_", X1) == T | 
      grepl("_PL_", X1) == T | 
      grepl("_Q5_", X1) == T | 
      grepl("PMM_ZymoControl", X1) == T)

temp_data_no_samples <- temp_data %>% 
  filter(grepl("_DA10", X1) != T)


write_tsv(temp_data, "data/process/amp.files", col_names = F)

write_tsv(temp_data_no_samples, "data/process/mock_amp.files", col_names = F)
