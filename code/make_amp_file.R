### Create custom files file for error analysis within mothur
## Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse"))

# Load in needed data tables
file_data <- read_tsv("data/process/stability_v2.files", col_names = F)

seq_metadata <- read_tsv("data/process/seq_metadata.txt") %>% 
  select(Run, Sample_Name)

temp_data <- file_data %>% 
  filter(
    grepl("_ACC_", X1) == T | 
      grepl("_K_", X1) == T | 
      grepl("_PHU_", X1) == T | 
      grepl("_PL_", X1) == T | 
      grepl("_Q5_", X1) == T | 
      grepl("PMM_ZymoControl", X1) == T | 
      grepl("Zmock", X1) == T) %>% 
  left_join(seq_metadata, by = c("X1" = "Sample_Name")) %>% 
  mutate(R1 = paste("data/process/", Run, "_1.fastq", sep = ""), 
         R2 = paste("data/process/", Run, "_2.fastq", sep = "")) %>% 
  select(X1, R1, R2)

temp_data_no_samples <- temp_data %>% 
  filter(grepl("_DA10", X1) != T) %>% 
  left_join(seq_metadata, by = c("X1" = "Sample_Name")) %>% 
  mutate(R1 = paste("data/process/", Run, "_1.fastq", sep = ""), 
         R2 = paste("data/process/", Run, "_2.fastq", sep = "")) %>% 
  select(X1, R1, R2)


# Create files file for each Taq specifically
seperated_data <- list(
  acc_data = temp_data %>% 
    filter(grepl("_ACC_", X1) == T |
           grepl("Zmock", X1) == T) %>% 
    filter(grepl("_Water", X1) == F, grepl("_DA10", X1) != T), 
  k_data = temp_data %>% 
    filter(grepl("_K_", X1) == T) %>% 
    filter(grepl("_Water", X1) == F, grepl("_DA10", X1) != T), 
  phu_data = temp_data %>% 
    filter(grepl("_PHU_", X1) == T) %>% 
    filter(grepl("_Water", X1) == F, grepl("_DA10", X1) != T), 
  pl_data = temp_data %>% 
    filter(grepl("_PL_", X1) == T) %>% 
    filter(grepl("_Water", X1) == F, grepl("_DA10", X1) != T), 
  q5_data = temp_data %>% 
    filter(grepl("_Q5_", X1) == T) %>% 
    filter(grepl("_Water", X1) == F, grepl("_DA10", X1) != T))


# Write out needed files file for mothur processing
write_tsv(temp_data, "data/process/amp.files", col_names = F)

write_tsv(temp_data_no_samples, "data/process/mock_amp.files", col_names = F)

sapply(names(seperated_data), 
       function(x) write_tsv(seperated_data[[x]], 
                             paste("data/process/", x, "_only.files", sep = ""), col_names = F))
