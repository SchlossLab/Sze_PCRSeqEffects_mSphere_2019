### Make nice master meta data file
### Create a master file to glue all other files together based on the files file
## Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr"))

# Read in Data tables needed
raw_file <- read_tsv("data/process/amp.files", col_names = F) %>% 
  rename(full_name = X1) %>% select(full_name)

separated_mock_file <- raw_file %>% 
  filter(str_detect(full_name, "DA1") == F) %>% 
  separate(full_name, c("cycles", "taq", "sample_type", "sample_number")) %>% 
  bind_cols(filter(raw_file, str_detect(full_name, "DA1") == F)[, "full_name"]) %>% 
  mutate(sample_number = str_replace(sample_number, "DNA", "")) %>% 
  select(full_name, cycles, taq, sample_type, sample_number) %>% 
  mutate(cycles = str_replace(cycles, "PMM", "30x"), 
         taq = str_replace(taq, "ZymoControl(\\d)", "ACC"))


