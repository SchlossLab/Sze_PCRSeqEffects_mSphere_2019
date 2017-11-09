### Make nice master meta data file
### Create a master files based on information in the used files file
## Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr"))

# Read in Data tables needed
raw_file <- read_tsv("data/process/amp.files", col_names = F) %>% 
  rename(full_name = X1) %>% select(full_name)

separated_mock_file <- raw_file %>% 
  filter(str_detect(full_name, "DA1") == F) %>% 
  separate(full_name, c("cycles", "taq", "sample_type", "sample_name")) %>% 
  bind_cols(filter(raw_file, str_detect(full_name, "DA1") == F)[, "full_name"]) %>% 
  select(full_name, cycles, taq, sample_type, sample_name) %>% 
  mutate(cycles = str_replace(cycles, "PMM", "30x"), 
         taq = str_replace(taq, "ZymoControl(\\d)", "ACC"))


separated_samples_file <- raw_file %>% 
  filter(str_detect(full_name, "DA1") == T) %>% 
  separate(full_name, c("cycles", "taq", "ext_kit", "sample_type", "sample_name")) %>% 
  bind_cols(filter(raw_file, str_detect(full_name, "DA1") == T)[, "full_name"]) %>% 
  select(full_name, cycles, taq, sample_type, sample_name)


combined_file <- separated_mock_file %>% bind_rows(separated_samples_file)

write_csv(combined_file, "data/process/tables/meta_data.csv")



