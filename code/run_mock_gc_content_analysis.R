### Run tests on GC content
## Does GC content and AT content change with cycle number
# Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test"))

# Load in needed data

error_summary <- read_tsv("data/process/mock_error.summary")

error_table <- read_tsv("data/process/mock_error.count_table") %>% 
  gather(sample_id, abund, -Representative_Sequence, -total) %>% 
  separate(sample_id, c("x1", "x2", "x3", "x4", "x5")) %>% 
  filter(is.na(x5) | x3 == "Water") %>% 
  mutate(
    x4 = case_when(
      x2 == "A" ~ "DNA1", 
      x2 == "B" ~ "DNA2", 
      x2 == "C" ~ "DNA3", 
      x2 == "D" ~ "DNA4", 
      TRUE ~ x4), 
    x2 = case_when(
      x1 == "Zmock" ~ "ACC", 
      TRUE ~ x2),
    x3 = case_when(
      x1 == "Zmock" ~ "Mock", 
      TRUE ~ x3), 
    x1 = case_when(
      x1 == "Zmock" ~ "30x", 
      TRUE ~ x1)) %>% 
  rename(cycles = x1, taq = x2, sample_type = x3, sample_number = x4) %>% 
  select(-x5) %>% 
  filter(sample_type != "Water") %>% 
  mutate(cycles = case_when(
    cycles == "15x" ~ "fifteen", 
    cycles == "20x" ~ "twenty", 
    cycles == "25x" ~ "twenty_five", 
    cycles == "30x" ~ "thirty", 
    cycles == "35x" ~ "thirty_five", 
    TRUE ~ "uh oh"))


# Reduce down the number of columns and take only the most abundant from each taxa
summarized_table <- error_summary %>% 
  group_by(reference) %>% 
  filter(weight == max(weight)) %>% 
  ungroup() %>% 
  select(query, reference, AA, TT, CC, GG, total) %>% 
  mutate(at_percent = (AA + TT)/total * 100, 
         gc_percent = (CC + GG)/total * 100, 
         subs_percent = 100 - (at_percent + gc_percent), 
         higher_gc = case_when(
           gc_percent > median(gc_percent) ~ "yes", 
           gc_percent <= median(gc_percent) ~ "no", 
           TRUE ~ "uh oh"))


# Summarise the counts from each sequence by polymerase and cycle
seq_counts_table <- error_table %>% 
  group_by(Representative_Sequence, taq, cycles) %>% 
  summarise(median_abund = median(abund), 
            iqr25 = quantile(abund)["25%"],
            iqr75 = quantile(abund)["75%"],
            min_abund = min(abund),
            max_abund = max(abund))
  
test <- seq_counts_table %>% 
  inner_join(select(summarized_table, query, reference, higher_gc), by = c("Representative_Sequence" = "query")) %>% 
  group_by(higher_gc, cycles, taq) %>% 
  summarise(median_of_median_abund = median(median_abund), 
            iqr25 = quantile(median_abund)["25%"],
            iqr75 = quantile(median_abund)["75%"],
            min_abund = min(median_abund),
            max_abund = max(median_abund))

# Write out the table 
write_csv(test, "data/process/tables/gc_content_amp_summary.csv")







