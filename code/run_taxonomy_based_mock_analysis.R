### Run The taxonomy centered analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret"))


shared <- read_tsv("data/process/all_amp.0.03.subsample.1000.shared")

metadata <- shared %>% select(Group) %>% 
  separate(Group, c("x1", "x2", "x3", "x4", "x5")) %>% 
  mutate(sample_id = pull(shared, Group)) %>% 
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
  select(-x5)

taxonomy <- read_tsv("data/process/all_amp.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")


# Remve the nzv data
nzv <- nearZeroVar(as.data.frame(shared))

test_shared <- as.data.frame(shared)[, -nzv]

# Combine the otu and metadata together and tidy it up
combined_shared <- metadata %>% 
  inner_join(test_shared, by = c("sample_id" = "Group")) %>% 
  gather(otu, abund, contains("Otu"))


otu_test_data <- combined_shared %>% 
  group_by(cycles, otu) %>% 
  nest() %>% 
  mutate(kruskal_analysis = map(data, ~kruskal.test(abund ~ factor(taq), data = .x)), 
         summary_data = map(kruskal_analysis, broom::tidy)) %>% 
  select(cycles, otu, summary_data) %>% 
  unnest(summary_data) %>% 
  ungroup() %>% 
  filter(!is.na(statistic)) %>% 
  mutate(bh = p.adjust(p.value, method = "BH")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("otu" = "OTU")) %>% 
  select(cycles, otu, genus, statistic, p.value, bh)

  





