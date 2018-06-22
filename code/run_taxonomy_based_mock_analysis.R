### Run The taxonomy centered analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test"))


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

# Run the kruskal wallis test for every OTU and cycle number 
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

# Pull out only the significant kruskal wallis measures
sig_only <- otu_test_data %>% 
  filter(bh < 0.05) %>% 
  arrange(cycles) %>% 
  unite(extra_name, cycles, otu, sep = "_")

# Use a Dunn's post hoc test to find out which polymerase were different
post_hoc_otu_data <- combined_shared %>% 
  unite(extra_name, cycles, otu, sep = "_") %>% 
  filter(extra_name %in% pull(sig_only, extra_name)) %>% 
  group_by(extra_name) %>% 
  nest() %>% 
  mutate(dunn_analysis = map(data, ~with(data = .x, expr = t(as.data.frame.list(dunn.test(abund, factor(taq), method = "bh"))))), 
         summary_data = map(dunn_analysis, broom::tidy)) %>%  
  select(extra_name, summary_data) %>% 
  unnest(summary_data) %>% 
  rename(row_names = ".rownames")

# tidy up the results for easier summary table generation
test <- post_hoc_otu_data %>% 
  gather(group, measure_result, X1:X10) %>% 
  select(-group)

test_ids <- test %>% pull(extra_name)

# Create a final summary data table with only significant differences
summary_test_data <- data_frame(
  chi2 = filter(test, row_names == "chi2") %>% pull(measure_result), 
  Z = filter(test, row_names == "Z") %>% pull(measure_result), 
  pvalue = filter(test, row_names == "P") %>% pull(measure_result),
  bh = filter(test, row_names == "P.adjusted") %>% pull(measure_result), 
  comparison = filter(test, row_names == "comparisons") %>% pull(measure_result), 
  extra_name = test_ids[seq(1, length(test_ids), 5)]) %>% 
  filter(bh < 0.05) %>% 
  separate(extra_name, c("cycles", "otu"), sep = "_") %>% 
  left_join(select(taxonomy, OTU, genus), by = c("otu" = "OTU")) %>% 
  mutate(present_in_standard = case_when(
    genus %in% c("Pseudomonas", "Escherichia/Shigella", "Salmonella", "Lactobacillus", 
                 "Enterococcus", "Staphylococcus", "Listeria", "Bacillus") ~ "yes", 
    TRUE ~ "no"))

# Create a file with only the OTUs that are supposed to be in the mock
mock_members_summary <- summary_test_data %>% 
  filter(present_in_standard == "yes") %>% 
  arrange(cycles, otu)

# Write out specific tables
write_csv(otu_test_data, "data/process/tables/kruskal_otu_polymerase_mock_diffs_summary.csv")

write_csv(mock_members_summary, "data/process/tables/dunns_sig_otu_polymerase_mock_diffs_summary.csv")




