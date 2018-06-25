### Run The taxonomy centered analysis
## For the fecal samples
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test"))


shared <- read_tsv("data/process/all_amp.0.03.subsample.1000.shared")

metadata <- shared %>% select(Group) %>% 
  separate(Group, c("x1", "x2", "x3", "x4", "x5")) %>% 
  mutate(sample_id = pull(shared, Group)) %>% 
  filter(x5 %in% c("DA10001", "DA10016", "DA10029", "DA10040")) %>% 
  rename(cycles = x1, taq = x2, sample_type = x3, sample_number = x4, sample_name = x5)


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


# Run the kruskal wallis test for every OTU and cycle number from 25x and above (20x and below will need custom comparison)
otu_test_data <- combined_shared %>% 
  filter(cycles != "15x",  cycles != "20x") %>% 
  group_by(cycles, otu) %>% 
  arrange(sample_name, .by_group = TRUE) %>% 
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
  # None are significant also going to check pvalue under 0.05


# Pull out only the kruskal wallis measures with pvalue < 0.05
pvalue_under_0.05_only <- otu_test_data %>% 
  filter(p.value < 0.05) %>% 
  arrange(cycles) %>% 
  unite(extra_name, cycles, otu, sep = "_")
  # There are only 2 significant OTUs and one is a mock community OTU


# Write out specific tables
write_csv(otu_test_data, "data/process/tables/kruskal_otu_polymerase_fecal_sample_diffs_summary.csv")













