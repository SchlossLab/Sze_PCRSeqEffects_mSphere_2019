### Run tests on GC content
## Does GC content and AT content change with cycle number
# Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test"))

# Load in needed data

error_summary <- read_tsv("data/process/mock_error.summary")

error_table <- read_csv("data/process/tables/non_summarised_nucleotide_mock_error_summary.csv") %>% 
  filter(str_detect(full_name, "30x_ACC_Mock") != T) %>%
  group_by(taq, cycles, sample_name.y, query) %>%
  summarise(total_seqs = length(query)) %>% 
  ungroup() %>% 
  filter(sample_name.y != "Water") %>% 
  mutate(cycles = case_when(
    cycles == "15x" ~ "fifteen", 
    cycles == "20x" ~ "twenty", 
    cycles == "25x" ~ "twenty_five", 
    cycles == "30x" ~ "thirty", 
    cycles == "35x" ~ "thirty_five", 
    TRUE ~ "uh oh"))


# Reduce down the number of columns and take only the most abundant from each taxa

# Staph aureus GC content 32.70% (PMID 22815443)
# Listeria monocytogenes GC content 38.04% (PMID 22689239)
# Escherichia coli GC content 50.6% (NCBI database)
# Salmonella enterica GC content 52.1% (PMID 22461552)
# Pseudomonas aeruginosa GC content 66.6% (PMID 17038190)
# Enterococcus faecalis GC content 38% (PMID 22366422)
# Lactobacillus fermentum GC content 49.7% (PMID 24072868)
# Bacillus subtilis GC content 43.61% (PMID 23105055)

wg_median_cutoff <- median(c(32.7, 38.04, 50.6, 52.1, 66.6, 38, 49.7, 43.61))


summarized_table <- error_summary %>% 
  group_by(reference) %>% 
  filter(weight == max(weight)) %>% 
  ungroup() %>% 
  select(query, reference, AA, TT, CC, GG, total) %>% 
  mutate(
    whole_genome_gc_percent = c(32.7, 38.04, 50.6, 52.1, 66.6, 38, 49.7, 43.61), 
    at_percent = (AA + TT)/total * 100, 
         gc_percent = (CC + GG)/total * 100, 
         subs_percent = 100 - (at_percent + gc_percent), 
         higher_gc = case_when(
           gc_percent > median(gc_percent) ~ "yes", 
           gc_percent <= median(gc_percent) ~ "no", 
           TRUE ~ "uh oh"), 
         whole_genome_higher_gc = case_when(
           whole_genome_gc_percent > wg_median_cutoff ~ "yes", 
           whole_genome_gc_percent <= wg_median_cutoff ~ "no", 
           TRUE ~ "uh oh"))



# Summarise the counts from each sequence by polymerase and cycle
seq_counts_table <- error_table %>% 
  group_by(query, taq, cycles) %>% 
  summarise(median_abund = median(total_seqs), 
            iqr25 = quantile(total_seqs)["25%"],
            iqr75 = quantile(total_seqs)["75%"],
            min_abund = min(total_seqs),
            max_abund = max(total_seqs))
  
test <- seq_counts_table %>% 
  inner_join(select(summarized_table, query, reference, higher_gc), by = "query") %>% 
  group_by(higher_gc, cycles, taq) %>% 
  summarise(median_of_median_abund = median(median_abund), 
            iqr25 = quantile(median_abund)["25%"],
            iqr75 = quantile(median_abund)["75%"],
            min_abund = min(median_abund),
            max_abund = max(median_abund))

whole_test <- seq_counts_table %>% 
  inner_join(select(summarized_table, query, reference, whole_genome_higher_gc), by = "query") %>% 
  group_by(whole_genome_higher_gc, cycles, taq) %>% 
  summarise(median_of_median_abund = median(median_abund), 
            iqr25 = quantile(median_abund)["25%"],
            iqr75 = quantile(median_abund)["75%"],
            min_abund = min(median_abund),
            max_abund = max(median_abund))

# Write out the table 
write_csv(test, "data/process/tables/gc_content_amp_summary.csv")

write_csv(whole_test, "data/process/tables/gc_content_whole_genome_amp_summary.csv")





