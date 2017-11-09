### Pull Taq and PCR data results
### Create seperate data files for this component
## Marc Sze

#### Need to incorporate the 30x accuprime into the workflow since these will be named differently ####

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "viridis"))

# Read in Data tables needed
meta_data <- read_csv("data/process/final.metadata.csv")
shared <- read_tsv("data/process/final.shared")
susampled_shared <- read_tsv("data/process/final.0.03.subsample.shared")


# Keep only those that are for the PCR component
dna_meta_data <- meta_data %>% 
  filter(sample_type == "DNA", is.na(bead_beat) | bead_beat == "10min")

stool_meta_data <- meta_data %>% 
  slice(1:146) %>% 
  filter(sample_type == "Stool")

# Filter the shared file to only contain samples of interest
# Remove OTUs without any values in them
dna_shared <- shared %>% filter(Group %in% dna_meta_data$sample_id)
sub_dna_shared <- susampled_shared %>% filter(Group %in% dna_meta_data$sample_id)

stool_shared <- shared %>% filter(Group %in% stool_meta_data$sample_id)
sub_stool_shared <- susampled_shared %>% filter(Group %in% stool_meta_data$sample_id)

# Get subsampled amount
sub_amount <- unique(rowSums(sub_shared[, -c(1:3)]))

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


# Add the number of OTUs to the meta data
dna_meta_data <- dna_meta_data %>% 
  mutate(numOTUs = rowSums(dna_pres_absen_test),  
         sub_numOTUs = rowSums(sub_dna_pres_absen_test))

stool_meta_data <- stool_meta_data %>% 
  mutate(numOTUs = rowSums(stool_pres_absen_test),  
    sub_numOTUs = rowSums(sub_stool_pres_absen_test))


# Generate graph of Mock DNA samples (not subsampled)
dna_meta_data %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(amp_cycles, numOTUs, color = taq, group = taq)) + 
  geom_smooth(size = 1, se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = "Number of OTUs") + 
  ggtitle("Mock DNA")

  
# Generate graph of Mock DNA samples (subsampled)
dna_meta_data %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(amp_cycles, sub_numOTUs, color = taq, group = taq)) + 
  geom_smooth(size = 1, se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = "Number of OTUs") + 
  ggtitle("Mock DNA (Subsampled)")


# Generate graph of Actual DNA samples (not subsampled)
stool_meta_data %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(amp_cycles, numOTUs, color = taq, group = taq)) + 
  geom_smooth(size = 1, se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + coord_cartesian(ylim = c(0, 1000)) + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = "Number of OTUs") + 
  ggtitle("Actual Samples")



# Generate graph of Actual DNA samples (subsampled)
stool_meta_data %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(amp_cycles, sub_numOTUs, color = taq, group = taq)) + 
  geom_smooth(size = 1, se = FALSE) + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + coord_cartesian(ylim = c(0, 700)) + 
  scale_color_viridis(name = "Taq Used", discrete = TRUE) + 
  labs(x = "Amplification Cycles", y = "Number of OTUs") + 
  ggtitle("Actual Samples (Subsampled)")






