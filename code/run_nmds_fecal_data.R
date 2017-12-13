### Create graphs of the different numOTUs based on cycle number and sub sampling
### Fecal Sample Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "vegan"))


###########################################################################################################################
############################### List of functions to make analysis work ###################################################
###########################################################################################################################


# Function to read in the data
read_data <- function(pathing, start_name, end_name, differentiator){
  # pathing is the variable that controls the path
  # start_name is the beginning of the file
  # end_name is how the file ends (e.g. .csv)
  # differentiator is the unique component of the file name
  
  # reads in the respective file 
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = "")) %>% 
    filter(sample_name %in% c("DA10001", "DA10016", "DA10029", "DA10040"))
  # returns the file to the global environment
  return(tempData)
  
}


















###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000", "15000", "20000")


# Read in distance matrix
braycurtis_dist <- read.dist('data/process/all_amp.braycurtis.0.03.lt.1000.dist')


# Read in meta data
metadata <- read_data("data/process/tables/", "meta_data", ".csv", "") %>% 
  filter(full_name %in% rownames(braycurtis_dist))

sample_to_keep <- metadata$full_name

braycurtis_dist <- braycurtis_dist[sample_to_keep, sample_to_keep]

# create select Taq distance table
phu <- filter(metadata, taq == "PHU")[, "full_name"]
q5 <- filter(metadata, taq == "Q5")[, "full_name"]

cycles <- filter(metadata, taq == "PHU")[, "cycles"]
q5_cycles <- filter(metadata, taq == "Q5")[, "cycles"]

phu_braycurtis <- braycurtis_dist[phu$full_name, phu$full_name]
q5_braycurtis <- braycurtis_dist[q5$full_name, q5$full_name]

# get PERMANOVA
set.seed(12345)
tempAnalysisData <- 
  adonis(as.dist(phu_braycurtis) ~ factor(cycles$cycles))

set.seed(12345)
q5_AnalysisData <- 
  adonis(as.dist(q5_braycurtis) ~ factor(q5_cycles$cycles))

# store important values
beta_diver_summary <- c(tempAnalysisData$aov.tab$F.Model[1], 
                             tempAnalysisData$aov.tab$R2[1], tempAnalysisData$aov.tab$`Pr(>F)`[1])

# Create NMDS values
set.seed(12345)
phu_bray.mds.data <- metaMDS(as.dist(phu_braycurtis), trymax = 3000, trace = 0) %>% 
  scores() %>% as.data.frame() %>% mutate(cycles = factor(cycles$cycles))

set.seed(12345)
q5_bray.mds.data <- metaMDS(as.dist(q5_braycurtis), trymax = 3000, trace = 0) %>% 
  scores() %>% as.data.frame() %>% mutate(cycles = factor(q5_cycles$cycles))

# Create Graph
ggplot(phu_bray.mds.data, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=cycles), size = 3) + 
  theme_bw()
    
ggplot(q5_bray.mds.data, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(color=cycles), size = 3) + 
  theme_bw() 
  
  

  
