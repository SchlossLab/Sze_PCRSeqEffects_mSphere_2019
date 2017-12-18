### Create graphs of the different numOTUs based on cycle number and sub sampling
### Mock Sample Community only
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
    filter(sample_type != "FS")
  # returns the file to the global environment
  return(tempData)
  
}


# Create unique metadata files for each sub-sampling
paredown_meta <- function(depth, metafile, distfile){
  
  tempvalues <- rownames(distfile[[depth]])
  
  metafile <- metafile %>% filter(full_name %in% tempvalues)
  
  return(metafile)
  
}

# Pare down distance matrix to match metadata file for each sub-sampling
paredown_dist <- function(depth, metaList, distList){
  
  tempMeta_values <- as.data.frame(metaList[[depth]])[, "full_name"]
  
  tempData <- distList[[depth]][tempMeta_values, tempMeta_values]
  
  
  return(tempData)
}


# Function to seperate each DNA polymerase within each sub sampling
separate_polymerase <- function(depth, metaList, distList, 
                                taq_used = c("ACC", "K", "PHU", "PL", "Q5"), meta_readout = F){
  
  tempMeta <- sapply(taq_used, 
                     function(x) filter(metaList[[depth]], taq == x), simplify = F)
  
  namesList <- sapply(taq_used, function(x) as.data.frame(tempMeta[[x]])[, "full_name"])
  
  tempDist <- sapply(taq_used, 
                     function(x) distList[[depth]][namesList[[x]], namesList[[x]]], simplify = F)
  
  if(meta_readout == F){
    
    return(tempDist)
  } else{
    
    return(tempMeta)
  }
  
}


# Function that pulls distances from the same individual but 5 cycles later from table 
get_next_distance <- function(depth, metaList, distList, 
                              taq_used = c("ACC", "K", "PHU", "PL", "Q5")){
  
  tempMeta <- metaList[[depth]]
  tempDist <- distList[[depth]]
  
  final_data <- sapply(taq_used, 
                       function(x) pull_dist_data(x, tempMeta, tempDist), simplify = F)
  
  final_data <- try(final_data %>% bind_rows())
  
  return(as.data.frame(final_data))
  
}


# Function that actuall gets the distances of the between current and previous cycle
pull_dist_data <- function(taq_of_int, metaList, distList){
  
  tempMeta <- metaList[[taq_of_int]]
  tempDist <- distList[[taq_of_int]]
  
  unique_names <- unique(tempMeta$sample_name)
  
  tempData <- sapply(unique_names, 
                     function(x) run_each_sample(x, taq_of_int, tempMeta, tempDist), simplify = F)
  
  tempData <- try(tempData %>% bind_rows())
  
  return(tempData)
  
}


# Function that controls the flow within a single sample 
run_each_sample <- function(u_names, taq_used, meta_file, distList){
  
  temp_names <- as.data.frame(filter(meta_file, sample_name == u_names))[, "full_name"]
  
  row_values <- temp_names[1:(length(temp_names)-1)]
  
  col_values <- temp_names[2:length(temp_names)]
  
  tempData <- t(sapply(c(1:length(row_values)), 
                       function(x) 
                         c(initial = row_values[x], 
                           final = col_values[x],
                           sample_name = u_names, 
                           taq = taq_used, 
                           distance = try(distList[row_values[x], col_values[x]])), simplify = T))
  
  if("distance" %in% colnames(tempData)){
    
    tempData <- as.data.frame(tempData)
  } else {
    
    tempData <- as.data.frame(tempData) %>% 
      mutate(distance = NA)
  }
  
  if(grepl("Error", tempData$distance) == T){
    
    tempData <- as.data.frame(tempData) %>% 
      mutate(distance = NA)
  }
  
  
  return(tempData)
  
  
}


# Function to expand meta data tables into a nice table
make_nice_table <- function(depth, dataList){
  
  tempData <- dataList[[depth]] %>% 
    select(-taq, -sample_name) %>% 
    separate(initial, c("lower_cycle", "taq", "kit", "sample_type", "sample_name")) %>% 
    select(-taq, -kit, -sample_type, -sample_name) %>% 
    separate(final, c("upper_cycle", "taq", "kit", "sample_type", "sample_name")) %>% 
    select(-kit, -sample_type) %>% 
    mutate(
      cycle_compare = paste(
        str_replace(lower_cycle, "x", "") , "to", str_replace(upper_cycle, "x", ""), sep = ""), 
      distance = as.numeric(distance)) %>% 
    select(cycle_compare, taq, sample_name, distance)
  
  
  return(tempData)
}


# Function to run the krukal wallis test and store the P-values
run_kruskal <- function(depth, dataList){
  
  taq_used <- c("ACC", "K", "PHU", "PL", "Q5")
  
  tempTest <- sapply(taq_used, 
                     function(x) try(kruskal.test(distance ~ factor(cycle_compare), 
                                                  data = filter(dataList[[depth]], taq == x))$p.value), simplify = F)
  
  tempTest <- as.data.frame(t(tempTest %>% 
                                bind_rows())) %>% 
    mutate(taq = rownames(.), 
           pvalue = V1, 
           pvalue = ifelse(grepl("Error", pvalue) == T, invisible(NA), 
                           invisible(as.numeric(as.character(pvalue)))), 
           bh = p.adjust(pvalue, method = "BH"), 
           sub_sample_level = depth) %>% 
    select(taq, sub_sample_level, pvalue, bh)
  
  
  return(tempTest)
}

# Create function to run PERMANOVA analysis by subsampling by Taq
get_permanova <- function(depth, metaList, dataList){
  
  taq_used <- c("ACC", "K", "PHU", "PL", "Q5")
  
  tempMeta <- metaList[[depth]]
  tempData <- dataList[[depth]]
  
  set.seed(12345)
  tempTest <- sapply(taq_used, 
                     function(x) 
                       try(adonis(as.dist(tempData[[x]]) ~ factor(tempMeta[[x]]$cycles))$aov.tab$`Pr(>F)`[1]), 
                     simplify = F) %>% bind_rows() %>% t() %>% as.data.frame() %>% 
    mutate(taq = rownames(.), 
           pvalue = V1, 
           pvalue = ifelse(grepl("Error", pvalue) == T, invisible(NA), invisible(as.character(pvalue))), 
           pvalue = as.numeric(pvalue), 
           bh = p.adjust(pvalue, method = "BH"), 
           sub_sample_level = depth) %>% 
    select(taq, sub_sample_level, pvalue, bh)
  
  return(tempTest)
  
}


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Vector of sub samples used
sub_sample_level <- c("1000", "5000", "10000")


# Read in distance matrix
braycurtis_dist <- sapply(sub_sample_level, 
                          function(x) read.dist(paste("data/process/all_amp.braycurtis.0.03.lt.", 
                                                      x, ".dist", sep = "")), simplify = F)


# Read in meta data
metadata <- read_data("data/process/tables/", "meta_data", ".csv", "")

# Generate pared down meta matching each sub-sampling
meta_list <- sapply(sub_sample_level, 
                    function(x) paredown_meta(x, metadata, braycurtis_dist), simplify = F)

# Generate pared down dist list for each sub-sampling that matches meta-data
red_bray_dist <- sapply(sub_sample_level, 
                        function(x) paredown_dist(x, meta_list, braycurtis_dist), simplify = F)


# Separate each polymerase within each sub-sampling
sep_red_bray_dist <- sapply(sub_sample_level, 
                            function(x) separate_polymerase(x, meta_list, braycurtis_dist), simplify = F)

sep_meta_list <- sapply(sub_sample_level, 
                        function(x) separate_polymerase(x, meta_list, 
                                                        braycurtis_dist, meta_readout = T), simplify = F)


next_distance_values <- sapply(sub_sample_level, 
                               function(x) get_next_distance(x, sep_meta_list, sep_red_bray_dist), simplify = F)

finalized_tables <- sapply(sub_sample_level, 
                           function(x) make_nice_table(x, next_distance_values), simplify = F)

# Get kruskal-wallis results on similarity to last cycle group
kruskal_test_data <- sapply(sub_sample_level, 
                            function(x) run_kruskal(x, finalized_tables), simplify = F) %>% 
  bind_rows()


# Get PERMANOVA results
permanova_results <- sapply(sub_sample_level, 
                            function(x) get_permanova(x, sep_meta_list, sep_red_bray_dist), simplify = F) %>% 
  bind_rows()

# Write out results 
write_csv(kruskal_test_data, "data/process/tables/bray_sim_to_prev_cycle_kruskal_results.csv")
write_csv(permanova_results, "data/process/tables/bray_permanova_by_taq_results.csv")

sapply(c(1:length(finalized_tables)), 
       function(x) write_csv(finalized_tables[[x]], 
                             paste("data/process/tables/bray_5_cycle_dist_", sub_sample_level[x], "_data.csv", sep = "")))