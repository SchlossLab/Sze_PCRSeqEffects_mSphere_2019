### Correlations between error and chimera with number of OTUs
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra"))


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
  tempData <- read_csv(paste(pathing, start_name, differentiator, end_name, sep = ""))
  # returns the file to the global environment
  return(tempData)
  
}


# FUnction to create a new common column
make_common_column <- function(i, dataList){
  
  tempData <- dataList[[i]] %>% 
    mutate(common_column = paste(cycles, "_", taq, "_", sample_name, sep = ""))
  
  return(tempData)
  
}

# Function to join data together
combine_data <- function(i, dataList, to_match_list){
  
  tempData <- dataList[[i]] %>% 
    left_join(to_match_list[[i]], by = "common_column")
  
  return(tempData)
    
}

# Function to make label for R2 (found on stack overflow)
lm_eqn = function(df){
  m = lm((chimera_prevalence*1000) ~ total_otus, df);
  eq <- substitute(~~R^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



# Read in subsample.shared files
error_files <- c("mock_precluster_error", "mock_chimera_error", "mock_error") 

# # Read in the count data
error_data <- sapply(error_files, function(x) read_csv(paste("data/process/tables/", x, "_summary.csv", sep = "")) %>% 
                       filter(total_seqs >= 1000), simplify = F)

# Read in the count data
test <- read_tsv("data/process/all_amp.shared") %>% 
  select(-label, -numOtus) %>% 
  filter(str_detect(Group, "FS") == FALSE)

test2 <- select(test, contains("Otu"))
test <- test2[, colSums(test2) > 0] %>% 
  mutate(Group = test$Group)

rm(test2)

test <- test %>%   
  mutate_at(vars(contains("Otu")), function(x) ifelse(x > 0, invisible(1), invisible(0))) %>% 
  mutate(total_otus = select(., -Group) %>% rowSums(.)) %>% 
  select(Group, total_otus)


num_otu_data <- sapply(error_files, 
                       function(x) 
                         error_data[[x]] %>% 
                         select(full_name, cycles, taq, sample_name) %>% 
                         left_join(test, by = c("full_name" = "Group")), simplify = F)
  

# Create new common column
up_error_data <- sapply(error_files, 
                        function(x) make_common_column(x, error_data), simplify = F)

up_numOTU_data <- sapply(error_files, 
                         function(x) make_common_column(x, num_otu_data), simplify = F)

# Combine the data together
combined_list <- sapply(error_files, 
                        function(x) combine_data(x, up_error_data, up_numOTU_data), simplify = F)


eqns_m_pc_error <- by(combined_list[["mock_precluster_error"]], 
               combined_list[["mock_precluster_error"]]$taq.x, lm_eqn)
df_m_pc_error<- data.frame(eq = unclass(eqns_m_pc_error), taq = rownames(eqns_m_pc_error)) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))

mock_pc <- combined_list[["mock_precluster_error"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(chimera_prevalence*100, total_otus, color = cycles.x, group = taq)) + 
  geom_smooth(size = 1, method = "lm", se = FALSE, color = "black") + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(taq ~.) + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#0000FF", "#00C957", "#CD8500", "#FF1493")) + 
  labs(x = "Chimera Sequence Prevalence (%)", y = "Number of OTUs") + 
  geom_text(data = df_m_pc_error, aes(x = 32, y = 50, label = eq), 
            color = 'black',  parse = TRUE) + 
  ggtitle("A") + coord_cartesian(ylim = c(0, 500)) + scale_x_continuous(limits = c(0, 40)) + 
  annotate("text", label = paste("Before Precluster Step"), x = 32, y = 110, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))



eqns_m_ch_error <- by(combined_list[["mock_chimera_error"]], 
               combined_list[["mock_chimera_error"]]$taq.x, lm_eqn)

df_m_ch_error <- data.frame(eq = unclass(eqns_m_ch_error), taq = rownames(eqns_m_ch_error)) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))

mock_ch <- combined_list[["mock_chimera_error"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(chimera_prevalence*100, total_otus, color = cycles.x, group = taq)) + 
  geom_smooth(size = 1, method = "lm", se = FALSE, color = "black") + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(taq ~.) + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#0000FF", "#00C957", "#CD8500", "#FF1493")) + 
  labs(x = "Chimera Sequence Prevalence (%)", y = "Number of OTUs") + 
  geom_text(data = df_m_ch_error, aes(x = 32, y = 50, label = eq), 
            color = 'black',  parse = TRUE) + 
  ggtitle("B") + coord_cartesian(ylim = c(0, 500)) + scale_x_continuous(limits = c(0, 40)) + 
  annotate("text", label = paste("Before Chimera Removal Step"), x = 32, y = 110, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


eqns_m_error <- by(combined_list[["mock_error"]], 
                      combined_list[["mock_error"]]$taq.x, lm_eqn)

df_m_error <- data.frame(eq = unclass(eqns_m_ch_error), taq = rownames(eqns_m_ch_error)) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))

mock_error <- combined_list[["mock_error"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(chimera_prevalence*100, total_otus, color = cycles.x, group = taq)) + 
  geom_smooth(size = 1, method = "lm", se = FALSE, color = "black") + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(taq ~.) + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#0000FF", "#00C957", "#CD8500", "#FF1493")) + 
  labs(x = "Chimera Sequence Prevalence (%)", y = "Number of OTUs") + 
  geom_text(data = df_m_error, aes(x = 32, y = 50, label = eq), 
            color = 'black',  parse = TRUE) + 
  ggtitle("C") +coord_cartesian(ylim = c(0, 500)) + scale_x_continuous(limits = c(0, 40)) + 
  annotate("text", label = paste("Full Pipeline"), x = 32, y = 110, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))

###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

combined_graph <- grid.arrange(mock_pc, mock_ch, mock_error, ncol = 3)

ggsave("results/figures/Figure6.pdf", combined_graph, width = 15, height = 8, dpi = 300)



