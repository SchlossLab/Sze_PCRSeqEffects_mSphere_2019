### Correlations between error and chimera with number of OTUs
### Mock Community only
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis"))


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
  m = lm((chimera_prevalence*1000) ~ numOTUs, df);
  eq <- substitute(~~R^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



# Read in subsample.shared files
sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in the count data
error_data <- sapply(sub_sample_level, 
                     function(x) read_data("data/process/tables/error_", "", "_summary.csv", x), 
                     simplify = F)


# Read in the count data
numOTU_data <- sapply(sub_sample_level, 
                      function(x) read_data("data/process/tables/", "mock_sub_sample_", "_count_table.csv", x), 
                      simplify = F)


# Create new common column
up_error_data <- sapply(sub_sample_level, 
                        function(x) make_common_column(x, error_data), simplify = F)

up_numOTU_data <- sapply(sub_sample_level, 
                         function(x) make_common_column(x, numOTU_data), simplify = F)

# Combine the data together
combined_list <- sapply(sub_sample_level, 
                        function(x) combine_data(x, up_error_data, up_numOTU_data), simplify = F)




eqns1000 <- by(combined_list[["1000"]], combined_list[["1000"]]$taq.x, lm_eqn)
df1000 <- data.frame(eq = unclass(eqns1000), taq = rownames(eqns1000)) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))

combined_list[["1000"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(chimera_prevalence*100, numOTUs, color = cycles.x, group = taq)) + 
  geom_smooth(size = 1, method = "lm", se = FALSE, color = "black") + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(taq ~.) + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#0000FF", "#00C957", "#CD8500", "#FF1493")) + 
  labs(x = "Percent Chimera Prevalence", y = "Number of OTUs") + 
  geom_text(data = df1000, aes(x = 0.5, y = 65, label = eq, family = "arial"), 
            color = 'black',  parse = TRUE) + 
  ggtitle("A") + coord_cartesian(ylim = c(0, 80)) + 
  annotate("text", label = paste("Sub-sampled to 1000 Sequences"), x = 24, y = 8.2, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))



eqns5000 <- by(combined_list[["5000"]], combined_list[["5000"]]$taq.x, lm_eqn)
df5000 <- data.frame(eq = unclass(eqns5000), taq = rownames(eqns5000)) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5")))

combined_list[["5000"]] %>% 
  mutate(taq = factor(taq.x, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  ggplot(aes(chimera_prevalence*100, numOTUs, color = cycles.x, group = taq)) + 
  geom_smooth(size = 1, method = "lm", se = FALSE, color = "black") + 
  geom_point(size = 2, alpha = 0.7) + theme_bw() + 
  facet_grid(taq ~.) + 
  scale_color_manual(name = "Cycle Number", 
                     values = c("#0000FF", "#00C957", "#CD8500", "#FF1493")) + 
  labs(x = "Percent Chimera Prevalence", y = "Number of OTUs") + 
  geom_text(data = df5000, aes(x = 0.5, y = 130, label = eq, family = "arial"), 
            color = 'black',  parse = TRUE) + 
  ggtitle("B") + coord_cartesian(ylim = c(0, 150)) + 
  annotate("text", label = paste("Sub-sampled to 5000 Sequences"), x = 14.5, y = 8.2, size = 2.5) + 
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





