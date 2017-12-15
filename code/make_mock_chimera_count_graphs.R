### Create graphs of the number of sequences with a chimera based on cycle number and sub sampling
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


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################


error_files <- c("mock_precluster_error", "mock_chimera_error", "mock_error") 

# # Read in the count data
error_data <- sapply(error_files, function(x) read_csv(paste("data/process/tables/", x, "_summary.csv", sep = "")) %>% 
                       filter(total_seqs >= 1000), simplify = F)


# Generate graph of Mock DNA samples 
before_precluster <- error_data[["mock_precluster_error"]]  %>% 
  filter(total_seqs >= 1000) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  group_by(taq, cycles) %>% 
  summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
            group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
  ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
  geom_line(show.legend = T) + 
  geom_errorbar(aes(ymin=group_iqr25, 
                    ymax=group_iqr75), show.legend = T,  
                width = 0.1, size = 0.5, alpha = 0.4) + 
  geom_point(size = 2, alpha = 0.7, show.legend = T) + 
  theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
  coord_cartesian(ylim = c(0, 0.30)) + 
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("A") + 
  annotate("text", label = paste("Before Pre-Cluster Step"), x = 2.5, y = 0.31, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.25, 0.83), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


before_chimera_vsearch <- error_data[["mock_chimera_error"]]  %>% 
  filter(total_seqs >= 1000) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  group_by(taq, cycles) %>% 
  summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
            group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
  ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
  geom_line(show.legend = F) + 
  geom_errorbar(aes(ymin=group_iqr25, 
                    ymax=group_iqr75), show.legend = F, 
                width = 0.1, size = 0.5, alpha = 0.4) + 
  geom_point(size = 2, alpha = 0.7, show.legend = F) + 
  theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
  coord_cartesian(ylim = c(0, 0.30)) + 
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("B") + 
  annotate("text", label = paste("Before Chimera VSEARCH Step"), x = 2.5, y = 0.31, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.30, 0.83), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


full_pipeline <- error_data[["mock_error"]]  %>% 
  filter(total_seqs >= 1000) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  group_by(taq, cycles) %>% 
  summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
            group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
  ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
  geom_line(show.legend = F) + 
  geom_errorbar(aes(ymin=group_iqr25, 
                    ymax=group_iqr75), show.legend = F, 
                width = 0.1, size = 0.5, alpha = 0.4) + 
  geom_point(size = 2, alpha = 0.7, show.legend = F) + 
  theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
  coord_cartesian(ylim = c(0, 0.30)) + 
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("C") + 
  annotate("text", label = paste("Full Pipeline"), x = 2.5, y = 0.31, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.10, 0.83), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))


combined_graph <- grid.arrange(before_precluster, before_chimera_vsearch, full_pipeline, ncol = 3)

ggsave("results/figures/Figure4.pdf", combined_graph, width = 11, height = 7, dpi = 300)

