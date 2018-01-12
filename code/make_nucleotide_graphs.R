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


# Function to gatger the data
run_data_gather <- function(i, dataList){
  
  tempData <- dataList[[i]] %>% 
    gather("sub_type", "rate", contains("_rate")) %>% 
    mutate(sub_type = str_replace(sub_type, "_rate", "")) %>% 
    gather("prev_subtype", "seq_prevalence", contains("_seq_prev")) %>% 
    mutate(prev_subtype = str_replace(prev_subtype, "_seq_prev", ""))
  
  return(tempData)
  
}


###########################################################################################################################
############################### Run actual analysis programs  #############################################################
###########################################################################################################################

# Read in the count data
nucleotide_data <- read_csv("data/process/tables/full_nucleotide_error_summary.csv") %>% 
  filter(total_seqs >= 1000) %>% 
  gather("sub_type", "rate", contains("_rate")) %>% 
  mutate(sub_type = str_replace(sub_type, "_rate", "")) %>% 
  gather("prev_subtype", "seq_prevalence", contains("_seq_prev")) %>% 
  mutate(prev_subtype = str_replace(prev_subtype, "_seq_prev", ""))


full_graph <- nucleotide_data %>% 
     mutate(sub_type = factor(sub_type, 
                          levels = c("ac", "ag", "at", "ca", "cg", "ct", "ga", "gc", "gt", "ta", "tc", "tg"), 
                          labels = c("A => C", "A => G", "A => T", "C => A", "C => G", "C => T", "G => A", "G => C", 
                                     "G => T", "T => A", "T => C", "T => G"))) %>% 
     ggplot(aes(sub_type, log10(rate + 0.000001), color = taq, group = taq)) + 
     geom_point(size = 2, alpha = 0.2, position = position_dodge(width = 0.7), show.legend = F, color = "gray") + 
     theme_bw() + 
     stat_summary(aes(color = taq), fun.y = median, fun.ymin = median, fun.ymax = median, 
                  position = position_dodge(width = 0.7), geom = "crossbar", size = 0.5, width = 0.7) + 
     facet_grid(cycles~.) + 
     scale_color_manual(name = "Taq Used", 
                      values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
   labs(x = "", y = expression(Log["10"]~Substitution~Rate)) + 
   coord_cartesian(ylim = c(0, -6)) + 
   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 10), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         axis.text.y = element_text(size = 10), 
         legend.position = "bottom", 
         legend.title = element_blank(), 
         legend.key = element_blank(), 
         legend.background = element_rect(color = "black"))


# # Vector of sub samples used
# sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")
# 
# # Read in the count data
# nucleotide_data <- sapply(sub_sample_level, 
#                      function(x) read_data("data/process/tables/nucleotide_error_", "", "_summary.csv", x), 
#                      simplify = F)
# 
# # Transform data for better graphing
# nucleotide_data <- sapply(sub_sample_level, 
#                function(x) run_data_gather(x, nucleotide_data), simplify = F)
# 
# 
# 
# # Generate graph of Mock DNA samples 
# thousand <- nucleotide_data[["1000"]] %>% 
#   mutate(sub_type = factor(sub_type, 
#                        levels = c("ac", "ag", "at", "ca", "cg", "ct", "ga", "gc", "gt", "ta", "tc", "tg"), 
#                        labels = c("A => C", "A => G", "A => T", "C => A", "C => G", "C => T", "G => A", "G => C", 
#                                   "G => T", "T => A", "T => C", "T => G"))) %>% 
#   ggplot(aes(sub_type, log10(rate + 0.000001), color = taq, group = taq)) + 
#   geom_point(size = 2, alpha = 0.2, position = position_dodge(width = 0.7), show.legend = F, color = "gray") + 
#   theme_bw() + 
#   stat_summary(aes(color = taq), fun.y = median, fun.ymin = median, fun.ymax = median, 
#                position = position_dodge(width = 0.7), geom = "crossbar", size = 0.5, width = 0.7) + 
#   facet_grid(cycles~.) + 
#   scale_color_manual(name = "Taq Used", 
#                      values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
#   labs(x = "", y = expression(Log["10"]~Substitution~Rate)) + 
#   coord_cartesian(ylim = c(0, -6)) + 
#   ggtitle("Sub-Sampled to 1000 Sequences") +  
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 10), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         legend.position = "bottom", 
#         legend.title = element_blank(), 
#         legend.key = element_blank(), 
#         legend.background = element_rect(color = "black"))
#   
#   
#   
# five_thousand <- nucleotide_data[["5000"]] %>% 
#   mutate(sub_type = factor(sub_type, 
#                            levels = c("ac", "ag", "at", "ca", "cg", "ct", "ga", "gc", "gt", "ta", "tc", "tg"), 
#                            labels = c("A => C", "A => G", "A => T", "C => A", "C => G", "C => T", "G => A", "G => C", 
#                                       "G => T", "T => A", "T => C", "T => G"))) %>% 
#   ggplot(aes(sub_type, log10(rate + 0.000001), color = taq, group = taq)) + 
#   geom_point(size = 2, alpha = 0.2, position = position_dodge(width = 0.7), show.legend = F, color = "gray") + 
#   theme_bw() + 
#   stat_summary(aes(color = taq), fun.y = median, fun.ymin = median, fun.ymax = median, 
#                position = position_dodge(width = 0.7), geom = "crossbar", size = 0.5, width = 0.7) + 
#   facet_grid(cycles~.) + 
#   scale_color_manual(name = "Taq Used", 
#                      values = c("#440154FF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
#   labs(x = "", y = expression(Log["10"]~Substitution~Rate)) + 
#   coord_cartesian(ylim = c(0, -6)) + 
#   ggtitle("Sub-Sampled to 5000 Sequences") + 
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 10), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         legend.position = "bottom", 
#         legend.title = element_blank(), 
#         legend.key = element_blank(), 
#         legend.background = element_rect(color = "black"))
# 
# 
# ten_thousand <- nucleotide_data[["10000"]] %>% 
#   mutate(sub_type = factor(sub_type, 
#                            levels = c("ac", "ag", "at", "ca", "cg", "ct", "ga", "gc", "gt", "ta", "tc", "tg"), 
#                            labels = c("A => C", "A => G", "A => T", "C => A", "C => G", "C => T", "G => A", "G => C", 
#                                       "G => T", "T => A", "T => C", "T => G"))) %>% 
#   ggplot(aes(sub_type, log10(rate + 0.000001), color = taq, group = taq)) + 
#   geom_point(size = 2, alpha = 0.2, position = position_dodge(width = 0.7), show.legend = F, color = "gray") + 
#   theme_bw() + 
#   stat_summary(aes(color = taq), fun.y = median, fun.ymin = median, fun.ymax = median, 
#                position = position_dodge(width = 0.7), geom = "crossbar", size = 0.5, width = 0.7) + 
#   facet_grid(cycles~.) + 
#   scale_color_manual(name = "Taq Used", 
#                      values = c("#440154FF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
#   labs(x = "", y = expression(Log["10"]~Substitution~Rate)) + 
#   coord_cartesian(ylim = c(0, -6)) + 
#   ggtitle("Sub-Sampled to 10000 Sequences") + 
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 10), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         legend.position = "bottom", 
#         legend.title = element_blank(), 
#         legend.key = element_blank(), 
#         legend.background = element_rect(color = "black"))


ggsave("results/figures/FigureS2.pdf", full_graph, width = 8, height = 7, dpi = 300)




