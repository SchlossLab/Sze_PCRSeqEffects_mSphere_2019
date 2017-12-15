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


# Vector of sub samples used
#sub_sample_level <- c("50", "100", "500", "1000", "5000", "10000")

# Read in the count data
error_data <- read_csv("data/process/tables/full_error_summary.csv")

#error_data <- sapply(sub_sample_level, 
#                     function(x) read_data("data/process/tables/error_", "", "_summary.csv", x), 
#                     simplify = F)


# Generate graph of Mock DNA samples 
full_graph <- error_data %>% 
  filter(total_seqs >= 1000) %>% 
  mutate(taq = factor(taq, 
                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
  group_by(taq, cycles) %>% 
  summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
            group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
  ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
  geom_line() + 
  geom_errorbar(aes(ymin=group_iqr25, 
                    ymax=group_iqr75), 
                width = 0.1, size = 0.5, alpha = 0.4) + 
  geom_point(size = 2, alpha = 0.7) + 
  theme_bw() + 
  scale_color_manual(name = "Taq Used", 
                     values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
  coord_cartesian(ylim = c(0, 0.10)) + 
  scale_y_continuous(labels = scales::percent) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        legend.position = c(0.10, 0.83), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(color = "black"))






#thousand <- error_data[["1000"]] %>% 
#  mutate(taq = factor(taq, 
#                      levels = c("ACC", "K", "PHU", "PL", "Q5"), 
#                      labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
#  group_by(taq, cycles) %>% 
#  summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
#            group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
#  ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
#  geom_line() + 
  # geom_errorbar(aes(ymin=group_iqr25, 
  #                   ymax=group_iqr75), 
  #               width = 0.1, size = 0.5, alpha = 0.4) + 
  # geom_point(size = 2, alpha = 0.7) + 
  # theme_bw() + 
  # scale_color_manual(name = "Taq Used", 
  #                    values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
  # labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
  # coord_cartesian(ylim = c(0, 0.10)) + 
  # scale_y_continuous(labels = scales::percent) + 
  # ggtitle("A") + 
  # annotate("text", label = paste("Sub-sampled to 1000 Sequences"), x = 1.7, y = 0.102, size = 2.5) + 
  # theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(), 
  #       axis.text.y = element_text(size = 10), 
  #       legend.position = c(0.25, 0.83), 
  #       legend.title = element_blank(), 
  #       legend.key = element_blank(), 
  #       legend.background = element_rect(color = "black"))


# five_thousand <- error_data[["5000"]] %>% 
#   mutate(taq = factor(taq, 
#                       levels = c("ACC", "K", "PHU", "PL", "Q5"), 
#                       labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
#   group_by(taq, cycles) %>% 
#   summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
#             group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
#   ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
#   geom_line() + 
#   geom_errorbar(aes(ymin=group_iqr25, 
#                     ymax=group_iqr75), 
#                 width = 0.1, size = 0.5, alpha = 0.4) + 
#   geom_point(size = 2, alpha = 0.7) + 
#   theme_bw() + 
#   scale_color_manual(name = "Taq Used", 
#                      values = c("#440154FF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
#   labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
#   coord_cartesian(ylim = c(0, 0.10)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   ggtitle("B") + 
#   annotate("text", label = paste("Sub-sampled to 5000 Sequences"), x = 1.7, y = 0.102, size = 2.5) + 
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         legend.position = c(0.25, 0.85), 
#         legend.title = element_blank(), 
#         legend.key = element_blank(), 
#         legend.background = element_rect(color = "black"))



# ten_thousand <- error_data[["10000"]] %>% 
#   mutate(taq = factor(taq, 
#                       levels = c("ACC", "K", "PHU", "PL", "Q5"), 
#                       labels = c("Accuprime", "Kappa", "Phusion", "Platinum", "Q5"))) %>% 
#   group_by(taq, cycles) %>% 
#   summarise(group_chimera = median(chimera_prevalence), group_iqr25 = quantile(chimera_prevalence)["25%"], 
#             group_iqr75 = quantile(chimera_prevalence)["75%"]) %>% 
#   ggplot(aes(cycles, group_chimera, color = taq, group = taq)) + 
#   geom_line() + 
#   geom_errorbar(aes(ymin=group_iqr25, 
#                     ymax=group_iqr75), 
#                 width = 0.1, size = 0.5, alpha = 0.4) + 
#   geom_point(size = 2, alpha = 0.7) + 
#   theme_bw() + 
#   scale_color_manual(name = "Taq Used", 
#                      values = c("#440154FF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + 
#   labs(x = "Amplification Cycles", y = "Chimeric Sequences") + 
#   coord_cartesian(ylim = c(0, 0.1)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   ggtitle("C") + 
#   annotate("text", label = paste("Sub-sampled to 10000 Sequences"), x = 1.7, y = 0.102, size = 2.5) + 
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.y = element_text(size = 10), 
#         legend.position = c(0.25, 0.85), 
#         legend.title = element_blank(), 
#         legend.key = element_blank(), 
#         legend.background = element_rect(color = "black"))


#combined_graph <- grid.arrange(thousand, five_thousand, ten_thousand, ncol = 3)

ggsave("results/figures/Figure4.pdf", full_graph, width = 7, height = 6, dpi = 300)

