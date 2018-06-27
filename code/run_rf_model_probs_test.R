### Test RF model probs obtained from logloss values
## Test classification of fecal sample versus polymerase and cycles
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test", "MLmetrics", "magrittr"))


combined_data <- read_csv("data/process/tables/sample_group_rf_model_summary.csv") %>% 
  mutate(model_type = "fecal_sample_group") %>% 
  bind_rows(read_csv("data/process/tables/polymerase_group_rf_model_summary.csv") %>% 
              mutate(model_type = "fecal_polymerase"), 
            read_csv("data/process/tables/cycles_group_rf_model_summary.csv") %>% 
              mutate(model_type = "fecal_cycles"), 
            read_csv("data/process/tables/polymerase_group_mock_rf_model_summary.csv") %>% 
              mutate(model_type = "mock_polymerase"), 
            read_csv("data/process/tables/cycles_group_mock_rf_model_summary.csv") %>% 
              mutate(model_type = "mock_cycles")) %>% 
  mutate(probs = exp(-logLoss)) %>% 
  select(model_type, probs)


# Run the test to compare how good the models are
initial_kruskal_test <- kruskal.test(probs ~ factor(model_type), data = combined_data) %>% 
  broom::tidy()

# Run the Dunn post hoc test 
post_hoc_test <- t(as.data.frame.list(dunn.test(combined_data$probs, factor(combined_data$model_type), method = "bh"))) %>% 
  broom::tidy() %>% t() %>% tbl_df() %>% 
  rename(chi2 = V1, Z = V2, pvalue = V3, bh = V4, comparison = V5) %>% 
  slice(2:length(rownames(.))) %>% 
  separate(comparison, c("var1", "var2"), sep = "-")

# Write out the tests
write_csv(initial_kruskal_test, "data/process/tables/kruskal_rf_model_test.csv")
write_csv(post_hoc_test, "data/process/tables/dunn_rf_model_test_summary.csv")


  















