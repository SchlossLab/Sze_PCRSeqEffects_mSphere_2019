### Run a random forest analysis
## Test classification of fecal sample versus polymerase
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse", "stringr", "viridis", "gridExtra", "caret", "dunn.test", "MLmetrics"))


# Load in the needed the data
shared <- read_tsv("data/process/all_amp.0.03.subsample.1000.shared")

metadata <- shared %>% select(Group) %>% 
  separate(Group, c("x1", "x2", "x3", "x4", "x5")) %>% 
  mutate(sample_id = pull(shared, Group)) %>% 
  filter(x5 %in% c("DA10001", "DA10016", "DA10029", "DA10040")) %>% 
  rename(cycles = x1, taq = x2, sample_type = x3, sample_number = x4, sample_name = x5) %>% 
  mutate(cycles = case_when(
    cycles == "15x" ~ "fifteen", 
    cycles == "20x" ~ "twenty", 
    cycles == "25x" ~ "twenty_five", 
    cycles == "30x" ~ "thirty", 
    cycles == "35x" ~ "thirty_five", 
    TRUE ~ "uh oh"))

taxonomy <- read_tsv("data/process/all_amp.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

# Remove the nzv data
nzv <- nearZeroVar(as.data.frame(shared))

test_shared <- as.data.frame(shared)[, -nzv]


# create function to run the model
run_rf <- function(shared_file, metadata_file, selection){
  
  # Create a random forest run to identify key or important taxa
  run_data <- shared_file %>% 
    inner_join(select(metadata_file, sample_id, selection), by = c("Group" = "sample_id")) %>% 
    select(selection, contains("Otu"))
  
  
  fitControl <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    p = 0.8,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final")

  number_try <- round(sqrt(ncol(run_data)))

  # Set the mtry hyperparameter for the training model
  tunegrid <- expand.grid(.mtry = number_try)

  imp_vars_results <- c()
  model_results <- c()

  for(i in 1:100){

    print(paste("Running model ", i, "...", sep = ""))

    #Train the model
    training_model <-
      train(formula(paste(selection, " ~ .", sep = "")), data = run_data,
            method = "rf",
            ntree = 500,
            importance = TRUE,
            trControl = fitControl,
            tuneGrid = tunegrid,
            metric = "logLoss",
            na.action = na.omit,
            verbose = FALSE)

    model_results <- bind_rows(model_results, training_model$results %>% mutate(model_run = i))

    imp_vars_results <- bind_rows(imp_vars_results, varImp(training_model, scale = F)$importance %>%
                                    mutate(otu = rownames(.),
                                           model_run = i))
  }


  summary_imp_vars <- imp_vars_results %>% 
    gather(groups, mda, -otu, -model_run) %>% 
    group_by(otu) %>%
    summarise(median_mda = median(mda),
              iqr25 = quantile(mda)["25%"],
              iqr75 = quantile(mda)["75%"],
              min_mda = min(mda),
              max_mda = max(mda)) %>%
    arrange(desc(median_mda)) %>%
    left_join(taxonomy, by = c("otu" = "OTU")) %>%
    select(otu, genus, median_mda, iqr25, iqr75, min_mda, max_mda)

  temp_list <- list(model_results = model_results,
                    raw_imp_data = imp_vars_results,
                    summary_imp_data = summary_imp_vars)

  return(temp_list)

}

# Run the RF runs for each respective class separation
sample_groups <- run_rf(test_shared, metadata, "sample_name")

polymerase_groups <- run_rf(test_shared, metadata, "taq")

cycle_groups <- run_rf(test_shared, metadata, "cycles")

# Pull out the summary important variables for respective files from the analysis
pg_imp_vars_summary <- polymerase_groups$summary_imp_data %>% 
  slice(1:10)
cg_imp_vars_summary <- cycle_groups$summary_imp_data %>% 
  slice(1:10)

sg_imp_vars_summary <- sample_groups$summary_imp_data %>% 
  mutate(other_model = case_when(
    otu %in% pull(pg_imp_vars_summary, otu) & !(otu %in% pull(cg_imp_vars_summary, otu)) ~ "polymerase", 
    otu %in% pull(cg_imp_vars_summary, otu) & !(otu %in% pull(pg_imp_vars_summary, otu)) ~ "cycles", 
    otu %in% pull(cg_imp_vars_summary, otu) & otu %in% pull(pg_imp_vars_summary, otu) ~ "both", 
    !(otu %in% pull(pg_imp_vars_summary, otu)) & !(otu %in% pull(cg_imp_vars_summary, otu)) ~ "neither", 
    TRUE ~ "uh oh"))

#Pull out the summary model runs for each model
pg_run_summary <- polymerase_groups$model_results
cg_run_summary <- cycle_groups$model_results
sg_run_summary <- sample_groups$model_results

  # AUC in this case refers to area under the precision-recall curve
  # Logloss is used for training instead of AUC (ROC) 
    # Logloss closer to 0 is good - exp(-logloss) = average probablilty of attribution to the right class 

# Write out the results
write_csv(sg_imp_vars_summary, "data/process/tables/sample_group_rf_imp_vars_summary.csv")
write_csv(pg_imp_vars_summary, "data/process/tables/polymerase_group_top10_rf_imp_vars_summary.csv")
write_csv(cg_imp_vars_summary, "data/process/tables/cycle_group_top10_rf_imp_vars_summary.csv")

write_csv(pg_run_summary, "data/process/tables/polymerase_group_rf_model_summary.csv")
write_csv(cg_run_summary, "data/process/tables/cycles_group_rf_model_summary.csv")
write_csv(sg_run_summary, "data/process/tables/sample_group_rf_model_summary.csv")

