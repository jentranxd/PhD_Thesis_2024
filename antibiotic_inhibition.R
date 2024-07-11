require("pacman")

p_load(data.table, dplyr, tidyr, purrr, DescTools)

setwd("C:/Users/jentr/OneDrive - UW-Madison/lab things/data/growth curves")

df <- fread("antibiotic_inhibition_data.csv", header=TRUE)


# Split the dataframe into a list of dataframes based on the plate number
split_df_by_plate <- function(df) {
  plate_numbers <- unique(gsub("(_.*)|(_time)", "", names(df)))
  plate_dfs <- map(plate_numbers, function(plate) {
    df %>%
      select(matches(paste0(plate, "_"))) %>%
      rename_with(~ gsub(paste0(plate, "_"), "", .), -matches(paste0(plate, "_time"))) %>%
      rename(time = matches(paste0(plate, "_time")))
  })
  names(plate_dfs) <- plate_numbers
  plate_dfs
}

# Subtract blanks from each corresponding antibiotic | concentration | replicate at each time point
subtract_blanks <- function(plate_df) {
  blank_cols <- grep("blank", names(plate_df), value = TRUE)
  plate_df <- plate_df %>%
    mutate(across(-time, ~ . - rowMeans(across(all_of(blank_cols)))))
  plate_df <- plate_df %>%
    select(-matches("blank"))
  plate_df
}


# Calculate AUC for each replicate
calculate_auc <- function(plate_df) {
  plate_df %>%
    pivot_longer(-time, names_to = c("antibiotic", "concentration", "replicate"), names_sep = "_") %>%
    group_by(antibiotic, concentration, replicate) %>%
    summarise(auc = AUC(time, value), .groups = 'drop')
}

# Calculate percent inhibition
calculate_percent_inhibition <- function(auc_df) {
  control_auc <- auc_df %>% filter(concentration == 0) %>% select(antibiotic, replicate, auc)
  auc_df <- auc_df %>%
    left_join(control_auc, by = c("antibiotic", "replicate"), suffix = c("", "_control")) %>%
    mutate(percent_inhibition = 100 * (1 - (auc / auc_control))) %>%
    select(-auc_control)
  auc_df
}

# Calculate average % inhibition and standard deviation
summarize_percent_inhibition <- function(inhibition_df, plate_number) {
  inhibition_df %>%
    group_by(antibiotic, concentration) %>%
    summarise(
      mean_percent_inhibition = mean(percent_inhibition),
      sd_percent_inhibition = sd(percent_inhibition),
      .groups = 'drop'
    ) %>%
    mutate(plate_number = plate_number)
}

# Main function to process the dataframe
process_growth_curves <- function(df) {
  plate_dfs <- split_df_by_plate(df)
  inhibition_summaries <- map2_dfr(plate_dfs, names(plate_dfs), function(plate_df, plate_number) {
    plate_df <- subtract_blanks(plate_df)
    auc_df <- calculate_auc(plate_df)
    inhibition_df <- calculate_percent_inhibition(auc_df)
    summarize_percent_inhibition(inhibition_df, plate_number)
  })
  return(inhibition_summaries)
}

# Assuming df is your dataframe
result_df <- process_growth_curves(df)
result_df

# View(result_df)

cleaned_results <- result_df %>% filter(concentration != 0)

fwrite(cleaned_results, "percent_inhibition_with_error.csv", sep = ",")
