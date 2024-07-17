###GENERATING GROWTH CURVES###
library(pacman)

p_load(data.table, tidyverse, hrbrthemes, growthcurver, dplyr, stringr, ggplot2, reshape2)



media_1 <- fread("media_growth_1.csv", header=TRUE)

media_1$Time_hours <- media_1$`Time [s]`/3600

# Identify all unique conditions (excluding the _1, _2, _3 suffixes)
conditions <- unique(str_remove(colnames(media_1)[str_detect(colnames(media_1), "_[1-3]$")], "_[1-3]$"))

# Subtract blanks
media_1_adjusted <- media_1 %>%
  mutate(across(contains("_1"), ~ . - get(paste0(str_remove(cur_column(), "_1"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(contains("_2"), ~ . - get(paste0(str_remove(cur_column(), "_2"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(contains("_3"), ~ . - get(paste0(str_remove(cur_column(), "_3"), "_blank")), .names = "adjusted_{col}")) %>%
  select(Time_hours, starts_with("adjusted_"))

colnames(media_1_adjusted) <- str_replace(colnames(media_1_adjusted), "^adjusted_", "")

media_long <- media_1_adjusted %>%
  pivot_longer(
    cols = -Time_hours, 
    names_to = c("condition", "replicate"), 
    names_pattern = "^(.*)_(\\d)$", 
    values_to = "adjusted_value"
  )

ggplot(media_long, aes(x = Time_hours, y = adjusted_value, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = condition), alpha = 0.2, color = NA) +
  labs(title = "Growth Curve with standard error",
       x = "Time (hours)",
       y = "OD (mean ± 95% CI)",
       color = "Condition",
       fill = "Condition") +
  theme_minimal()

################################

media_succ_acid <- fread("media_succinic_acid_1.csv", header=TRUE)


media_succ_acid$Time_hours <- as.numeric(media_succ_acid$`Time [s]`)/3600

# Identify all unique conditions (excluding the _1, _2, _3 suffixes)
conditions <- unique(str_remove(colnames(media_succ_acid)[str_detect(colnames(media_succ_acid), "_[1-3]$")], "_[1-3]$"))

# Subtract blanks
media_succ_acid_adjusted <- media_succ_acid %>%
  mutate(across(ends_with("_1"), ~ . - get(paste0(str_remove(cur_column(), "_1"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(ends_with("_2"), ~ . - get(paste0(str_remove(cur_column(), "_2"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(ends_with("_3"), ~ . - get(paste0(str_remove(cur_column(), "_3"), "_blank")), .names = "adjusted_{col}")) %>%
  select(Time_hours, starts_with("adjusted_"))


colnames(media_succ_acid_adjusted) <- str_replace(colnames(media_succ_acid_adjusted), "^adjusted_", "")

media_long <- media_succ_acid_adjusted %>%
  pivot_longer(
    cols = -Time_hours, 
    names_to = c("condition", "replicate"), 
    names_pattern = "^(.*)_(\\d)$", 
    values_to = "adjusted_value"
  )

ggplot(media_long, aes(x = Time_hours, y = adjusted_value, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = condition), alpha = 0.2, color = NA) +
  labs(title = "Growth Curve with SE",
       x = "Time (hours)",
       y = "OD (mean ± SE)",
       color = "Condition",
       fill = "Condition") +
  theme_minimal()

#######################################

media_sod_succ <- fread("media_sod_succ_1.csv", header=TRUE)


media_sod_succ$Time_hours <- as.numeric(media_sod_succ$`Time [s]`)/3600

# Identify all unique conditions (excluding the _1, _2, _3 suffixes)
conditions <- unique(str_remove(colnames(media_sod_succ)[str_detect(colnames(media_sod_succ), "_[1-3]$")], "_[1-3]$"))

# Subtract blanks
media_sod_succ_adjusted <- media_sod_succ %>%
  mutate(across(ends_with("_1"), ~ . - get(paste0(str_remove(cur_column(), "_1"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(ends_with("_2"), ~ . - get(paste0(str_remove(cur_column(), "_2"), "_blank")), .names = "adjusted_{col}")) %>%
  mutate(across(ends_with("_3"), ~ . - get(paste0(str_remove(cur_column(), "_3"), "_blank")), .names = "adjusted_{col}")) %>%
  select(Time_hours, starts_with("adjusted_"))


colnames(media_sod_succ_adjusted) <- str_replace(colnames(media_sod_succ_adjusted), "^adjusted_", "")

media_long <- media_sod_succ_adjusted %>%
  pivot_longer(
    cols = -Time_hours, 
    names_to = c("condition", "replicate"), 
    names_pattern = "^(.*)_(\\d)$", 
    values_to = "adjusted_value"
  )

ggplot(media_long, aes(x = Time_hours, y = adjusted_value, color = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = condition), alpha = 0.2, color = NA) +
  labs(title = "Growth Curve with SE",
       x = "Time (hours)",
       y = "OD (mean ± SE)",
       color = "Condition",
       fill = "Condition") +
  theme_minimal()
