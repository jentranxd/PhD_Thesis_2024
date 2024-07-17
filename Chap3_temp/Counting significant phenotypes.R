####COUNTING SIGNIFICANT PHENOTYPES#######
#######################################################
##########################################

#load
p_load(stringr, dplyr)


# # #filter for perfect guides on the gene-level
filtered_results <- median_melted_results %>% rename(locus_tag = AB19606) %>%
  filter(type == "perfect" & (medLFC >= 1 | medLFC <= -1) & FDR < 0.05)

# #filter for perfect guides on the gene-level with a positive LFC
filtered_results_pos <- median_melted_results %>% rename(locus_tag = AB19606) %>%
  filter(type == "perfect" & (medLFC >= 1) & FDR < 0.05)

# #filter for perfect guides on the gene-level with a negative LFC
filtered_results_neg <- median_melted_results %>% rename(locus_tag = AB19606) %>%
  filter(type == "perfect" & (medLFC <= -1) & FDR < 0.05)


#group by gene-level and count how many unique conditions there are
unique_condition_count <- filtered_results %>%
  group_by(locus_tag) %>%
  summarise(unique_condition_count = n_distinct(condition))

# Join unique_condition_count with curated_names based on locus_tag
merged_df <- unique_condition_count %>%
  left_join(select(curated_names, AB19606, unique_name, STRING_names), by = c("locus_tag" = "AB19606"))

# Create a dataframe with missing locus_tag and unique_name with count 0
missing_df <- anti_join(select(curated_names, AB19606, unique_name, STRING_names), unique_condition_count, by = c("AB19606" = "locus_tag")) %>%
  mutate(count = 0) %>%  rename(locus_tag = AB19606, unique_condition_count = count)

# Bind rows to get the final dataframe
sig_final_df <- bind_rows(merged_df, missing_df)

#View significant phenotype histogram
ggplot(sig_final_df, aes(x = `unique_condition_count`)) +
  geom_histogram(binwidth = 1, fill = "maroon", color = "grey30", alpha = 0.7) + 
  labs(title = "Significant phenotypes per gene \n(gene-level; FDR<0.05 and |LFC|>1)",
       x = "Number of significant phenotypes across conditions compared to induction",
       y = "Gene Count (out of 406)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12)
  ) 

#same thing but graph it as a proportion
sig_final_df$proportion <- sig_final_df$unique_condition_count/46

ggplot(sig_final_df, aes(x = `proportion`)) +
  geom_density(fill = "grey", color = "grey30", alpha = 0.7) + 
  labs(title = "Essential gene chemical responses",
       x = "Chemical response (proportion of significant chemical phenotypes across conditions)",
       y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12)
  ) 


#################################################
##combine pos and neg results
# Add a new column 'direction' to each dataframe
filtered_results_pos$direction <- "pos"
filtered_results_neg$direction <- "neg"
filtered_results <- rbind(filtered_results_pos, filtered_results_neg)

# # Calculate the number of unique conditions for each gene and direction
unique_condition_counts <- filtered_results %>%
  group_by(locus_tag, unique_name, direction) %>%
  summarise(count = n_distinct(condition)) %>%
  ungroup()

######
#include zeroes

# Creating a distinct list of locus_tag and unique_name combinations
distinct_combinations <- median_melted_results %>%
  select(locus_tag = AB19606, unique_name) %>%
  distinct()

# Creating a two-row dataframe for directions
directions_df <- data.frame(direction = c("pos", "neg"))

# Cross joining to get all combinations
all_combinations <- crossing(distinct_combinations, directions_df)

# Now, proceed with the join with your summarized data and fill in zeros
unique_cond_counts_w_zeroes <- all_combinations %>%
  full_join(unique_condition_counts, by = c("locus_tag", "unique_name", "direction")) %>%
  replace_na(list(count = 0))

# final_results now includes all combinations, with zero fill for missing counts


# #plot histogram
unique_condition_counts %>%
  filter(direction == "pos") %>%
  ggplot(aes(x = count, fill = direction)) +
  geom_histogram(binwidth = 1, color="black") +
  scale_fill_manual(values = c("pos" = "#6699cc", "neg" = "#c94c4c")) +
  geom_vline(xintercept = 0, size = 1.5, color = "black") +
  geom_hline(yintercept = 0, size = 1.5, color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none") +
  labs(y = "Gene count (total = 406)", x = "Chemical resistance across conditions (total = 46)",
       title = "Distribution of Gene-level Chemical Resistance") +
  coord_cartesian(ylim = c(-5, 60), xlim = c(0, 18))

unique_condition_counts %>%
  filter(direction == "neg") %>%
  ggplot(aes(x = count, fill = direction)) +
  geom_histogram(binwidth = 1, color="black") +
  scale_fill_manual(values = c("pos" = "#1466dd", "neg" = "#c94c4c")) +
  geom_vline(xintercept = 0, size = 1.5, color = "black") +
  geom_hline(yintercept = 0, size = 1.5, color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none") +
  labs(y = "Gene count (total = 406)", x = "Chemical sensitivities across conditions (total = 46)",
       title = "Distribution of Chemical Sensitivites Per Gene") +
  geom_segment(aes(x = 31.5, y = -2, xend = 24.5, yend = -2), size = 2, color = "black") + # Bold line represents the lpt genes!
  coord_cartesian(ylim = c(-2, 25), xlim = c(0, 32))


#gene-level positive and negative fitness changes, plotted
unique_cond_counts_w_zeroes %>%
  group_by(locus_tag) %>%
  mutate(total_count = sum(count)) %>%
  filter(total_count !=0) %>%
  ungroup() %>%
  arrange(total_count, locus_tag) %>%
  mutate(direction = factor(direction, levels = c("pos", "neg"))) %>%
  ggplot(aes(x = reorder(locus_tag, total_count), y = count, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("neg" = "red", "pos" = "blue")) +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom") +
  labs(title = "Significant phenotypes for each gene",
       x = "A. baumannii essential gene",
       y = "Significant chemical phenotype",
       fill = "Fitness change direction")


#lpx vs lpt phenotypes (for extended view figure)
lpx_lpt <- unique_cond_counts_w_zeroes %>% filter(direction == "neg") %>% 
  filter(unique_name %in% c("lpxA", "lpxB", "lpxC", "lpxD", "lpxH", "lpxK", "lpxL", "lptA", "lptB", "lptC", "lptD", "lptF", "lptG")) %>%
  mutate(ID = if_else(startsWith(unique_name, "lpx"), "lpx",
                      if_else(startsWith(unique_name, "lpt"), "lpt", NA_character_)))

lpx_lpt %>%
  ggplot(aes(x=ID, y=count, fill = ID)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lpt" = "#b78395", "lpx" = "gray80")) +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none") +
  labs(y = "Chemical sensitivity count (total = 45)", x = "Gene group") +
  coord_cartesian(ylim = c(0, 35))

# Subset the counts for lpx
lpx_counts <- lpx_lpt$count[lpx_lpt$ID == "lpx"]

# Subset the counts for lpt
lpt_counts <- lpx_lpt$count[lpx_lpt$ID == "lpt"]

# Perform a t-test
t_test_result <- t.test(lpx_counts, lpt_counts)

# Print the results of the t-test
print(t_test_result)
