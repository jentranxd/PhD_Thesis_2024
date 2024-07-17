#########
###Mapping conservation (percent presence across representative species)###
###vs chemical response###


#load required packages
require("pacman")

p_load(data.table, dplyr, tidyr, ggplot2, ggrepel, ggforce)

#load in fraction/presence based on fDOG output and phenotype counts

phenotype_cons <- fread("./tables_phenotypes_conservation.csv", header = TRUE)
sig_phenotypes_updated <- fread("./sig_phenotypes_final.csv", header = TRUE)

# Remove old phenotype number column/match formatting
phenotype_cons_updated <- phenotype_cons %>% select(-c(`Phenotype count`))
colnames(phenotype_cons_updated)[colnames(phenotype_cons_updated) == 'Locus Tag'] <- "locus_tag"

#merge tables together
sig_pheno_final <- merge(phenotype_cons_updated, sig_phenotypes_updated, by = "locus_tag")

# Convert unique_condition_count to numeric for comparison and set as factor for later graphing
sig_pheno_final$unique_condition_count_numeric <- as.numeric(as.character(sig_pheno_final$unique_condition_count))
sig_pheno_final$unique_condition_count <- as.factor(sig_pheno_final$unique_condition_count)

# Create a placeholder row for unique_condition_count that are not represented 
#(33 here) with all columns set to NA
placeholder <- data.frame(matrix(NA, nrow = 1, ncol = ncol(sig_pheno_final)))
colnames(placeholder) <- colnames(sig_pheno_final)
placeholder$unique_condition_count_numeric <- 33
placeholder$unique_condition_count <- as.factor(33)

# Add the placeholder to the data
sig_pheno_final <- rbind(sig_pheno_final, placeholder)

# Define the correct order for the levels
levels_order <- sort(unique(as.numeric(as.character(sig_pheno_final$unique_condition_count))))
sig_pheno_final$unique_condition_count <- factor(sig_pheno_final$unique_condition_count, levels = levels_order)


# Create the boxplot for Non-Abau_Ac_Fraction
sig_pheno_final$label <- ifelse(sig_pheno_final$`Non-Abau_Ac_Fraction` < 0.4 & 
                                  sig_pheno_final$unique_condition_count_numeric >= 5, 
                                sig_pheno_final$unique_name, NA)

sig_pheno_final$is_labeled <- ifelse(is.na(sig_pheno_final$label), "No", "Yes")


ggplot(sig_pheno_final, aes(x = unique_condition_count, y = `Non-Abau_Ac_Fraction`, 
                            color = is_labeled)) +
  geom_sina(size = 2, alpha = 0.6) +  
  geom_text_repel(aes(label = label), size = 2, force = 0.2, min.segment.length = 0.3, nudge_x = -0.1, max.overlaps = 10) + 
  scale_color_manual(values = c("No" = "black", "Yes" = "#332288")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.05)) +
  labs(x = "Chemical-specific phenotypes", y = "Non-baumannii Acinetobacter Conservation") +
  ggtitle("Non-baumannii Acinetobacter Conservation vs. Significant Chemical Phenotypes") +
  theme_minimal() +
  theme(legend.position = "none")  

# Create the plot for Non-Ac_Moraxellales_Fraction
sig_pheno_final$label <- ifelse(sig_pheno_final$`Non-Ac_Moraxellales_Fraction` < 0.4 & 
                                  sig_pheno_final$unique_condition_count_numeric >= 5, 
                                sig_pheno_final$unique_name, NA)

sig_pheno_final$is_labeled <- ifelse(is.na(sig_pheno_final$label), "No", "Yes")


ggplot(sig_pheno_final, aes(x = unique_condition_count, y = `Non-Ac_Moraxellales_Fraction`, 
                            color = is_labeled)) +
  geom_sina(size = 2, alpha = 0.6) +  
  geom_text_repel(aes(label = label), size = 2, force = 0.2, min.segment.length = 0.3, nudge_x = -0.1, max.overlaps = 10) + 
  scale_color_manual(values = c("No" = "black", "Yes" = "#332288")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.05)) +
  labs(x = "Chemical-specific phenotypes", y = "Non-Acinetobacter Moraxellales Conservation") +
  ggtitle("Non-Acinetobacter Moraxellales Conservation vs. Significant Chemical Phenotypes") +
  theme_minimal() +
  theme(legend.position = "none")

# Create the plot for Non-Moraxellales_Gamma_Fraction
sig_pheno_final$label <- ifelse(sig_pheno_final$`Non-Moraxellales_Gamma_Fraction` < 0.4 & 
                                  sig_pheno_final$unique_condition_count_numeric >= 5, 
                                sig_pheno_final$unique_name, NA)

sig_pheno_final$is_labeled <- ifelse(is.na(sig_pheno_final$label), "No", "Yes")


ggplot(sig_pheno_final, aes(x = unique_condition_count, y = `Non-Moraxellales_Gamma_Fraction`, 
                            color = is_labeled)) +
  geom_sina(size = 2, alpha = 0.6) +  
  geom_text_repel(aes(label = label), size = 2, force = 0.2, min.segment.length = 0.3, nudge_y = -0.05, max.overlaps = 10) + 
  scale_color_manual(values = c("No" = "black", "Yes" = "#332288")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.05)) +
  labs(x = "Chemical-specific phenotypes", y = "Non-Moraxellales Gamma Conservation") +
  ggtitle("Non-Moraxellales Gammaproteobacteria Conservation vs. Significant Chemical Phenotypes") +
  theme_minimal() +
  theme(legend.position = "none")
