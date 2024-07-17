require("pacman")

p_load(ggplot2, dplyr)


#load ATP assay data

atp <- fread("./atp_assay.csv", header = TRUE)

# Calculate relative luminescence
atp <- atp %>%
  group_by(strain) %>%
  mutate(relative_luminescence = luminescence / OD_spec) %>%
  ungroup() %>%
  mutate(relative_luminescence = relative_luminescence / 
           mean(relative_luminescence[strain == "non-targeting"]))

# Set the order of strains
atp$strain <- factor(atp$strain, levels = c("non-targeting", "nuoB", "sdhB", "cyoA", "atpB"))

ggplot(atp, aes(x = strain, y = relative_luminescence)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color = "black", width = 0.2, size = 4, alpha = 0.7) +
  labs(title = "ATP assay",
       x = "Strain",
       y = "Relative ATP (normalized luminescence)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

####################

