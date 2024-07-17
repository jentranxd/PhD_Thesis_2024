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

#load NAD+/NADH assay data

nad <- fread("C:/Users/jentr/OneDrive - UW-Madison/lab things/writing/dissertation/appendix c/230210_Abau_NAD_NADH.csv", header = TRUE)


# Subtract the corresponding replicate blanks from all the values
nad_adjusted <- nad %>%
  group_by(replicate) %>%
  mutate(lum_NAD = lum_NAD - lum_NAD[strain == "blank"],
         lum_NADH = lum_NADH - lum_NADH[strain == "blank"]) %>%
  ungroup()

# Calculate NAD/NADH ratio
nad_adjusted <- nad_adjusted %>%
  mutate(ratio = lum_NAD / lum_NADH)

# Set the order of strains
nad_adjusted$strain <- factor(nad_adjusted$strain, levels = c("non-targeting", "nuoB", "sdhB", "cyoA", "atpB"))

# Filter out the blank rows
nad_adjusted <- nad_adjusted %>%
  filter(strain != "blank")

# Plot the boxplots with replicates as dots
ggplot(nad_adjusted, aes(x = strain, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(color = "black", width = 0.2, size = 4, alpha = 0.7) +
  labs(title = "NAD+/NADH Ratio of Strains",
       x = "Strain",
       y = "NAD+/NADH Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))  
