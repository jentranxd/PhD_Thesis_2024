####Graphing membrane permeability graphs####

#load necessary packages
require("pacman")
p_load(ggplot2, tidyverse, data.table, dplyr)

#read in data; preprocessed in excel - fluorescence (time and values only, wells as columns) or OD (single timepoint, wide format)
cyo_etbr <- fread("../data/lpt cyo experiments/230901_EtBr_OxPhos_fluorescence.csv")
cyo_etbr_map <- fread("../data/lpt cyo experiments/230901_EtBr_well_map.csv")

#cyo ethidium bromide assay
# Convert fluorescence from wide to long format
cyo_etbr_long <- cyo_etbr %>% 
  pivot_longer(
    -`Time [s]`, # Exclude the first column (time) from melting
    names_to = "well", # Name of the new column for the 'melted' variable names
    values_to = "fluorescence" # Name of the new column for the 'melted' values
  )
# convert the time column to numeric if it's stored as characters
cyo_etbr_long$`Time [s]` <- as.numeric(cyo_etbr_long$`Time [s]`)

# Joining strain and rep columns based on 'well'
cyo_etbr_long_joined <- left_join(cyo_etbr_long, cyo_etbr_map, by = "well")

# Calculate the average fluorescence for blanks at each timepoint
average_blanks <- cyo_etbr_long_joined %>%
  filter(strain == "blank") %>%
  group_by(`Time [s]`) %>%
  summarize(average_fluorescence = mean(fluorescence))

# Adjust the fluorescence values by subtracting the average blank fluorescence
cyo_etbr_adjusted <- cyo_etbr_long_joined %>%
  left_join(average_blanks, by = "Time [s]") %>%
  mutate(adjusted_fluorescence = fluorescence - average_fluorescence) %>%
  select(`Time [s]`, well, strain, rep, adjusted_fluorescence) %>%
  filter(strain != "blank") #remove blanks

#convert ODs and normalize
# Set proper column names for clarity
colnames(cyo_etbr_OD) <- c("Row", paste0("Col", 2:7))

# Convert from wide to long format
cyo_etbr_OD_long <- cyo_etbr_OD %>%
  pivot_longer(cols = starts_with("Col"), names_to = "Column", values_to = "OD600") %>%
  mutate(well = paste0(Row, gsub("Col", "", Column))) %>%
  select(well, OD600)

# Joining strain and rep columns based on 'well'
cyo_etbr_OD_joined <- left_join(cyo_etbr_OD_long, cyo_etbr_map, by = "well")

# Calculate the average fluorescence for blanks at each timepoint
OD_blanks <- cyo_etbr_OD_joined %>%
  filter(strain == "blank")

average_OD600_blanks <- mean(OD_blanks$OD600)

# Adjust the fluorescence values by subtracting the average blank fluorescence
cyo_etbr_OD_adjusted <- cyo_etbr_OD_joined %>%
  mutate(adjusted_OD600 = OD600 - average_OD600_blanks) %>%
  select(well, strain, rep, adjusted_OD600) %>%
  filter(strain != "blank") #remove blanks

#normalize fluorescence to OD
normalized_cyo_etbr <- left_join(cyo_etbr_adjusted, cyo_etbr_OD_adjusted, by = c("well", "strain", "rep"))

normalized_cyo_etbr$normalized_fluorescence <- normalized_cyo_etbr$adjusted_fluorescence/normalized_cyo_etbr$adjusted_OD600

normalized_cyo_etbr$Time_min <- normalized_cyo_etbr$`Time [s]`/60

# Plotting
normalized_cyo_etbr %>% filter(strain %in% c("cyoA", "atpB", "nuoB", "sdhB", "nontargeting control")) %>% #select appropriate strains for figure
  ggplot(aes(x = Time_min, y = normalized_fluorescence, color = strain)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = strain), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("cyoA" = "navyblue", "atpB" = "green4", "nuoB" = "orange3", sdhB = "slateblue", "nontargeting control" = "gray20")) + 
  scale_fill_manual(values = c("cyoA" = "navyblue", "atpB" = "green4", "nuoB" = "orange3", sdhB = "slateblue", "nontargeting control" = "gray20")) +  
  labs(x = "Time (min)", y = "Fluorescence (AU)", 
       title = "EtBr permeability assay") +
  theme_minimal() 
