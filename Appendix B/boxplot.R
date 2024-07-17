require("pacman")

p_load(ggplot2, dplyr)

#Load OD at saturation (media) data

sat_lib <- fread("./OD-saturation_library.csv", header = TRUE)

sat_lib$`OD600 at saturation` <- sat_lib$`OD600 at 20 hours (1:10 dilution)`*10


ggplot(sat_lib, aes(x = Media, y = `OD600 at saturation`)) +
  geom_boxplot() +
  geom_jitter(color = "black", width = 0.2, size = 4) +
  labs(title = "Boxplot of cell density at 20 hours",
       x = "Media",
       y = "OD600 at saturation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 4))
