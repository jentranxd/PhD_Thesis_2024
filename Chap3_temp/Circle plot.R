require(pacman)


##fluorescence graph over time/cycles

pacman::p_load(drc, hrbrthemes, purrr, lmtest, tidyverse, data.table, data.table, tidyverse, broom, modelr, Hmisc, R.utils, ggplot2,
               dplyr, packcircles, ggforce, ggrepel)



##############
##############
#chemical mode of action chart

# Create the data frame
data <- data.frame(
  antibiotic_class = c("cell wall/division", "translation", "DNA/replication", 
                       "folate metabolism", "membrane disruption", "transcription",
                       "tRNA ligase inhibition", "heavy metals/chelators",
                       "respiration/ROS", "poorly understood",
                       "protein denaturation"),
  count = c(12, 7, 3, 2, 6, 1, 1, 6, 1, 5, 1)
)

# Convert count to radii
data$radii <- sqrt(data$count / pi)

packing <- circleProgressiveLayout(data$radii, sizetype = 'radius')
data <- cbind(data, packing)

#plot


ggplot(data) +
  geom_circle(aes(x0 = x, y0 = y, r = radii, fill = antibiotic_class), color = NA) +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(aes(x = x, y = y, label = antibiotic_class),
                   box.padding = 0.35, 
                   point.padding = 1, 
                   segment.color = NA, 
                   force = 0.5,
                   max.overlaps = Inf,
                   max.iter = 2000) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


