require('pacman')

p_load(clipr, data.table, dplyr, ggallin, scales, edgeR, statmod, poolr, pheatmap, purrr, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer, tidyr)


setwd("./appendix A")


mariner_tests <- fread("./191213_mariner_time_ratio.csv", header = TRUE)

# Convert ratio to factor for proper ordering in plot
mariner_tests$ratio <- factor(mariner_tests$`donor/recipient ratio`, levels = c(1, 0.333333333, 0.2))

# Create the plot
ggplot(mariner_tests, aes(x = ratio, y = Efficiency, color = factor(`Incubation time (hrs)`))) +
  geom_point(stat = "identity", position = position_dodge(), size = 4) +  # Adjust the size for bigger circles
  facet_wrap(~ Strain, scales = "free_y") +
  labs(x = "Donor/Recipient Ratio", y = "Efficiency", color = "Incubation Time (hrs)") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_y_log10(limits = c(1e-8, 1))

##########################################################
##########################################################
##########################################################

tn_tests <- fread("./200116_mariner_tn5.csv", header = TRUE)

# Convert ratio to factor for proper ordering in plot
tn_tests$tn <- factor(tn_tests$Transposon, levels = c('Tn5', 'mariner'))

# Create the plot
ggplot(tn_tests, aes(x = Strain, y = Efficiency, color = factor(tn))) +
  geom_point(stat = "identity", size = 4) +  # Adjust the size for bigger circles
  labs(x = "Strain", y = "Transposition Efficiency", color = "Transposon system") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_y_log10(limits=c(1e-7, 1))
