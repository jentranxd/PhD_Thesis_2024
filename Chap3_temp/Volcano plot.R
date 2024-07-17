if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(ggplot2)
library(data.table)
library(gridExtra)

# Convert your data.table to a data.frame if it's not already one, as ggplot generally works with data.frames
melted_results_plot <- as.data.frame(melted_results)

# Apply the new significance criteria
melted_results_plot$significant <- ifelse(abs(melted_results_plot$LFC) > 1 & melted_results_plot$FDR < 0.05, "Significant", "Not Significant")

# Convert FDR to -log10 scale for plotting
melted_results_plot$negLogFDR <- -log10(melted_results_plot$FDR)

melted_results_plot <- merge(melted_results_plot, curated_names, by.x = "locus_tag", by.y = "AB19606")

melted_results_plot %>% filter(condition %like% "condition") %>% #filter specific condition
  filter(type == "perfect") %>% #filter guide type
  ggplot(aes(x = LFC, y = negLogFDR, color = significant)) +
  geom_point(alpha = 0.5) +  # Adjust point transparency as needed
  theme_minimal() +  # Use a clean theme
  scale_color_manual(values = c("Significant" = "black", "Not Significant" = "grey")) +  # Color code
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  ggrepel::geom_label_repel(aes(label=unique_name)) +
  theme(plot.title = element_text(hjust = 0.5))  # Center the plot title
