setwd("./chemical clustering/Pathway effects limma-voom/")

source("enrichment_application.r")


library(dplyr)
library(ggplot2)
library(scales)
library(Hmisc)

score_set <- function(selected_contrast_pattern, pattern = TRUE) {

  if(pattern) {
    these_contrasts <- all_sets |>
      filter(contrast %like% selected_contrast_pattern) |>
      pull(contrast)
  }

  else(
    these_contrasts <- all_sets |>
      filter(contrast %in% selected_contrast_pattern) |>
      pull(contrast)
  )

  all_sets |>
    filter(contrast %like% selected_contrast_pattern) |> # could relace with c("contrast1", "contrast2")
    select(Direction, term, description, FDR, genes_targeted, gene_count, contrast) |>
    mutate(FDR = signif(FDR, 2)) |>
    mutate(prop_targ = signif(genes_targeted / gene_count, 2)) |>
    mutate(sign = ifelse(Direction == "Up", 1, -1)) |>
    mutate(score = sign * -log10(FDR * prop_targ) * sqrt(genes_targeted)) |>
    arrange(desc(abs(score)))
}

# to score everything you could do
score_set("*") |>
  mutate(score = signif(score, 2)) |>
  dcast(term + description + genes_targeted + gene_count ~ contrast, value.var = "score")
# then write to clipboard and paste into excel with clipr::write_clip()

create_enrichment_plot <- function(selected_term, selected_contrast_patterns) {
  # Filter 'term_stats' using the 'selected_term' argument
  title <- term_stats %>%
    filter(term == selected_term) %>%
    pull(description)
  title <- paste(title, " (", selected_term, ")", sep = "")
  
  # Helper function to check if any of the patterns match the contrast
  matches_any_pattern <- function(contrast, patterns) {
    any(sapply(patterns, function(pattern) grepl(pattern, contrast)))
  }
  
  # Construct the plot with dynamic filtering based on 'selected_term' and 'selected_contrast_patterns'
  enrichment_data <- contrast_assignments %>%
    inner_join(group_assignments) %>%
    inner_join(
      all_sets %>%
        filter(term == selected_term) %>%
        filter(sapply(contrast, matches_any_pattern, selected_contrast_patterns))
    ) %>%
        inner_join(contrast_assignments) %>%
        mutate(contrast = factor(contrast, levels = unique(contrast))) %>%
    mutate(
      label = case_when(
        FDR <= 0.05 ~ paste(
          contrast,
          paste0(Direction, " (", signif(FDR, 2), ")"),
          paste(paste(paste(genes_targeted, gene_count, sep = "/"), " genes present", sep = "")),
          sep = "\n"
      ), 
      FDR > 0.05 ~ paste(
        contrast,
        paste0("No change", " (", signif(FDR, 2), ")"),
        paste(paste(paste(genes_targeted, gene_count, sep = "/"), " genes present", sep = "")),
        sep = "\n"
      ))
    ) %>%
    mutate(label = factor(label, levels = unique(label))) %>%
    inner_join(enrichments) %>%
    inner_join(targets) %>%
    inner_join(annotated_data) %>%
    inner_join(v_targets) %>%
    arrange(assignment) %>%
    mutate(assignment = case_when(
      assignment == 1 ~ "Treatment",
      assignment == -1 ~ "Control"
    )) %>%
    mutate(group = factor(group, levels = unique(group))) %>%
    mutate(`Predicted Efficacy` = rescale(as.numeric(weight), to = c(1, 100)))

  quantiles <- enrichment_data %>% 
    group_by(contrast, assignment, label) %>% 
    summarise(cpm = wtd.quantile(cpm, weights = weight, probs = 0.5))
  
  print(quantiles)




  enrichment_plot <- enrichment_data %>%
    ggplot(aes(x = as.character(assignment), y = cpm)) +
    geom_tile(aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))), width = Inf, height = Inf, fill = "light grey") +
    geom_sina(aes(weight = weight, size = `Predicted Efficacy`, color = group), shape = 20, alpha = 0.5) +
    # geom_violin(aes(weight = weight), alpha = 0.0, draw_quantiles = c(0.25, 0.5, 0.75), lwd = 1.25) +
    geom_boxplot(aes(weight = weight), alpha = 0.0, lwd = 1.25) +
    scale_alpha_manual(values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(10^(0:5)), labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    facet_wrap(~label) +
    labs(x = NULL, y = "Counts per Million") +
    ggtitle(title) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    scale_size_continuous(range = c(0.5, 5), limits = c(0, 100)) +
    geom_label(data = quantiles, aes(label = round(cpm, 0)), alpha = 0.75) + # Add labels
    geom_label(data = quantiles, aes(label = round(cpm, 0)), fill = NA) + # Add labels

    theme_minimal()

  print(enrichment_plot)
}

# variety of things for the cyclines
create_enrichment_plot("CL:1517", "cycline") # heme copper terminal oxidase


enrichment_cheat_sheet <- curated_names %>% rename('locus_tag' = 'AB19606') %>% 
  select(c('locus_tag', 'unique_name')) %>% full_join(enrichments, by = 'locus_tag')
               
View(enrichment_cheat_sheet)

CL113 <- enrichment_cheat_sheet %>% filter(term == "CL:113") %>% select(unique_name)

# Filter contrasts that contain "- none_T2"
filtered_contrast_sets <- all_sets %>% 
  filter(grepl("- none_T2", contrast)) %>%
  filter(FDR <= 0.05)

# Calculate the percentage of unique terms for each contrast
termhits <- filtered_contrast_sets %>%
  group_by(contrast) %>%
  summarise(unique_terms_count = n_distinct(term)) %>%
  ungroup() %>%
  mutate(percentage = (unique_terms_count / n_distinct(all_sets$term)) * 100)

new_labels <- c(
  'acriflavine_T2_7 - none_T2_5' = 'ACR',
  'amikacin_T2_1 - none_T2_1a' = 'AMK',
  'apramycin_T2_8 - none_T2_8' = 'APR',
  'azithromycin_T2_14 - none_T2_14' = 'AZI',
  'aztreonam_T2_14 - none_T2_14' = 'AZT',
  'carvacrol_T2_6 - none_T2_5' = 'CRV',
  'cefaclor_T2_1 - none_T2_1a' = 'CEF',
  'chlorhexidine_T2_14 - none_T2_14' = 'CHX-hi',
  'chlorhexidine_T2_3 - none_T2_3' = 'CHX-lo',
  'ciprofloxacin_T2_6 - none_T2_5' = 'CIP',
  'colistin_T2_10 - none_T2_10' = 'COL',
  'copper_II_sulfate_T2_9 - none_T2_8' = 'CuS',
  'daptomycin_T2_2 - none_T2_1b' = 'DAP',
  'D_cycloserine_T2_3 - none_T2_3' = 'DCS',
  'DMSO_T2_17 - none_T2_17' = 'DMSO',
  'doripenem_T2_9 - none_T2_8' = 'DRP',
  'EDTA_T2_6 - none_T2_5' = 'EDTA',
  'ethidium_bromide_T2_6 - none_T2_5' = 'EtBr',
  'fosfomycin_T2_1 - none_T2_1a' = 'FOS',
  'imipenem_T2_3 - none_T2_3' = 'IMP',
  'indole_T2_14 - none_T2_14' = 'IND',
  'isoniazid_T2_14 - none_T2_14' = 'ISO',
  'lactoferrin_T2_12 - none_T2_12' = 'LCF',
  'levofloxacin_T2_6 - none_T2_5' = 'LVF',
  'LL_37_T2_16 - none_T2_16' = 'LL37',
  'lysozyme_T2_9 - none_T2_8' = 'LYS',
  'mecillinam_T2_1 - none_T2_1a' = 'MEC',
  'meropenem_T2_13 - none_T2_13' = 'MRP',
  'minocycline_T2_15 - none_T2_15' = 'MIN',
  'mupirocin_T2_4 - none_T2_3' = 'MUP',
  'myricetin_T2_4 - none_T2_3' = 'MYR',
  'Nickel_II_chloride_T2_10 - none_T2_10' = 'NiCl',
  'phenazine_T2_7 - none_T2_5' = 'PHZ',
  'polymyxin_b_T2_15 - none_T2_15' = 'PMB',
  'rifampicin_T2_13 - none_T2_13' = 'RIF',
  'SDS_T2_5 - none_T2_5' = 'SDS',
  'tazobactam_piperacillin_T2_8 - none_T2_8' = 'TZ-PIP',
  'tetracycline_T2_1 - none_T2_1a' = 'TET',
  'thymol_T2_6 - none_T2_5' =  'THY',
  'tigecycline_T2_2 - none_T2_1b' = 'TIG',
  'tobramycin_T2_5 - none_T2_5' = 'TOB',
  'TPEN_T2_5 - none_T2_5' = 'TPEN',
  'trimethoprim_T2_2 - none_T2_1b' = 'TMP',
  'trimethoprim_sulfamethoxazole_T2_7 - none_T2_5' = 'TMP-SMX',
  'vancomycin_T2_5 - none_T2_5' = 'VAN',
  'Zinc_sulfate_T2_16 - none_T2_16' = 'ZnSO4')

termhits$contrast <- new_labels[termhits$contrast]

# Plot the data
ggplot(termhits, aes(x = reorder(contrast, -percentage), y = percentage, fill = contrast)) +
  geom_bar(stat = "identity", aes(fill=percentage), color = "black") +
  scale_fill_gradient(low = "gray90", high = "black") +
  theme_minimal() +
  labs(title = "Percentage of Unique Terms in Each Contrast",
       x = "Contrast",
       y = "Percentage of Unique Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
