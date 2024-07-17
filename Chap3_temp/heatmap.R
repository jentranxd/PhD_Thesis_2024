library(heatmaply)
################
#set up requirements here

generate_breaks = function(x, n, center = F) {
  
  if (center) {
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)}
  
  else {
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)}
  
  return(res)
}

plot_colors <- c(
  colorRampPalette(
    c("firebrick", "white"))(144^3) %>%
    `[`((1:144)^3) %>%
    `[`(-144),
  "white",
  rev(colorRampPalette(
    c("royalblue4", "white"))(144^3) %>%
      `[`((1:144)^3) %>%
      `[`(-144)))

# # Define the custom color palette
# plot_colors <- colorRampPalette(c("firebrick", "white", "royalblue4"))(100)
# 
# # Define the breaks to ensure white is at 0
# breaks <- seq(-max_value, max_value, length.out = 101)

p_load(hrbrthemes, pheatmap)
#############################

#########################
#heatmap#

melted_results <- curated_names[melted_results, on = .(AB19606 == locus_tag)]


selected_LFC <- 
  melted_results[type == "perfect"]


#select appropriate conditions here if needed
selected_LFC <-
  melted_results[
    type == "perfect"
    & condition %in% c(
      'acriflavine_T2_7 - none_T2_5',
      'amikacin_T2_1 - none_T2_1a',
      'apramycin_T2_8 - none_T2_8',
      'azithromycin_T2_14 - none_T2_14',
      'aztreonam_T2_14 - none_T2_14',
      'carvacrol_T2_6 - none_T2_5',
      'cefaclor_T2_1 - none_T2_1a',
      'chlorhexidine_T2_14 - none_T2_14',
      'chlorhexidine_T2_3 - none_T2_3',
      'ciprofloxacin_T2_6 - none_T2_5',
      'colistin_T2_10 - none_T2_10',
      'copper_II_sulfate_T2_9 - none_T2_8',
      'daptomycin_T2_2 - none_T2_1b',
      'D_cycloserine_T2_3 - none_T2_3',
      # 'DMSO_T2_17 - none_T2_17',
      'doripenem_T2_9 - none_T2_8',
      'EDTA_T2_6 - none_T2_5',
      'ethidium_bromide_T2_6 - none_T2_5',
      'fosfomycin_T2_1 - none_T2_1a',
      'imipenem_T2_3 - none_T2_3',
      'indole_T2_14 - none_T2_14',
      'isoniazid_T2_14 - none_T2_14',
      'lactoferrin_T2_12 - none_T2_12',
      'levofloxacin_T2_6 - none_T2_5',
      'LL_37_T2_16 - none_T2_16',
      'lysozyme_T2_9 - none_T2_8',
      'mecillinam_T2_1 - none_T2_1a',
      'meropenem_T2_13 - none_T2_13',
      'minocycline_T2_15 - none_T2_15',
      'mupirocin_T2_4 - none_T2_3',
      'myricetin_T2_4 - none_T2_3',
      'Nickel_II_chloride_T2_10 - none_T2_10',
      # 'none_T2_10 - none_T0_10',
      # 'none_T2_11 - none_T0_11',
      # 'none_T2_12 - none_T0_12',
      # 'none_T2_13 - none_T0_13',
      # 'none_T2_14 - none_T0_14',
      # 'none_T2_15 - none_T0_15',
      # 'none_T2_16 - none_T0_16',
      # 'none_T2_17 - none_T0_17',
      # 'none_T2_1a - none_T0_1',
      # 'none_T2_1b - none_T0_2',
      # 'none_T2_3 - none_T0_3',
      # 'none_T2_5 - none_T0_5',
      # 'none_T2_8 - none_T0_8',
      'phenazine_T2_7 - none_T2_5',
      'polymyxin_b_T2_15 - none_T2_15',
      'rifampicin_T2_13 - none_T2_13',
      'SDS_T2_5 - none_T2_5',
      'tazobactam_piperacillin_T2_8 - none_T2_8',
      'tetracycline_T2_1 - none_T2_1a',
      'thymol_T2_6 - none_T2_5',
      'tigecycline_T2_2 - none_T2_1b',
      'tobramycin_T2_5 - none_T2_5',
      'TPEN_T2_5 - none_T2_5',
      'trimethoprim_T2_2 - none_T2_1b',
      'trimethoprim_sulfamethoxazole_T2_7 - none_T2_5',
      'vancomycin_T2_5 - none_T2_5',
      'Zinc_sulfate_T2_16 - none_T2_16'
    )]
###################################################
###################################################
median_LFC <- dcast(
  selected_LFC, 
  AB030 + unique_name ~ condition, 
  value.var = "LFC", 
  fun.aggregate = median)

###################################################
#rename with short names for antibiotics
median_LFC <- median_LFC %>% rename(
  'ACR' = 'acriflavine_T2_7 - none_T2_5',
  'AMK' = 'amikacin_T2_1 - none_T2_1a',
  'APR' = 'apramycin_T2_8 - none_T2_8',
  'AZI' = 'azithromycin_T2_14 - none_T2_14',
  'AZT' = 'aztreonam_T2_14 - none_T2_14',
  'CRV' = 'carvacrol_T2_6 - none_T2_5',
  'CEF' = 'cefaclor_T2_1 - none_T2_1a',
  'CHX (high)' = 'chlorhexidine_T2_14 - none_T2_14',
  'CHX (low)' = 'chlorhexidine_T2_3 - none_T2_3',
  'CIP' = 'ciprofloxacin_T2_6 - none_T2_5',
  'COL' = 'colistin_T2_10 - none_T2_10',
  'CuS' = 'copper_II_sulfate_T2_9 - none_T2_8',
  'DAP' = 'daptomycin_T2_2 - none_T2_1b',
  'DCS' = 'D_cycloserine_T2_3 - none_T2_3',
  'DRP' = 'doripenem_T2_9 - none_T2_8',
  'EDTA' = 'EDTA_T2_6 - none_T2_5',
  'EtBr' = 'ethidium_bromide_T2_6 - none_T2_5',
  'FOS' = 'fosfomycin_T2_1 - none_T2_1a',
  'IMP' = 'imipenem_T2_3 - none_T2_3',
  'IND' = 'indole_T2_14 - none_T2_14',
  'ISO' = 'isoniazid_T2_14 - none_T2_14',
  'LCF' = 'lactoferrin_T2_12 - none_T2_12',
  'LVF' = 'levofloxacin_T2_6 - none_T2_5',
  'LL37' = 'LL_37_T2_16 - none_T2_16',
  'LYS' = 'lysozyme_T2_9 - none_T2_8',
  'MEC' = 'mecillinam_T2_1 - none_T2_1a',
  'MRP' = 'meropenem_T2_13 - none_T2_13',
  'MIN' = 'minocycline_T2_15 - none_T2_15',
  'MUP' = 'mupirocin_T2_4 - none_T2_3',
  'MYR' = 'myricetin_T2_4 - none_T2_3',
  'NiCl' = 'Nickel_II_chloride_T2_10 - none_T2_10',
  'PHZ' = 'phenazine_T2_7 - none_T2_5',
  'PMB' = 'polymyxin_b_T2_15 - none_T2_15',
  'RIF' = 'rifampicin_T2_13 - none_T2_13',
  'SDS' = 'SDS_T2_5 - none_T2_5',
  'TZ-PIP' = 'tazobactam_piperacillin_T2_8 - none_T2_8',
  'TET' = 'tetracycline_T2_1 - none_T2_1a',
  'THY' = 'thymol_T2_6 - none_T2_5',
  'TIG' = 'tigecycline_T2_2 - none_T2_1b',
  'TOB' = 'tobramycin_T2_5 - none_T2_5',
  'TPEN' = 'TPEN_T2_5 - none_T2_5',
  'TMP' = 'trimethoprim_T2_2 - none_T2_1b',
  'TMP-SMX' = 'trimethoprim_sulfamethoxazole_T2_7 - none_T2_5',
  'VAN' = 'vancomycin_T2_5 - none_T2_5',
  'ZnSO4' = 'Zinc_sulfate_T2_16 - none_T2_16')

#####################
##alternatively, keep the long name antibiotics for the heatmap
# median_LFC <- median_LFC %>% rename(
#   'acriflavine' = 'acriflavine_T2_7 - none_T2_5',
#   'amikacin' = 'amikacin_T2_1 - none_T2_1a',
#   'apramycin' = 'apramycin_T2_8 - none_T2_8',
#   'azithromycin' = 'azithromycin_T2_14 - none_T2_14',
#   'aztreonam' = 'aztreonam_T2_14 - none_T2_14',
#   'carvacrol' = 'carvacrol_T2_6 - none_T2_5',
#   'cefaclor' = 'cefaclor_T2_1 - none_T2_1a',
#   'chlorhexidine (high)' = 'chlorhexidine_T2_14 - none_T2_14',
#   'chlorhexidine (low)' = 'chlorhexidine_T2_3 - none_T2_3',
#   'ciprofloxacin' = 'ciprofloxacin_T2_6 - none_T2_5',
#   'colistin' = 'colistin_T2_10 - none_T2_10',
#   'copper(II) sulfate' = 'copper_II_sulfate_T2_9 - none_T2_8',
#   'daptomycin' = 'daptomycin_T2_2 - none_T2_1b',
#   'D-cycloserine' = 'D_cycloserine_T2_3 - none_T2_3',
#   'DMSO' = 'DMSO_T2_17 - none_T2_17',
#   'doripenem' = 'doripenem_T2_9 - none_T2_8',
#   'EDTA' = 'EDTA_T2_6 - none_T2_5',
#   'ethidium bromide' = 'ethidium_bromide_T2_6 - none_T2_5',
#   'fosfomycin' = 'fosfomycin_T2_1 - none_T2_1a',
#   'imipenem' = 'imipenem_T2_3 - none_T2_3',
#   'indole' = 'indole_T2_14 - none_T2_14',
#   'isoniazid' = 'isoniazid_T2_14 - none_T2_14',
#   'lactoferrin' = 'lactoferrin_T2_12 - none_T2_12',
#   'levofloxacin' = 'levofloxacin_T2_6 - none_T2_5',
#   'LL-37' = 'LL_37_T2_16 - none_T2_16',
#   'lysozyme' = 'lysozyme_T2_9 - none_T2_8',
#   'mecillinam' = 'mecillinam_T2_1 - none_T2_1a',
#   'meropenem' = 'meropenem_T2_13 - none_T2_13',
#   'minocycline' = 'minocycline_T2_15 - none_T2_15',
#   'mupirocin' = 'mupirocin_T2_4 - none_T2_3',
#   'myricetin' = 'myricetin_T2_4 - none_T2_3',
#   'nickel(II) chloride' = 'Nickel_II_chloride_T2_10 - none_T2_10',
#   # 'none_T2_10 - none_T0_10',
#   # 'none_T2_11 - none_T0_11',
#   # 'none_T2_12 - none_T0_12',
#   # 'none_T2_13 - none_T0_13',
#   # 'none_T2_14 - none_T0_14',
#   # 'none_T2_15 - none_T0_15',
#   # 'none_T2_16 - none_T0_16',
#   # 'T2 no drug - T0' = 'none_T2_17 - none_T0_17',
#   # 'none_T2_1a - none_T0_1',
#   # 'none_T2_1b - none_T0_2',
#   # 'none_T2_3 - none_T0_3',
#   # 'none_T2_5 - none_T0_5',
#   # 'none_T2_8 - none_T0_8',
#   'phenazine' = 'phenazine_T2_7 - none_T2_5',
#   'polymyxin B' = 'polymyxin_b_T2_15 - none_T2_15',
#   'rifampicin' = 'rifampicin_T2_13 - none_T2_13',
#   'SDS' = 'SDS_T2_5 - none_T2_5',
#   'tazobactam-piperacillin' = 'tazobactam_piperacillin_T2_8 - none_T2_8',
#   'tetracycline' = 'tetracycline_T2_1 - none_T2_1a',
#   'thymol' = 'thymol_T2_6 - none_T2_5',
#   'tigecycline' = 'tigecycline_T2_2 - none_T2_1b',
#   'tobramycin' = 'tobramycin_T2_5 - none_T2_5',
#   'TPEN' = 'TPEN_T2_5 - none_T2_5',
#   'trimethoprim' = 'trimethoprim_T2_2 - none_T2_1b',
#   'trimethoprim-sulfamethoxazole' = 'trimethoprim_sulfamethoxazole_T2_7 - none_T2_5',
#   'vancomycin' = 'vancomycin_T2_5 - none_T2_5',
#   'zinc sulfate' = 'Zinc_sulfate_T2_16 - none_T2_16')
#   
#   

LFC_grid <- data.matrix(
  median_LFC[, c(-1:-2)])

rownames(LFC_grid) <- median_LFC$unique_name

# format = Pathway + italic(Gene Name)
newnames <- lapply(
  median_LFC[, .I],
  function(x) bquote(.(median_LFC[x]$Pathway) ~ italic(.(median_LFC[x]$unique_name))))

plot_matrix <- LFC_grid


#############select certain chemicals only
# selected_columns <- c("chlorhexidine (high)", "chlorhexidine (low)", "ciprofloxacin", "thymol", "EDTA")
# plot_matrix <- LFC_grid[, colnames(LFC_grid) %in% selected_columns]


##################select certain genes only
# gene_groups <- c("gtrOC1", "GO593_15125")
# # #
# filtered_rows <- grepl(paste(gene_groups, collapse="|"), rownames(LFC_grid))
# plot_matrix <- LFC_grid[filtered_rows, ]
# 


break_halves <- length(unique(as.vector(plot_matrix)))


breaks <- generate_breaks(plot_matrix, n = 286, center = T)


to_plot_title <- paste("Log Fold Change Fitness")
# 
# #######
# # #to transpose
# # plot_matrix <- t(plot_matrix)
# 
# ######
# 
to_plot <- pheatmap(plot_matrix,
                    col = plot_colors,
                    breaks = breaks,
                    border_color = NA,
                    cellwidth = 8,
                    cellheight = 8,
                    cutree_rows = 2,
                    # cutree_cols = 2,
                    main = to_plot_title,
                    angle_col = 90,
                    fontsize_col = 8,
                    fontsize_row = 8,
                    clustering_method = "ward.D2",
                    clustering_distance_rows = "canberra",
                    clustering_distance_cols = "canberra",
                    # labels_row = as.expression(newnames),
                    show_colnames = TRUE,
                    cluster_cols = TRUE
                    )

#####################
#####################

#select a section/cluster from the full heatmap

pheatmap_object <- pheatmap(plot_matrix, clustering_distance_rows = "canberra", clustering_method = "ward.D2", silent = TRUE)
row_dendrogram <- pheatmap_object$tree_row
original_column_order <- colnames(plot_matrix)


fraction_of_height = 0.2
h <- max(row_dendrogram$height) * fraction_of_height
clusters <- cutree(as.hclust(row_dendrogram), h = h)

# Find the cluster containing gene: GO593_10640 as an example
target_cluster <- clusters["GO593_10640"]

# Get all genes in that cluster
genes_in_target_cluster <- names(clusters[clusters == target_cluster])

sub_data <- plot_matrix[genes_in_target_cluster,]


to_plot <- pheatmap(sub_data,
                    col = plot_colors,
                    breaks = breaks,
                    border_color = NA,
                    cellwidth = 8,
                    cellheight = 8,
                    cutree_rows = 2,
                    # cutree_cols = 2,
                    main = to_plot_title,
                    angle_col = 90,
                    fontsize_col = 8,
                    fontsize_row = 8,
                    clustering_method = "ward.D2",
                    clustering_distance_rows = "canberra",
                    clustering_distance_cols = "canberra",
                    # labels_row = as.expression(newnames),
                    show_colnames = TRUE,
                    cluster_cols = TRUE
)

