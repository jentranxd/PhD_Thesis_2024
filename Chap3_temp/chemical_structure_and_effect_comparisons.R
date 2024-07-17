#####################
#####################
# Dendrogram/clustering of chemical relationships/similarities 
# based on genetic pathway perturbations
#####################
#####################

#Load several packages from CRAN and Bioconductor
require("pacman")

p_load(
  data.table,
  scales,
  edgeR,
  statmod,
  poolr,
  pheatmap,
  svglite,
  ggplot2,
  ggrepel,
  Rtsne,
  pracma,
  colourpicker,
  RColorBrewer,
  vegan,
  tidyverse,
  magrittr,
  ggtext,
  ggforce,
  reshape2,
  dplyr,
  dendextend,
  rcdk,
  mclust
)

dge <- data_y #this input is direct from EdgeR script


#adjust column names
colnames(dge$design) <- gsub("[^[:alnum:]]", "_", colnames(dge$design))
colnames(dge$design) <- gsub("___", " - ", colnames(dge$design))

dge$samples$group <- gsub("[^[:alnum:]]", "_", dge$samples$group)
dge$samples$group <- gsub("___", " - ", dge$samples$group)

contrast_levels <- gsub("[^[:alnum:]]", "_", contrast_levels)
contrast_levels <- gsub("___", " - ", contrast_levels)

contrast_list <- setNames(as.list(contrast_levels), contrast_levels)

# Add the 'levels' argument to the list
contrast_list$levels <- dge$design

# Then use these contrast levels in the makeContrasts function
contrasts <- do.call(makeContrasts, contrast_list)

#read in STRING groups
all_string <- fread("STRG0060QIE.protein.enrichment.terms.v11.5.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()

targets <- fread("Ab_library.tsv", na.strings = "None")


gene_groups <- all_string %>%
  # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
  group_by(category, term, description) %>%
  summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
  mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))


term_stats <- gene_groups %>%
  unnest(locus_tag) %>%
  inner_join(targets %>%
               select(locus_tag) %>% unique()) %>%
  group_by(term, gene_count, description) %>%
  summarize(genes_targeted = n())

complete_terms <- term_stats %>%
  filter(gene_count == genes_targeted)

# only perform enrichments where all genes are available
# gene_groups <- complete_terms %>% inner_join(gene_groups)

repeated_gene_groups <- gene_groups %>%
  group_by(locus_tag) %>%
  mutate(times_listed = n()) %>%
  arrange(locus_tag) %>%
  ungroup()


# pick the best annotation for each locus_tag_group, i.e., highest in term, and the lowest in the category_rank
ranked_annotations <- repeated_gene_groups %>%
  group_by(locus_tag_group, category) %>%
  arrange(versionsort::ver_sort(term)) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(category_rank = case_when(
    category == "Biological Process (Gene Ontology)" ~ 1,
    category == "Molecular Function (Gene Ontology)" ~ 2,
    category == "Cellular Component (Gene Ontology)" ~ 3,
    category == "Protein Domains and Features (InterPro)" ~ 4,
    category == "Protein Domains (SMART)" ~ 5,
    category == "Protein Domains (Pfam)" ~ 6,
    category == "Annotated Keywords (UniProt)" ~ 7,
    category == "Reactome Pathways" ~ 8,
    category == "Subcellular localization (COMPARTMENTS)" ~ 9,
    category == "Local Network Cluster (STRING)" ~ 10,
    TRUE ~ NA_integer_
  )) %>%
  group_by(locus_tag_group) %>%
  filter(category_rank == min(category_rank))

enrichments <- ranked_annotations %>%
  ungroup() %>%
  distinct(locus_tag_group, .keep_all = TRUE) %>%
  select(-locus_tag_group) %>%
  unnest(locus_tag) %>%
  inner_join(term_stats)


# Get the unique terms
unique_terms <- unique(enrichments$term)

target_spacers_for_terms <- term_stats %>%
  inner_join(enrichments, relationship = "many-to-many") %>%
  inner_join(targets, relationship = "many-to-many")


#########################################################################################

most_representative_sets <- enrichments %>%
  inner_join(targets %>% select(spacer, locus_tag)) %>%
  group_by(term) %>%
  arrange(spacer) %>%
  mutate(spacers = paste(spacer, collapse = ",")) %>%
  select(spacers, term, genes_targeted) %>%
  unique() %>%
  group_by(spacers) %>%
  inner_join(enrichments %>% select(term, gene_count)) %>%
  unique() %>%
  mutate(pct_targ = genes_targeted / gene_count) %>%
  ungroup() %>%
  group_by(spacers) %>%
  filter(pct_targ == max(pct_targ)) %>%
  ungroup() %>%
  pull(term)

#########################################################################################

# Split the spacer column by term
sets_to_locus_tags <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
sets_to_locus_tags_indices <- lapply(sets_to_locus_tags, function(locus_tags) which(rownames(dge) %in% locus_tags))

v <- voomWithQualityWeights(dge, dge$design, plot = TRUE)

# filter out guides that do not have bona-fide targets, but keep non-targeting guides
# i.e. guides that have no target, or guides that have a target in the same (wrong) direction
# overlap of spacer to gene should be maximal, i.e. 20

v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>%
      filter(locus_tag %in% all_string$locus_tag) %>%
      group_by(spacer) %>%
      filter(
        is.na(target) |
          target == "None" |
          (
            sp_dir != tar_dir &
              abs(as.numeric(offset)) == min(abs(as.numeric(offset))) &
              overlap == max(overlap)
          )
      ) %>%
      group_by(target)
  )

# Assign the weight to the guides based on y_pred to be between 1 and 100
v_targets[y_pred == "None", y_pred := NA_integer_]

v_targets$y_pred <- as.numeric(v_targets$y_pred)

v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

v_targets[mismatches >= 1, weight := y_pred]

v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))

# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = sets_to_locus_tags_indices,
    design = dge$design,
    weights = v_targets$weight,
    inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)

all_sets <- all_sets %>%
  inner_join(enrichments) %>%
  inner_join(v_targets) %>%
  inner_join(term_stats) %>%
  group_by(contrast, term, description) %>%
  nest(locus_tags = locus_tag) %>%
  group_by(locus_tags, contrast) %>%
  mutate(missing_genes = gene_count - genes_targeted) %>%
  arrange(FDR, missing_genes) %>%
  ungroup() %>%
  rename(guide_count = NGenes) %>%
  select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count) %>%
  unique() %>%
  data.table()

contrast_assignments <- contrasts %>%
  data.table(keep.rownames = "group") %>%
  melt(
    id.vars = "group",
    variable.name = "contrast",
    value.name = "assignment"
  ) %>%
  filter(assignment != 0)

group_assignments <- dge$design %>%
  data.table(keep.rownames = "sample") %>%
  melt(id.vars = "sample", variable.name = "group") %>%
  filter(value != 0) %>%
  select(-value)

original_data <- dge$counts %>%
  data.table(keep.rownames = "spacer") %>%
  melt(
    value.name = "count",
    id.vars = "spacer",
    variable.name = "sample"
  )

annotated_data <- dge$samples %>%
  data.table(keep.rownames = "sample") %>%
  inner_join(original_data) %>%
  group_by(sample) %>%
  mutate(cpm = 1e6 * count / sum(count))


scores_wide <- all_sets %>%
  mutate(score = case_when(Direction == "Up" ~ -log10(FDR), Direction == "Down" ~ log10(FDR))) %>%
  dcast(term + description + gene_count + genes_targeted + guide_count ~ contrast, value.var = "score")


all_sets %>%
  mutate(score = case_when(Direction == "Up" ~ -log10(FDR), Direction == "Down" ~ log10(FDR))) %>%
  dcast(term + description + gene_count + genes_targeted + guide_count ~ contrast, value.var = "score") %>%
  inner_join(target_spacers_for_terms %>% select(term, spacer))


###############
# Clustering
repreSetsAndScores <- all_sets %>%
  mutate(score = case_when(Direction == "Up" ~ -log10(FDR), Direction == "Down" ~ log10(FDR))) %>%
  dcast(term + description + gene_count + genes_targeted + guide_count ~ contrast,
        value.var =
          "score"
  ) %>%
  filter(term %in% most_representative_sets)

repreSetsAndScoresMatrix <- repreSetsAndScores %>%
  select(
    -term,
    -description,
    -gene_count,
    -genes_targeted,
    -guide_count
  ) %>%
  data.matrix()

rownames(repreSetsAndScoresMatrix) <- repreSetsAndScores$term

pathwayheatmap <- repreSetsAndScoresMatrix %>% 
  pheatmap()

# Clustering
# Pull the row dendrogram from the heatmap

row_dend <- pathwayheatmap$tree_row
cutree(row_dend, k = 10) %>%
  as.list() %>%
  as_tibble() %>%
  data.table() %>%
  melt(variable.name = "term", value.name = "group") %>%
  inner_join(repreSetsAndScores) %>%
  filter(group == 10)


# Define the color palette
my_palette <- colorRampPalette(c("red", "white", "blue"))

# Define the number of colors you want in your palette
num_colors <- 10001

# Generate the color gradient
color_gradient <- my_palette(num_colors)

color_gradient[num_colors] <- "white"

# Define the breaks to correspond to the colors
breaks <- seq(-1, 1, length.out = num_colors)

repreSetsAndScoresMatrixCor <- repreSetsAndScoresMatrix %>%
  cor()

diag(repreSetsAndScoresMatrixCor) <- 0

# Use the color gradient in the pheatmap function
geneSetHeatMap <- repreSetsAndScoresMatrixCor %>%
  pheatmap(
    clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = color_gradient, # Use the color gradient here
    breaks = breaks # Use the breaks here
  )

plot(pathwayheatmap$tree_col)

# Assuming pathwayheatmap$tree_col is a hclust
hclust <- pathwayheatmap$tree_col

# Convert hclust to dendrogram
dend <- as.dendrogram(hclust)


#plot without colors
plot(dend, main = "Chemical dendrogram based on effect on pathways, ward.d2)") # 'label = FALSE' if you don't want to label each leaf


# Adjust the size, add colors, and modify line properties in a single chain of operations
hc <- dend %>%
  set("labels_cex", 1) %>% # Change the font size
  color_branches(k = 5) %>% # Color branches based on clusters
  set("branches_lwd", 3) %>% # Line width
  # set("branches_lty", 2) %>% # Uncomment if you want to change line type
  color_labels(k = 5) # Color labels based on clusters

# Plot
plot(hc, dend_track_height = 0.5)

#


##################
##################
#clustering of chemicals based on fingerprinting 
##################
##################

# Read your CSV file, which contains your chemical names and SMILES from PubChem
# For drug combinations, like trimethoprim-sulfamethoxazole or tazo-piperacillin -  
# I used the SMILES for the main active ingredient, based on experimental growth defect with individual components
# For inherently mixed antibiotics (colistin - a mix of multiple polymyxin Es) - 
# I selected the SMILES string for one; these are chemically very similar so this method will work for this case
# Other approaches may be needed for different chemical combinations

chem_data <- read.csv("./antibiotics_in_samples_with_SMILES.csv", stringsAsFactors = FALSE, header = TRUE)
chem_data <- chem_data %>% filter(shortform.chemical.name != "DMSO")
# Assuming your CSV has a column named 'SMILES' with SMILES strings
smiles_vector <- chem_data$SMILES
names(smiles_vector) <- chem_data$comparison
molecules <- parse.smiles(smiles_vector)

#generate molecular fingerprints - unique binary identifiers based on chemical composition
fps <- lapply(molecules, get.fingerprint, type="standard")


# Convert an S4 fingerprint object to a numeric vector of bits on
convert_fp_to_vector <- function(fp) {
  bits_on <- fp@bits # 'bits' is a slot, representing positions where bits are on
  return(bits_on)
}

# Assuming 'fps' is your list of S4 fingerprint objects
fps_vectors <- lapply(fps, convert_fp_to_vector)

# Now 'fps_vectors' is a list of numeric vectors

# Function to calculate Tanimoto coefficient - measure of chemical similarity
tanimoto <- function(fp1, fp2) {
  intersection = length(intersect(fp1, fp2))
  union = length(union(fp1, fp2))
  
  return(intersection / union)
}

# Initialize an empty matrix to store the Tanimoto coefficients
num_fps <- length(fps_vectors)
similarity_matrix <- matrix(nrow = num_fps, ncol = num_fps, dimnames = list(names(fps_vectors), names(fps_vectors)))

# Calculate Tanimoto coefficients for all pairs
for (i in 1:num_fps) {
  for (j in i:num_fps) { 
      similarity_matrix[i, j] <- tanimoto(fps_vectors[[i]], fps_vectors[[j]])
      similarity_matrix[j, i] <- similarity_matrix[i, j] # Fill in the symmetric value
    }
  }

# `similarity_matrix` now contains the Tanimoto coefficients for all pairs
plot_sim_matrix <- similarity_matrix
diag(plot_sim_matrix) <- NA
color_palette <- colorRampPalette(c("white", "blue"))(100)
# 
# # Plotting
# pheatmap(plot_sim_matrix,
#          color = color_palette,
#          na_col = "gray", # This will set the color for NA values
#          clustering_distance_rows = "correlation",
#          clustering_distance_cols = "correlation",
#          clustering_method = "ward.D2")


###########
# Using Tanimoto coefficients for hierarchichal trees
tanimoto_distance_matrix <- 1 - similarity_matrix

# Convert the distance matrix to a dist object if necessary
distance_object <- as.dist(tanimoto_distance_matrix)

# Perform hierarchical clustering
hc_structure <- hclust(distance_object, method = "ward.D2") %>% as.dendrogram()


# Plot the dendrogram
plot(hc_structure, main = "Cluster dendrogram (1 - Tanimoto coefficient, ward.d2)") # 'label = FALSE' if you don't want to label each leaf



##################
##################
#Compare the pathway level to the tanimoto/structural clustering of the chemicals
##################
##################

#Bring in the genetic pathway distance values
pathway_dend <- dend

#Bring in the chemical structure dendrogram
structure_dend <- hc_structure

#plot side-by-side
par(mfrow=c(1,2))
plot(pathway_dend, main="Cluster dendrogram (pathway-level effects)")
plot(structure_dend, main="Cluster dendrogram (chemical structure)")

#side-by-side plot with color branches for visual effect
k <- 5
dendp_colored <- color_branches(pathway_dend, k=k)
dends_colored <- color_branches(structure_dend, k=k)

par(mfrow=c(1,2))
plot(dendp_colored, main="Cluster dendrogram (pathway-level effects)")
plot(dends_colored, main="Cluster dendrogram (chemical structure)")
################

#compare the two dendrograms


# Modify labels

# abbrevations
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

## full names
# new_labels <- c(
#   'acriflavine_T2_7 - none_T2_5' = 'acriflavine',
#   'amikacin_T2_1 - none_T2_1a' = 'amikacin',
#   'apramycin_T2_8 - none_T2_8' = 'apramycin',
#   'azithromycin_T2_14 - none_T2_14' = 'azithromycin',
#   'aztreonam_T2_14 - none_T2_14' = 'aztreonam',
#   'carvacrol_T2_6 - none_T2_5' = 'carvacrol',
#   'cefaclor_T2_1 - none_T2_1a' = 'cefaclor',
#   'chlorhexidine_T2_14 - none_T2_14' = 'chlorhexidine (high)',
#   'chlorhexidine_T2_3 - none_T2_3' =  'chlorhexidine (low)',
#   'ciprofloxacin_T2_6 - none_T2_5' =  'ciprofloxacin',
#   'colistin_T2_10 - none_T2_10' = 'colistin',
#   'copper_II_sulfate_T2_9 - none_T2_8' =  'copper (II) sulfate',
#   'daptomycin_T2_2 - none_T2_1b' = 'daptomycin',
#   'D_cycloserine_T2_3 - none_T2_3' = 'D-cycloserine',
#   'DMSO_T2_17 - none_T2_17' = 'DMSO',
#   'doripenem_T2_9 - none_T2_8' = 'doripenem',
#   'EDTA_T2_6 - none_T2_5' = 'EDTA',
#   'ethidium_bromide_T2_6 - none_T2_5' = 'ethidium bromide',
#   'fosfomycin_T2_1 - none_T2_1a' = 'fosfomycin',
#   'imipenem_T2_3 - none_T2_3' = 'imipenem',
#   'indole_T2_14 - none_T2_14' = 'indole',
#   'isoniazid_T2_14 - none_T2_14' = 'isoniazid',
#   'lactoferrin_T2_12 - none_T2_12' = 'lactoferrin',
#   'levofloxacin_T2_6 - none_T2_5' = 'levofloxacin',
#   'LL_37_T2_16 - none_T2_16' = 'LL-37',
#   'lysozyme_T2_9 - none_T2_8' = 'lysozyme',
#   'mecillinam_T2_1 - none_T2_1a' = 'mecillinam',
#   'meropenem_T2_13 - none_T2_13' = 'meropenem',
#   'minocycline_T2_15 - none_T2_15' = 'minocycline',
#   'mupirocin_T2_4 - none_T2_3' = 'mupirocin',
#   'myricetin_T2_4 - none_T2_3' = 'myricetin',
#   'Nickel_II_chloride_T2_10 - none_T2_10' = 'nickel (II) chloride',
#   'phenazine_T2_7 - none_T2_5' =  'phenazine',
#   'polymyxin_b_T2_15 - none_T2_15' = 'polymyxin B',
#   'rifampicin_T2_13 - none_T2_13' = 'rifampicin',
#   'SDS_T2_5 - none_T2_5' = 'SDS',
#   'tazobactam_piperacillin_T2_8 - none_T2_8' = 'tazobactam-piperacillin',
#   'tetracycline_T2_1 - none_T2_1a' = 'tetracycline',
#   'thymol_T2_6 - none_T2_5' = 'thymol',
#   'tigecycline_T2_2 - none_T2_1b' = 'tigecycline',
#   'tobramycin_T2_5 - none_T2_5' = 'tobramycin',
#   'TPEN_T2_5 - none_T2_5' = 'TPEN',
#   'trimethoprim_T2_2 - none_T2_1b' = 'trimethoprim',
#   'trimethoprim_sulfamethoxazole_T2_7 - none_T2_5' = 'trimethoprim-sulfamethoxazole',
#   'vancomycin_T2_5 - none_T2_5' = 'vancomycin',
#   'Zinc_sulfate_T2_16 - none_T2_16' = 'zinc sulfate')

# Function to apply label changes with debug output
apply_labels_debug <- function(dend, label_map) {
  if (is.leaf(dend)) {
    old_label <- attr(dend, "label")
    new_label <- label_map[old_label]
    if (!is.null(new_label)) {
      cat(sprintf("Replacing %s with %s\n", old_label, new_label)) # Debug message
      attr(dend, "label") <- new_label
    } else {
      cat(sprintf("No replacement found for %s\n", old_label)) # Debug message for no match
    }
  } else {
    dend[[1]] <- apply_labels_debug(dend[[1]], label_map)
    dend[[2]] <- apply_labels_debug(dend[[2]], label_map)
  }
  return(dend)
}

# Apply the modified function with debug messages
pathway_dend <- apply_labels_debug(pathway_dend, new_labels)
structure_dend <- apply_labels_debug(structure_dend, new_labels)


##########################
# score chemical pairs based on clustering within dendrograms and across dendrograms

# define a range of k values for clustering 
k_values <- seq(1, 45) #this is for min to max for these dendrograms

leaf_names <- labels(structure_dend)  # Reference leaf names from the pathway dendrogram

# Number of leaves
n_leaves <- length(labels(pathway_dend))  # Or use `structure_dend`; they should be the same

# Initialize an empty list to store scores for each k
scores_list <- vector("list", length(k_values))

# Assume leaf_names contains all unique leaf names from both dendrograms and is sorted if necessary

for (ki in seq_along(k_values)) {
  k <- k_values[ki]
  # Cut the dendrograms at this k to obtain cluster labels
  labels_pathway <- cutree(pathway_dend, k = k)
  labels_structure <- cutree(structure_dend, k = k)
  
  # Align the labels based on leaf_names
  # This assumes labels_pathway and labels_structure are named vectors where names correspond to leaf names
  aligned_labels_pathway <- labels_pathway[leaf_names]
  aligned_labels_structure <- labels_structure[leaf_names]
  
  # Initialize a temporary score matrix for this k
  temp_score_matrix <- matrix(0, nrow = n_leaves, ncol = n_leaves, dimnames = list(leaf_names, leaf_names))
  
  # Compare and update the score matrix based on aligned labels
  # If 2 chemicals are in the same cluster, add 1 
  # chemicals in the same cluster for both dendrograms get double the score
  for (i in 1:n_leaves) {
    for (j in 1:n_leaves) {
      if (aligned_labels_pathway[i] == aligned_labels_pathway[j]) {
        temp_score_matrix[i, j] <- temp_score_matrix[i, j] + 1
      }
      if (aligned_labels_structure[i] == aligned_labels_structure[j]) {
        temp_score_matrix[i, j] <- temp_score_matrix[i, j] + 1
      }
    }
  }
  
  # Store the temp_score_matrix for this k in the list
  scores_list[[ki]] <- temp_score_matrix
}


# Aggregate scores across all k values to produce a final agreement score matrix
# This involves summing up the scores from each k and then normalizing
agreement_matrix <- Reduce("+", scores_list) / (2 * length(k_values)) # Divide by 2*length(k_values) for normalization

#plot
pheatmap(agreement_matrix,
         color = colorRampPalette(c("white", "#332288"))(100),
         border_color = NA,
         main = "Cluster similarity within and across dendrograms",
         legend_title = "Agreement Score",
         cluster_rows = FALSE,
         cluster_cols = FALSE)

# Convert the agreement matrix to a dataframe
agreement_long <- as.data.frame(as.table(agreement_matrix))

# Rename the columns for clarity
names(agreement_long) <- c("Leaf1", "Leaf2", "Value")

# Convert factors to characters to avoid comparison issues
agreement_long$Leaf1 <- as.character(agreement_long$Leaf1)
agreement_long$Leaf2 <- as.character(agreement_long$Leaf2)

# Remove duplicates (keep pairs where Leaf1 is alphabetically before Leaf2)
# This step also removes rows where Leaf1 == Leaf2, as they are not needed for comparison
agreement_long <- agreement_long[which(agreement_long$Leaf1 < agreement_long$Leaf2), ]

# Sort the dataframe by the Value column in descending order
agreement_long <- agreement_long[order(-agreement_long$Value), ]

print(agreement_long)

############################
# In a similar manner, score within clusters only

# Initialize empty matrices for associations
association_matrix_pathway <- matrix(0, nrow = n_leaves, ncol = n_leaves, dimnames = list(leaf_names, leaf_names))
association_matrix_structure <- matrix(0, nrow = n_leaves, ncol = n_leaves, dimnames = list(leaf_names, leaf_names))

for (k in k_values) {
  # Cut the dendrograms into clusters
  labels_pathway_raw <- cutree(pathway_dend, k)
  labels_structure_raw <- cutree(structure_dend, k)
  
  # Align labels according to leaf_names
  # Convert to named vectors for alignment
  labels_pathway_aligned <- labels_pathway_raw[match(leaf_names, names(labels_pathway_raw))]
  labels_structure_aligned <- labels_structure_raw[match(leaf_names, names(labels_structure_raw))]
  
  # Update association matrices
  for (i in 1:(n_leaves-1)) {
    for (j in (i+1):n_leaves) {
      if (labels_pathway_aligned[i] == labels_pathway_aligned[j]) {
        association_matrix_pathway[i, j] <- association_matrix_pathway[i, j] + 1
        association_matrix_pathway[j, i] <- association_matrix_pathway[j, i] + 1
      }
      if (labels_structure_aligned[i] == labels_structure_aligned[j]) {
        association_matrix_structure[i, j] <- association_matrix_structure[i, j] + 1
        association_matrix_structure[j, i] <- association_matrix_structure[j, i] + 1
      }
    }
  }
}

# Normalize the association matrices by the number of k-values considered
association_matrix_pathway <- association_matrix_pathway / length(k_values)
association_matrix_structure <- association_matrix_structure / length(k_values)

# Calculate the difference in associations between the two dendrograms
# This will highlight if a chemical pair has a high association in one dendrogram
# but not in the other
association_difference <- abs(association_matrix_pathway - association_matrix_structure)

# plot
pheatmap(association_difference,
         color = colorRampPalette(c("white", "darkred"))(100),
         border_color = "lightgrey",
         main = "Association differences between dendrograms",
         legend_title = "Association score difference",
         cluster_rows = FALSE,
         cluster_cols = FALSE)

association_difference_bidirectional <- association_matrix_structure - association_matrix_pathway


# Define a custom color palette that transitions 
custom_colors <- colorRampPalette(c("#332288", "white", "#661100"))(100)
# custom_colors <- viridis::cividis(100)

# Define fixed breaks from -1 to 1
breaks <- seq(-1, 1, length.out = 101)

# Plotting the heatmap with the custom color palette and fixed breaks
pheatmap(association_difference_bidirectional,
         color = custom_colors,  # Use the custom color palette
         breaks = breaks,  # Use fixed breaks from -1 to 1
         border_color = NA,
         main = "Structure minus pathway",
         legend_title = "Association Score Difference",
         cluster_rows = FALSE,
         cluster_cols = FALSE)




