# Load necessary packages
# The pacman package will help us install and load packages seamlessly
install.packages('pacman', dependencies = TRUE)
library('pacman')

# Use p_load to install and load multiple packages in one step
p_load(clipr, data.table, dplyr, ggallin, scales, edgeR, statmod, poolr, pheatmap, 
       purrr, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer, tidyr)

# Load Gene Names and Annotations
# Curated names mapping gene numbers to gene names (User-provided)
curated_names <- fread("curated_names.tsv")

# aba_bed contains gene location information on the chromosome (User-provided)
aba_bed <- fread("CP046654.1.bed",
                 col.names = c("chromosome", "left", "right", "locus_tag", "gene_name",
                               "strand", "coding", "completeness"))

# aba_key maps guides to genes (User-provided)
aba_key <- fread("aba_key.tsv")

# Load Data
# Read in count data (User-provided)
data <- fread("count_by_position.tsv.gz",
              header = TRUE,
              col.names = c("condition", "spacer", "count")) %>%
  select(spacer, count, condition)

# Standardize Sample Names
# Harmonize condition names to a common format
data$condition <- dplyr::recode(data$condition, 
                                # Your specific recode rules here
)

# Remove unrelated lab samples (if any)
data <- data[!data$condition %like% "dJMP",]

# Load Experimental Design
# Read in custom experimental design details (User-provided)
data_design <- fread("all_experimental_design.tsv", na.strings = c("NA"))

# Rename column to match 'data'
data_design <- rename(data_design, condition = sample)

# Sort data and design by 'condition'
data.table::setorder(data, condition)
data.table::setorder(data_design, condition)

# Limit the design space to only 'tube' experiments
data_design <- data_design[data_design$experiment == "tube"]

# Quality Control and Filtering
# Load lists of samples that fail QC (User-provided files)
botneck <- fread("conditions_list_botneck.tsv", header = TRUE)
corrQC <- fread("condition_list_correlation_qc.tsv")

# Filter out samples that failed QC, except for the specific control 4AbJT67
data_design <- data_design %>% 
  filter(
    condition %in% corrQC$sample1 & 
      !(timing == "T2" & !(condition %in% botneck$condition))
  )

# Limit to counts that exist within the defined experimental design space
data <- data[condition %in% data_design$condition]

# Data Transformation and Integrity Check
# Reshape data from long format to wide format
data_grid <- data.table::dcast(
  data, 
  spacer ~ condition,
  value.var = "count", 
  fill = 0
)

# Melting data back to long format for integrity check
data_grid_remelted <- melt(data_grid, variable.name = "condition", value.name = "count", id.vars = c('spacer'))
data.table::setorder(data_grid_remelted, condition)

# Creating data matrix from the grid
data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])
row.names(data_grid_matrix) <- data_grid$spacer

# Defining Experimental Groups
# Create a factor combining drug, timing, and batch information
data_group <- factor(data_design[,  paste(drug, timing, batch, sep = "_")],
                     levels = unique(data_design[,  paste(drug, timing, batch, sep = "_")]))

# Create model matrix based on the group factor
data_permut <-  model.matrix(~ 0 + data_group)

# Rename columns and rows for clarity
colnames(data_permut) <- levels(data_group)
rownames(data_permut) <- colnames(data_grid_matrix)

# Statistical Analysis and Quality Control
# Prepare data for edgeR's DGEList
data_y <- DGEList(counts = data_grid_matrix, 
                  group = data_group, 
                  genes = row.names(data_grid_matrix))

# Filter out low counts across all samples
data_keep <- rowSums(cpm(data_y) > 0.5) >= 10

# Perform normalization and dispersion estimation
data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]
data_y <- calcNormFactors(data_y)
data_y <- estimateDisp(data_y, data_permut)

# Fit model using Quasi-Likelihood methods in edgeR
data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

# Compute Counts per Million (CPM) for further analysis
data_CPM <- cpm(data_y, prior.count = 1)
colnames(data_CPM) <- factor(data_design[,  paste(drug, timing, batch, rep, sep = "_")])

#combine keys
aba_genome <- aba_bed[aba_key[, .(spacer, type, locus_tag)], on = .(locus_tag)]

#make sure your sample names don't have spaces or periods
colnames(data_permut) <- gsub("\\.", "_", colnames(data_permut))
colnames(data_permut) <- gsub(" ", "_", colnames(data_permut))
colnames(data_permut) <- gsub("-", "_", colnames(data_permut))

# put your condition comparisons here
contrast_levels <- c('acriflavine_T2_7 - none_T2_5',
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
                     'Zinc_sulfate_T2_16 - none_T2_16')

#Make those contrasts and report both the LFC and FDRs
##########################

data_contrast <- makeContrasts(contrasts = contrast_levels, levels = data_permut)

########################

results_FDR <- aba_genome[, .(genes = unique(spacer))]
results_LFC <- aba_genome[, .(genes = unique(spacer))]

########################

for (i in 1:ncol(data_contrast)){
  
  results <- glmQLFTest(data_fit, contrast = data_contrast[,i])
  results <- topTags(results, n = Inf)
  results <- data.table(results$table)
  
  print(paste("Processing results for", contrast_levels[i], "..."))
  
  results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
  setnames(results_FDR, "FDR", contrast_levels[i])
  
  results_LFC <- results[, .(genes, logFC)][results_LFC, on = .(genes)]
  setnames(results_LFC, "logFC", contrast_levels[i])
  
}


#Put those guide-level results into a table

melted_results_FDR <-
  data.table::melt(
    results_FDR,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "FDR",
    measure.vars = contrast_levels
  )

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

################################################################################

melted_results_LFC <-
  data.table::melt(
    results_LFC,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "LFC",
    measure.vars = contrast_levels
  )

################################################################################

melted_results <-
  melted_results_LFC[
    melted_results_FDR,
    on = .(genes, condition)]

melted_results <- melted_results[!is.na(FDR)]

melted_results <-
  aba_genome[melted_results, on = .(spacer == genes)]



#Now put the gene-level results into a table
#Briefly, taking the median LFC value and recalculating Stouffer's p from guide FDRs


melted_results[
  , LFC := melted_results[
    i  = type == "control",
    j  = .(ctrl_medLFC = median(LFC, na.rm = TRUE)),
    by = .(condition)][
      i  = .SD,
      on = .(condition),
      j  = .(adj_medLFC = LFC - ctrl_medLFC),
      by = .EACHI]$adj_medLFC]



median_melted_results <-
  melted_results[, .(
    medLFC = median(LFC),
    FDR = stouffer(FDR)$p),
    by = .(locus_tag, gene_name, type, condition)]

################################################################################

median_melted_results <- curated_names[median_melted_results, on = .(AB19606 == locus_tag)]

################################################################################

setorder(median_melted_results, FDR)

################################################################################
################################################################################

CPM_melted <- melt(
  data.table(
    data_CPM, 
    keep.rownames = "spacer"), 
  id.vars = "spacer", 
  variable.name = "sample", 
  value.name = "CPM")

setorder(CPM_melted, sample)

CPM_melted <- aba_key[, .(spacer, type, locus_tag)][CPM_melted, on = .(spacer)]

CPM_melted <- aba_bed[, .(locus_tag, gene_name)][CPM_melted, on = .(locus_tag)]

CPM_melted[, condition := gsub("_[0-9]*$", "", sample)]

################################################################################


median_melted_results[, gene_name_stylized := paste0("italic('", unique_name, "')")]



#Make those gene-level results comprehensible

filter_optimized_results_LFC <- data.table::dcast(median_melted_results, AB19606 + gene_name + type  ~ condition , value.var = "medLFC")
filter_optimized_results_FDR <- data.table::dcast(median_melted_results, AB19606 + gene_name + type  ~ condition , value.var = "FDR")

#putting in columns for the string names
string_results_LFC <- curated_names[filter_optimized_results_LFC, on = .(AB19606 == AB19606)]
string_results_FDR <- curated_names[filter_optimized_results_FDR, on = .(AB19606 == AB19606)]