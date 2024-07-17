require('pacman')

library("ggallin", "stringr", "tidyverse", "janitor")


p_load(stringr, clipr, data.table, dplyr, ggallin, scales, edgeR, statmod, poolr, 
       pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, 
       RColorBrewer, purrr, naniar)


## #read data file into R
# data<-read.table("231018_LFC_t2_position_count.tsv", sep = "\t", header = TRUE)
##read in curated_names
# curated_names <- fread("curated_names.tsv")

#or alternatively take it straight from EdgeR script
data <- string_results_LFC


####pairwise_correlation_matrix_in_R###
#organize data file

data$condition <- str_c(data$unique_name, "_", data$type)
data <- subset(data, select = -c(1:9))
colnames(data) <- gsub("_T2.*", "", colnames(data))
data <- t(data)
colnames(data) <- data[c("condition"),]
data <- data[, colnames(data) != ""]
data <- data[rownames(data) != "condition", ]

#filter perfect guides only
perfectdf <- data[, colnames(data) %like% "_perfect"] 
colnames(perfectdf) <- gsub("_perfect", "", colnames(perfectdf))

#pull the numerical data from the data frame as a matrix
mdf<-data.matrix(perfectdf)
mdf <- matrix(as.numeric(mdf),   
              ncol = ncol(mdf))
colnames(mdf) <- colnames(perfectdf)
rownames(mdf) <- rownames(perfectdf)

##transpose for chemical network
# mdf <- t(mdf)

#create a Pearson correlation matrix from the data, can be calculated for columns with "NA" values
correlations <- cor(mdf, use="pairwise.complete.obs", method = "pearson")

# melt
correlations_melt <- reshape2::melt(correlations)


#randomize LFCs with 5000 permutations, makes list of data frames
random<-replicate(5000, apply(mdf,2,sample), simplify=FALSE)

#find correlations between random s-scores in each dataframe
random_cor<- lapply(random, cor)

#find the maximum correlation value in each dataframe that is not 1
max_test<-rep(0,5000)
for (i in 1:5000) {
  max_test[i] <- max(random_cor[[i]][random_cor[[i]]!=1])
}



#percentiles
# q75 <- quantile (max_test, 0.75)
# q80 <- quantile (max_test, 0,8)
q85 <- quantile (max_test, 0.85)
# q90 <- quantile(max_test, 0.9)
# q95 <- quantile(max_test, 0.95)
# q99 <- quantile(max_test, 0.99)

#85th percentile cutoff!
melt_correlations_85 <- correlations_melt %>% filter(value > q85) %>% filter(value != 1) 
melt_correlations_85_string <- melt_correlations_85 %>%
  left_join(curated_names, by = c("Var1" = "unique_name")) %>% select(`STRING_names`, AB19606, Var1, Var2, value) %>%
  left_join(curated_names, by = c("Var2" = "unique_name")) %>% select(`STRING_names.x`, AB19606.x, Var1, `STRING_names.y`, AB19606.y, Var2, value)

#write a table for input into Cytoscape
write.table(melt_correlations_85_string, "meltcorrelations85_chemicals.tsv", 
            sep="\t", row.names=FALSE, col.names= TRUE, quote=FALSE)

##############
##############
#add in operon information for Cytoscape network formatting

#Read in operon info
operon <- fread("curated_names_operons_pathways.tsv")

# Create a mapping table for AB19606 to operon
mapping_table <- operon %>% select(AB19606, operon)

# Perform the joins and create the operon column
melt_correlations_85_string <- melt_correlations_85_string %>%
  left_join(mapping_table, by = c("AB19606.x" = "AB19606")) %>%
  rename(operon_x = operon) %>%
  left_join(mapping_table, by = c("AB19606.y" = "AB19606")) %>%
  rename(operon_y = operon) %>%
  mutate(operon = ifelse(operon_x == operon_y, 1, 0)) %>%
  select(-operon_x, -operon_y)

operon_df <- melt_correlations_85_string %>%
  mutate(shared_name = paste(Var1, "(interacts with)", Var2)) %>%
  select(shared_name, operon)


