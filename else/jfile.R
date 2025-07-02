# Install required packages if not already installed
packages <- c("readxl", "dplyr", "tidyr", "Mfuzz", "BiocManager", "ggplot2", "tibble")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install)) install.packages(to_install)

# Mfuzz needs specific Bioconductor packages
if (!require("Mfuzz", quietly = TRUE)) {
  BiocManager::install("Mfuzz")
}

# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(Mfuzz)
library(tibble)

# Read Excel files
df_2_6 <- read_excel("miRNA_2M_6M.xlsx")
df_6_12 <- read_excel("miRNA_6M_12M.xlsx")

# Ensure 'sequence' column is used as gene ID
df_2_6 <- df_2_6 %>% filter(!is.na(sequence)) %>%
  select(sequence, log2FoldChange) %>%
  rename(logFC_6M = log2FoldChange)

df_6_12 <- df_6_12 %>% filter(!is.na(sequence)) %>%
  select(sequence, log2FoldChange) %>%
  rename(logFC_12M = log2FoldChange)

# Merge both by 'sequence'
merged <- full_join(df_2_6, df_6_12, by = "sequence")

# Add baseline timepoint (2M = log2FC = 0)
merged <- merged %>%
  mutate(logFC_2M = 0) %>%
  relocate(logFC_2M, .after = sequence)

# Remove rows with missing values
merged_clean <- merged %>%
  drop_na()

# Prepare expression matrix for Mfuzz
expr_matrix <- merged_clean %>%
  column_to_rownames("sequence") %>%
  as.matrix() %>%
  t()

# Create ExpressionSet for Mfuzz
library(Biobase)
eset <- ExpressionSet(assayData = expr_matrix)

# Standardize (normalize) the data
eset_std <- standardise(eset)

# Estimate optimal number of clusters (optional)
m <- mestimate(eset_std)

# Set number of clusters (you can experiment with this)
c <- 6  # try 4 to 8 clusters depending on data

# Run Mfuzz clustering
cl <- mfuzz(eset_std, c = c, m = m)

# Plot clusters
png("mfuzz_clusters.png", width = 1000, height = 800)
mfuzz.plot(eset_std, cl = cl, mfrow = c(2, 3), time.labels = c("2M", "6M", "12M"))
dev.off()

# Save cluster membership
cluster_membership <- tibble(sequence = rownames(merged_clean),
                             cluster = cl$cluster)
write.csv(cluster_membership, "mfuzz_cluster_membership.csv", row.names = FALSE)

