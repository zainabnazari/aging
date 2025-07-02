# Load required libraries
suppressMessages(library("Mfuzz"))
library(openxlsx)  # For reading .xlsx files
library(Biobase)

# Load your merged DataFrame (with GeneID, FC_2_6, FC_6_12)
merged_d <- "/Users/zainabnazari/mfuzz/data/merged_df.xlsx"
datafile <- read.xlsx(merged_d)

# Keep only the fold changes and set time points: 0, 1, 2
exprs <- datafile[, c("FC_2_6", "FC_6_12")]

#  Standardize ONLY the FC_2_6 and FC_6_12 columns
exprs_scaled <- scale(exprs)

# Insert a zero column for the baseline
exprs_scaled <- cbind(Baseline=rep(0, nrow(exprs_scaled)), exprs_scaled)

# Set rownames as Gene IDs
rownames(exprs_scaled) <- datafile$sequence

# Prepare the ExpressionSet object
exprSet <- ExpressionSet(assayData = as.matrix(exprs_scaled))

# Clustering
nb_clusters <- 5
m1 <- mestimate(exprSet)  # Now m is estimated on standardized data
cl <- mfuzz(exprSet, c=nb_clusters, m=m1)

# Output directory
output <- "/Users/zainabnazari/mfuzz/specific_st/cluster_output"
dir.create(output, showWarnings = TRUE)
membership_cutoff <- 1.0
dir <- file.path(output, paste("cluster_with_membership", membership_cutoff, sep="_"))
dir.create(dir, showWarnings = FALSE)

# Plot all clusters in a single page using mfuzz.plot
pdf(file.path(dir, "clusters_cutoff_1_n_5_try_0.pdf"), width=10, height=8)
mfuzz.plot(exprSet, cl=cl, mfrow=c(2, 2), time.labels=c("0", "2_6", "6_12"), new.window=FALSE)
dev.off()

