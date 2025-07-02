# Load required libraries
suppressMessages(library("Mfuzz"))
library(openxlsx)  # For reading .xlsx files
library(Biobase)

# Load your merged DataFrame (with GeneID, FC_2_6, FC_6_12)
merged_d <- "/Users/zainabnazari/mfuzz/data/merged_df.xlsx"
datafile <- read.xlsx(merged_d)

# Keep only the fold changes and set time points: 0, 1, 2
exprs <- datafile[, c("FC_2_6", "FC_6_12")]
# Insert a zero column for the baseline
exprs <- cbind(Baseline=rep(0, nrow(exprs)), exprs)
rownames(exprs) <- datafile$GeneID  # Set gene IDs as rownames

# Prepare the ExpressionSet object
exprSet <- ExpressionSet(assayData = as.matrix(exprs))

# Standardize the data
exprSet.s <- standardise(exprSet)

# Clustering
nb_clusters <- 3
m1 <- mestimate(exprSet.s)
cl <- mfuzz(exprSet.s, c=nb_clusters, m=m1)

# Output directory
output <- "/Users/zainabnazari/mfuzz/cluster_output"
dir.create(output, showWarnings = TRUE)
membership_cutoff <- 0.7
dir <- file.path(output, paste("cluster_with_membership", membership_cutoff, sep="_"))
dir.create(dir, showWarnings = FALSE)

# Plot all clusters in a single page using mfuzz.plot
pdf(file.path(dir, "clusters_all_in_one_page.pdf"), width=10, height=8)
mfuzz.plot(exprSet.s, cl=cl, mfrow=c(2, 2), time.labels=c("0", "1", "2"), new.window=FALSE)
dev.off()



