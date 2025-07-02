# Load required libraries
suppressMessages(library("Mfuzz"))

# Load your merged DataFrame (with GeneID, Sample1, Sample2)
merged_d <- "/Users/zainabnazari/mfuzz/data/merged_df.xlsx"
# Read the data
datafile <- read.xlsx(merged_d)
# Assume columns are named "GeneID", "Sample1", "Sample2" (adjust if necessary)
# For soft clustering, we only need the Sample1 and Sample2 columns
exprs <- datafile[, c("FC_2_6", "FC_6_12")]

# Prepare the ExpressionSet object for Mfuzz
exprSet <- ExpressionSet(assayData=as.matrix(exprs))

# Standardize the data (optional, but usually recommended for Mfuzz)
exprSet.s <- standardise(exprSet)

# Choose the number of clusters (nb_clusters) based on your needs (e.g., 4)
nb_clusters <- 4

# Perform soft clustering using Mfuzz
m1 <- mestimate(exprSet.s)
cl <- mfuzz(exprSet.s, c=nb_clusters, m=m1)
getwd()

# Create the output directory
output <- "cluster_output"
dir.create(output, showWarnings = TRUE)

# For each membership value, generate output
membership_cutoff <- 0.6
dir <- paste(output, paste("cluster_with_membership", membership_cutoff, sep="_"), sep="/")
dir.create(dir, showWarnings = FALSE)

# Plot clusters
pdf(paste(dir, "clusters_Mfuzz_membership_equals_0.7.pdf", sep="/"))
mfuzz.plot2(exprSet.s, cl=cl, time.labels=colnames(exprs), min.mem=membership_cutoff, colo="fancy", x11=FALSE)
dev.off()

# Generate gene lists per cluster
acore.list <- acore(exprSet.s, cl=cl, min.acore=membership_cutoff)

# Write gene lists for each cluster to files
for (cluster in 1:nb_clusters) {
  print(paste("Number of genes in cluster", cluster, ":", dim(acore.list[[cluster]])[1]))
  cluster_table <- merge(datafile, acore.list[[cluster]][2], by="row.names", all.y=TRUE)
  write.table(cluster_table, paste(dir, paste("list_of_genes_in_cluster", cluster, ".txt", sep="_"), sep="/"), sep="\t", row.names=FALSE, dec=".")
}

