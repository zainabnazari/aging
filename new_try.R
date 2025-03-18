# Load required library
library(multiMiR)

# Define the miRNA you want to search for
mirna <- "mmu-miR-21a-5p"  # Change this to your desired miRNA

# Define predicted target databases
predicted_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda", 
                   "mirdb", "pictar", "pita", "targetscan")

# Fetch predicted target genes using multiMiR
result <- get_multimir(org = "mmu", mirna = mirna, table = "predicted")

# Check if results exist
if (!is.null(result@data) && nrow(result@data) > 0) {
  # Extract unique target genes
  target_genes <- unique(result@data$mRNA_Gene)
  
  # Print the target genes
  print(target_genes)
  
  # Save to CSV
  write.csv(target_genes, "predicted_target_genes.csv", row.names = FALSE)
} else {
  print("No predicted target genes found for this miRNA.")
}

