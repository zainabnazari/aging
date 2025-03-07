# Load required library
library(multiMiR)

# Read the CSV file into a data frame
data <- read.csv('/Users/zainabnazari/data_work/only_6M_12M_miRNA_DEGs_FDR_5pc_1_5fc_Bowtie1_align_miRBase_v_0_m_1_only_R1_samtools_idxstat_filtered_Log2CPM_TMM_Norm.txt', 
                 stringsAsFactors = FALSE)

# Initialize an empty list to store results
results_list <- list()

# Define experimentally validated databases (correct lowercase names)
experimental_dbs <- c("mirecords", "mirtarbase", "tarbase")

# Define predicted databases (correct lowercase names)
predicted_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita", "targetscan")

# Function to extract top 20% targets from a given database
get_top_targets <- function(mirna, db, percentage = 0.2) {
  result <- get_multimir(org = "mmu", mirna = mirna, table = db)
  
  if (!is.null(result@data) && nrow(result@data) > 0) {
    top_n <- ceiling(nrow(result@data) * percentage)  # Select top 20% of targets
    top_targets <- result@data[order(-result@data$score), ][1:top_n, ]
    return(unique(top_targets$mRNA_Gene))  # Return unique mRNA target genes
  }
  return(character())  # Return an empty character vector if no results
}

# Iterate over each miRNA in the dataset
for (i in 1:nrow(data)) {
  mirna <- data$SampleID[i]
  
  # Step 1: Get validated targets (appear in at least 1 experimental database)
  validated_targets <- character()
  for (db in experimental_dbs) {
    validated_targets <- union(validated_targets, get_top_targets(mirna, db))
  }
  
  # Step 2: Get predicted targets (appear in at least 3 predicted databases)
  predicted_targets <- list()
  for (db in predicted_dbs) {
    targets <- get_multimir(org = "mmu", mirna = mirna, table = db)
    if (!is.null(targets@data) && nrow(targets@data) > 0) {
      predicted_targets[[db]] <- unique(targets@data$mRNA_Gene)
    }
  }
  
  # Pool all predicted targets into one list
  all_predicted_targets <- unlist(predicted_targets, use.names = FALSE)
  
  # Count occurrences of each predicted target across different databases
  target_counts <- table(all_predicted_targets)
  
  # Select only targets that appear in at least 3 different predicted databases
  likely_predicted_targets <- names(target_counts[target_counts >= 3])
  
  # Combine validated and predicted targets
  combined_targets <- unique(c(validated_targets, likely_predicted_targets))
  
  # Store results for the current miRNA
  results_list[[i]] <- data.frame(SampleID = mirna, 
                                  mRNA_Genes = paste(combined_targets, collapse = ", "), 
                                  stringsAsFactors = FALSE)
}

# Combine the results into a single data frame
final_results <- do.call(rbind, results_list)

# Save the combined results to a CSV file
write.csv(final_results, 'miRNA_targets_results.csv', row.names = FALSE)

# Display the first few rows of the final results
print(head(final_results))

print('done!')

