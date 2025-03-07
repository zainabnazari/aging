# Load required libraries
library(multiMiR)
library(readr)  # For reading CSV files

# Step 1: Read the CSV file into a data frame
# Replace 'your_file.csv' with the path to your CSV file
df <- read.delim('/Users/zainabnazari/data_work/only_6M_12M_miRNA_DEGs_FDR_5pc_1_5fc_Bowtie1_align_miRBase_v_0_m_1_only_R1_samtools_idxstat_filtered_Log2CPM_TMM_Norm.txt', stringsAsFactors = FALSE)

# Extract miRNA names from the 'SampleID' column
mirna_names <- df$SampleID

# Function to extract targets
extract_targets <- function(mirna_list, dbs_exp, dbs_pred) {
  results <- list()

  for (mirna in mirna_list) {
    # Initialize an empty set for experimental targets
    experimental_targets <- character()

    # Step 1: Extract experimental validated targets
    for (db in dbs_exp) {
      # Fetch targets from the database
      targets <- get_multimir(org = "mmu", mirna = mirna, table = db)
      
      # Get top 20% targets based on scores
      if (!is.null(targets@data) && nrow(targets@data) > 0) {
        top_n <- ceiling(nrow(targets@data) * 0.2)  # Calculate top 20%
        top_targets <- targets@data[order(-targets@data$score), ][1:top_n, ]
        
        # Pool targets using set union
        experimental_targets <- union(experimental_targets, top_targets$mRNA_Gene)
      }
    }

    # Step 2: Extract in silico predicted targets
    predicted_targets <- character()
    
    # Fetch targets from all in silico databases
    for (db in dbs_pred) {
      targets <- get_multimir(org = "mmu", mirna = mirna, table = db)
      if (!is.null(targets@data) && nrow(targets@data) > 0) {
        predicted_targets <- union(predicted_targets, targets@data$mRNA_Gene)
      }
    }

    # Count occurrences of predicted targets across databases
    target_counts <- table(predicted_targets)
    
    # Keep targets present in at least 3 databases
    likely_targets <- names(target_counts[target_counts >= 3])

    # Combine experimental and likely targets
    results[[mirna]] <- unique(c(experimental_targets, likely_targets))
  }

  return(results)
}

# Step 2: Define databases
experimental_dbs <- c("miRecords", "miRTarBase", "TarBase")
predicted_dbs <- c("DIANA-microT-CDS", "ElMMo", "MicroCosm", "miRanda", "miRDB", "PicTar", "PITA", "TargetScan")

# Step 3: Extract targets for the DEGs
miRNA_targets <- extract_targets(mirna_names, experimental_dbs, predicted_dbs)

# Step 4: Display the results
print(miRNA_targets)

# Save the results to a CSV file if needed
results_df <- do.call(rbind, lapply(names(miRNA_targets), function(mirna) {
  data.frame(mirna = mirna, targets = paste(miRNA_targets[[mirna]], collapse = ", "), stringsAsFactors = FALSE)
}))
write.csv(results_df, 'miRNA_targets_results.csv', row.names = FALSE)

