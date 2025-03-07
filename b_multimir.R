library(multiMiR)

# File paths
log_file <- "error_log.txt"
output_file <- "combined_results.csv"
mirna_file <- "Degs_mRNA_miRNA_Pearson_Correlation_df_negative_p_0_05.csv"
mrna_file <- "mRNA_file.csv"  # Change this to the actual file name

# Function to log messages
log_message <- function(message) {
    write(message, file = log_file, append = TRUE)
}

# Read the CSV files
mirna_data <- read.csv(mirna_file, stringsAsFactors = FALSE)
mrna_data <- read.csv(mrna_file, stringsAsFactors = FALSE)

# Merge the data on a common column (modify "Gene_ID" as needed)
merged_data <- merge(mirna_data, mrna_data, by = "Gene_ID", all.x = TRUE)

# Initialize an empty list for results
results_list <- list()

# Iterate over each row in the merged data
for (i in 1:nrow(merged_data)) {
    mirna <- merged_data$miRNA_Gene[i]
    mrna <- merged_data$mRNA_Gene[i]  # Now coming from the mRNA file
    
    # Fetch miRNA-mRNA interactions with error handling
    result <- tryCatch({
        get_multimir(org = "mmu", mirna = mirna, table = "validated")
    }, error = function(e) {
        log_message(paste("Error for miRNA:", mirna, "->", e$message))
        return(NULL)
    })
    
    # Check if results are found
    if (!is.null(result) && !is.null(result@data) && nrow(result@data) > 0) {
        results_list[[i]] <- result@data
    } else {
        msg <- paste("No results found for miRNA:", mirna)
        log_message(msg)
        results_list[[i]] <- data.frame(mRNA_Gene = mrna, miRNA_Gene = mirna, Message = "No results found")
    }
}

# Combine the results into a single data frame
combined_results <- do.call(rbind, results_list)

# Save the combined results
write.csv(combined_results, output_file, row.names = FALSE)

# Indicate completion
print("Processing complete. Check 'error_log.txt' for issues if any.")

