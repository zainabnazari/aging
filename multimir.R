#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("multiMiR")
library(multiMiR)


#tables <- all_tables()
#print(tables)



# Read the CSV file into a data frame
# Replace 'your_file.csv' with the path to your CSV file
data <- read.csv('/Users/zainabnazari/data_work/only_6M_12M_miRNA_DEGs_FDR_5pc_1_5fc_Bowtie1_align_miRBase_v_0_m_1_only_R1_samtools_idxstat_filtered_Log2CPM_TMM_Norm.txt', stringsAsFactors = FALSE)

# Initialize an empty list to store results
results_list <- list()

# Iterate over each row in the data frame
for (i in 1:nrow(data)) {
    mirna <- data$SampleID[i]
    
    # Fetch miRNA-mRNA interactions
    result <- get_multimir(org = "mmu", mirna = mirna, table = "predicted")  # Use other tables if needed

    # Check if results are found
    if (!is.null(result@data) && nrow(result@data) > 0) {
        # Append the results to the list
        results_list[[i]] <- result@data
    } else {
        # If no results, store a message or NA
        results_list[[i]] <- data.frame(SampleID = data$SampleID[i], SampleID = mirna, Message = "No results found")
    }
}

# Combine the results into a single data frame
combined_results <- do.call(rbind, results_list)

# Save the combined results to a CSV file if needed
write.csv(combined_results, 'combined_results.csv', row.names = FALSE)

# Display the first few rows of the combined results
print(head(combined_results))



print('done!')