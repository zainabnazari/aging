library(multiMiR)

rm(list = ls())
options("warn"=0)  # Print max 10 warnings on screen

# Define thresholds
percent_cutoff = 20  ## Take top 20% scored targets of each database
min_num_predicted_DBs_for_target = 3  ## Target must appear in at least this many predicted databases
min_num_validated_DBs_for_target = 1  ## Target must appear in at least this many validated databases

# File paths
work_path = "/Users/zainabnazari/work/"
name_file_input_miRNAs = "input_microRNAs_list.txt"

setwd(work_path)

# Read miRNA list
miRNA_list = read.table(file=name_file_input_miRNAs, sep = "\t", quote = "\"", header=FALSE, fill=TRUE)

# Initialize result lists
all_predicted_miRNA_targets = list()
all_validated_miRNA_targets = list()

# Process each miRNA
for (ii in 1:nrow(miRNA_list)) {
    curr_mirna = miRNA_list[ii, 1]

    print(paste("Processing Line:", ii, ", miRNA:", curr_mirna))

    # Fetch predicted targets with error handling
    curr_predicted_targets = tryCatch({
        get_multimir(org = "mmu", mirna = curr_mirna, table = "predicted",
                     predicted.cutoff = percent_cutoff, predicted.cutoff.type = "p",
                     predicted.site = "all", summary = FALSE, add.link = FALSE,
                     use.tibble = FALSE, limit = NULL, legacy.out = FALSE)@data
    }, error = function(e) NULL)

    # Fetch validated targets with error handling
    curr_validated_targets = tryCatch({
        get_multimir(org = "mmu", mirna = curr_mirna, table = "validated",
                     summary = FALSE, add.link = FALSE, use.tibble = FALSE,
                     limit = NULL, legacy.out = FALSE)@data
    }, error = function(e) NULL)

    # Handle NULL cases
    if (is.null(curr_predicted_targets)) curr_predicted_targets <- data.frame(target_symbol=character())
    if (is.null(curr_validated_targets)) curr_validated_targets <- data.frame(target_symbol=character())

    # Count target occurrences across databases
    predicted_counts = table(curr_predicted_targets$target_symbol)
    validated_counts = table(curr_validated_targets$target_symbol)

    # Filter targets by minimum database presence
    predicted_targets = names(predicted_counts[predicted_counts >= min_num_predicted_DBs_for_target])
    validated_targets = names(validated_counts[validated_counts >= min_num_validated_DBs_for_target])

    # Store results
    all_predicted_miRNA_targets[[ii]] = predicted_targets
    all_validated_miRNA_targets[[ii]] = validated_targets
}

# Assign names
names(all_predicted_miRNA_targets) = miRNA_list[,1]
names(all_validated_miRNA_targets) = miRNA_list[,1]

# Save predicted targets
max_length_pred = max(sapply(all_predicted_miRNA_targets, length))
predicted_padded = lapply(all_predicted_miRNA_targets, function(x) c(x, rep("\t", max_length_pred - length(x))))
write.table(as.data.frame(predicted_padded), "out_predicted_targets.txt", sep="\t", row.names=FALSE, col.names=TRUE)

# Save validated targets
max_length_valid = max(sapply(all_validated_miRNA_targets, length))
validated_padded = lapply(all_validated_miRNA_targets, function(x) c(x, rep("\t", max_length_valid - length(x))))
write.table(as.data.frame(validated_padded), "out_validated_targets.txt", sep="\t", row.names=FALSE, col.names=TRUE)

print("Processing complete. Results saved.")

