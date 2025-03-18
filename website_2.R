library(multiMiR)

predicted_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita", "targetscan")
mirna <- "mmu-miR-187-3p"

# Initialize a counter for the total number of targets
total_targets <- 0

# Loop through each database and count the targets
for (db in predicted_dbs) {
  example5 <- get_multimir(org = "mmu",
                              mirna = mirna,
                              table = db,
                              summary = FALSE, # set summary to FALSE to get the full data
                              use.tibble = TRUE)

  if (!is.null(example5@data) && nrow(example5@data) > 0) {
    num_targets <- nrow(example5@data)
    total_targets <- total_targets + num_targets
    print(paste("Database:", db, "- Targets found:", num_targets))
  } else {
    print(paste("Database:", db, "- No targets found."))
  }
}

# Print the total number of targets
print(paste("Total number of targets across all databases:", total_targets))
