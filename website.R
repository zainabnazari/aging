library(multiMiR)

predicted_dbs <- c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita", "targetscan")
mirna <- "mmu-miR-187-3p"

# Initialize an empty vector to store results
all_types <- character()

# Loop through each database and collect the 'type' data
for (db in predicted_dbs) {
  example5 <- get_multimir(org = "mmu",
                              mirna = mirna,
                              table = db,
                              summary = TRUE,
                              use.tibble = TRUE)

  if (!is.null(example5@data)) { # Check if data exists
    all_types <- c(all_types, example5@data$type)
  }
}

# Create a table of the 'type' data
table(all_types)

