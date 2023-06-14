args <- commandArgs(trailingOnly = TRUE)

# Read in the file
results <- readLines(args[1])

# Create a list to store the results
output <- list()

# Loop through each line of the results
for (i in 1:length(results)) {
    # If the line starts with "L", extract the module and cluster information
    if (startsWith(as.character(results[i]), "L")) {
        modClus <- unlist(strsplit(results[i], "_"))
        module <- modClus[1]
        cluster <- modClus[2]

        # Add the module to the output list if it doesn't exist
        if (!(module %in% names(output))) {
            output[[module]] <- list()
        }

        # Add the cluster to the output list for the current module
        output[[module]] <- c(output[[module]], cluster)
    }
}

# Write the output to a file
writeLines(paste0(names(output), ": ", sapply(output, paste, collapse = ", ")), args[2])
