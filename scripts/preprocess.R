# Load required libraries
library(data.table)
library(dplyr)
library(sva)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]  # Path to the counts directory
output_folder <- args[2]  # Output folder

# Initialize an empty data frame
df <- data.frame()

# Get all file paths matching the pattern
file_paths <- list.files(path = paste0(counts_path, "/counts/"), pattern = "*", full.names = TRUE)

# Initialize batch vector
batch <- c()

# Process each file
for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    tmp <- fread(file_path)
    
    # Assign batch number based on file index
    batch <- c(batch, rep(i, ncol(tmp) - 1))  # Exclude the geneID column
    
    # Merge data frames by geneID
    if (nrow(df) == 0) {
        df <- tmp
    } else {
        df <- merge(df, tmp, by = "geneID", all = TRUE)
    }
    #print(file_path)
    #print(tmp[1:5,1:5])
}

df <- data.frame(df)
# Set geneID as row names
rownames(df) <- df$geneID
df <- df[, -1]

# Filter out rows with zero variance
df <- df[apply(df, 1, var) != 0, ]

#print(batch)
#print(df[1:5,1:5])

# Perform batch correction using ComBat-Seq
data_corrected <- ComBat_seq(as.matrix(df), batch = batch, group = NULL)

# Save intermediate file: batch-corrected data
write.table(data.frame("geneID"=rownames(data_corrected),data_corrected, check.names = FALSE), 
            file = paste0(output_folder, "/expression.batch_corrected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Log2 transform and scale the data
data_final <- log2(data_corrected + 1)  # Exclude geneID column for transformation
data_final <- t(scale(t(data_final)))
#data_final <- cbind(geneID = data_corrected$geneID, data_final)  # Add geneID back to the transformed data

# Save final processed data
write.table(data.frame("geneID"=rownames(data_final), data_final, check.names = FALSE),
                        file = paste0(output_folder, "/expression.batch.logscale.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)