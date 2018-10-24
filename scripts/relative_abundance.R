#'
#' Calculate relative abundances for all OTUs based on the input OTU table
#'

#' Please give the file name of the original OTU-table with taxonomic classification 

otu_file <- "all.reformat.OTU.summary.taxonomy"

# Load the tab-delimited file containing the values to be be checked (rownames in the first column)
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")


# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Save taxonomy information in vector
taxonomy <- as.vector(otu_table$Taxonomy)

# Delete column with taxonomy information in dataframe
otu_table$Taxonomy <- NULL

# Calculate relative abundances for all OTUs over all samples
# Divide each value by the sum of the sample and multiply by 100
rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))

# Reinsert the taxonomy information in relative abundance table
rel_otu_table_tax <- cbind(rel_otu_table,taxonomy)

# Write the normalized relative abundance with taxonomy table in a file and copy in directory Taxonomic-Binning if existing
write.table(rel_otu_table_tax, "OTUs_Table-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table_tax, "OTUs_Table-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))
