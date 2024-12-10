
# Load DESeq2 library
library(DESeq2)

# Load count data and sample information
counts <- read.csv("count_data.csv", row.names = 1)
coldata <- read.csv("sample_data.csv", row.names = 1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file = "DESeq2_results.csv")
