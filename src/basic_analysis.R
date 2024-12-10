# Load required libraries
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

# Define function to load data
load_data <- function(file_path) {
  if (file.exists(file_path)) {
    return(read.csv(file_path))
  } else {
    stop("File not found!")
  }
}

# Load a sample dataset using the load_data function
data <- load_data("data/sample_data.csv")

# Quick summary of the data
print(summary(data))

# Create a basic scatter plot using ggplot2
ggplot(data, aes(x = Variable1, y = Variable2, color = Group)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Variable1 vs Variable2", x = "Variable 1", y = "Variable 2")

# Perform t-test to compare Variable1 between groups
result <- t.test(data$Variable1 ~ data$Group)
print(result)

# Create a heatmap (make sure the data is properly formatted for the heatmap)
# Assuming that the first column is a label, we exclude it
data_matrix <- as.matrix(data[, -1])  # Modify this as needed based on your data structure
Heatmap(data_matrix, name = "Expression", row_title = "Samples", column_title = "Features")


# Load required libraries
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(DESeq2)
library(clusterProfiler)

# Quick summary of the data
print(summary(data))

# Normalize data (example: Min-Max normalization)
normalize_data <- function(data) {
  return(as.data.frame(lapply(data, function(x) (x - min(x)) / (max(x) - min(x)))))
}

# Normalize the data (assuming first column is non-numeric)
data_normalized <- normalize_data(data[, -1])

# Perform PCA (Principal Component Analysis)
pca_result <- prcomp(data_normalized, center = TRUE, scale. = TRUE)

# Plot PCA
pca_df <- data.frame(pca_result$x)
ggplot(pca_df, aes(x = PC1, y = PC2, color = data$Group)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Normalized Data", x = "Principal Component 1", y = "Principal Component 2")

# Correlation Matrix Heatmap
cor_matrix <- cor(data_normalized)
Heatmap(cor_matrix, name = "Correlation")

# Differential Expression Analysis with DESeq2 (assuming countData and colData exist)
# Example: countData is your count matrix, colData contains metadata (e.g., condition labels)
# Assuming the data is in the form of countData and colData
# colData <- data.frame(condition = factor(c("Control", "Treatment", "Control", "Treatment")))  # Example

dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)  # Modify as needed
dds <- DESeq(dds)
res <- results(dds)

# View results of DESeq2
print(summary(res))

# Volcano plot for Differential Expression Results
volcano_data <- as.data.frame(res)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression", x = "Log2 Fold Change", y = "-log10(p-value)")

# Gene Ontology (GO) Enrichment Analysis (for DE genes with padj < 0.05)
gene_list <- rownames(res)[which(res$padj < 0.05)]
go_results <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")  # Modify based on your organism

# Plot GO enrichment results
dotplot(go_results)

# Save outputs (optional)
ggsave("pca_plot.png")
ggsave("volcano_plot.png")

