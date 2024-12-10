# Load required libraries
library(ggplot2)
library(dplyr)

# Load a sample dataset
data <- read.csv("data/sample_data.csv")

# Quick summary of the data
print(summary(data))

# Create a basic plot
ggplot(data, aes(x = Variable1, y = Variable2, color = Group)) +
  geom_point() +
  theme_minimal()
