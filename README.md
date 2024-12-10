## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


# Bioinformatics Tool

This tool performs a variety of bioinformatics analyses, including PCA, differential expression analysis, and GO enrichment.

## Prerequisites
- R
- Required libraries: `ggplot2`, `dplyr`, `ComplexHeatmap`, `DESeq2`, `clusterProfiler`.

## Installation
Install required libraries by running:
```R
install.packages(c("ggplot2", "dplyr", "ComplexHeatmap", "DESeq2", "clusterProfiler"))

## Shiny App

This tool includes a Shiny app for interactive data analysis. To run the app:

1. Install required R packages if you haven't already:
   ```r
   install.packages(c("shiny", "ggplot2", "dplyr"))


## Linux Analysis Scripts

This repository includes a set of Linux-based analysis scripts located in the `linux_analysis` folder.

### Running Data Preprocessing

To run the data preprocessing script:

```bash
bash linux_analysis/data_preprocessing.sh
