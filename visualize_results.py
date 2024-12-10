
import pandas as pd
import matplotlib.pyplot as plt

# Load DESeq2 results
res = pd.read_csv("DESeq2_results.csv")

# Create a Volcano plot
plt.figure(figsize=(10, 6))
plt.scatter(res['log2FoldChange'], -np.log10(res['pvalue']), c=(res['pvalue'] < 0.05), cmap='coolwarm')
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')

# Mark significant points
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')
plt.axvline(x=0, color='gray', linestyle='--')

# Save the plot
plt.savefig('volcano_plot.png')

# Show the plot
plt.show()



