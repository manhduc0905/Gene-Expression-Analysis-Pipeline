library(tidyverse)

#reading files
metadata <- read_csv("Processed_Data/metadata.csv")
expr <- read_csv("Processed_Data/normalized_expression.csv")
expr$...1 <- NULL

#differentiate between "leukemia" and "normal"
samples <- colnames(expr)[3:ncol(expr)]
is_normal <- metadata$condition == "normal"

#calculate mean & l2fc
#l2fc : leukemia/normal but since I normalized data already (applied log2), I
#subtract here

mean_normal <- rowMeans(expr[3:ncol(expr)][,is_normal])
mean_leukemia <- rowMeans(expr[3:ncol(expr)][,!is_normal])
log2_foldchange <- mean_leukemia - mean_normal

#rearrange the data
expr <- cbind(expr, log2_foldchange)
expr <- expr%>%
  select(ID, gene_symbol, log2_foldchange, everything())

#p-test
p_values <- apply(expr[, 3:ncol(expr)], 1, function (x) {t.test(x[!is_normal], x[is_normal])$p.value})

# Calculate adjusted p-values using the "fdr" (False Discovery Rate) method
adj_p_values <- p.adjust(p_values, method = "fdr")


# Check the first few
expr <- cbind(expr, p_values)
expr <- cbind(expr, adj_p_values)

expr <- expr%>%
  select(ID, gene_symbol, log2_foldchange, adj_p_values,p_values,  everything())

significant_genes <- expr%>%
  filter(adj_p_values < 0.05 & abs(log2_foldchange) > 1)

write.csv(significant_genes,
          "Processed_Data/significant_genes.csv",
          row.names = TRUE)

#sorting to order
up_regulated <- significant_genes[order(significant_genes$log2_foldchange, decreasing = TRUE), ]
down_regulated <- significant_genes[order(significant_genes$log2_foldchange), ]

#More expression of genes in normal sample
write.csv(down_regulated,
          "Processed_Data/l2fc_normal_sorting.csv",
          row.names = TRUE)
#More expression of genes in leukemia sample
write.csv(up_regulated,
          "Processed_Data/l2fc_leukemia_sorting.csv",
          row.names = TRUE)

write.csv(expr,
          "Processed_Data/expr_stats.csv",
          row.names = TRUE)
