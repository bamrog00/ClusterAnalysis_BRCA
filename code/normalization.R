filtered = data.frame(read.csv('../data/gene_expression/firebrowse_BRCA/BRCA_filtered_counts.csv'))
gene_lengths = read.table("../data/gene_expression/firebrowse_BRCA/BRCA.mRNAseq_median_length_normalized.txt",sep= '\t', header = TRUE, row.names = 1)
rpkm_firebrowse = read.table("../data/gene_expression/firebrowse_BRCA/BRCA.mRNAseq_RPKM.txt",sep= '\t', header = TRUE, row.names = 1)

genes_filtered = filtered[,1]
print(genes_filtered)

gene_lengths_filtered = gene_lengths[genes_filtered,]

genes_in_lengths = rownames(gene_lengths_filtered)

print(identical(genes_filtered,genes_in_lengths))

row.names(filtered) = filtered[,1]
filtered = filtered[-1]

#Replace NA
na_matrix <- is.na(filtered)

# Replace all NA values with 0
filtered[na_matrix] <- 0


## RPKM
pm = colSums(filtered) / 1000000
pm_matrix <- matrix(pm, nrow = nrow(filtered), ncol = ncol(filtered), byrow = TRUE)
rpm = filtered/ pm_matrix

print(rpm[1,2])
print(filtered[1,2] / pm_matrix[1,2])

rpkm_df = rpm / gene_lengths_filtered

## TPM
rpk_df = filtered / gene_lengths_filtered

# Change NA to 0 in rpk_df
rpk_sums = colSums(rpk_df)

print(rpk_sums)
# Sums are not the same because the lengths vary, other normalization ?
# TODO divide by the sum