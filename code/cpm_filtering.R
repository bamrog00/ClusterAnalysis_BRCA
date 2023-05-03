library(edgeR)

# read in count data
counts = read.table("../data/gene_expression/firebrowse_BRCA/BRCA.mRNAseq_raw_counts.txt",sep= '\t', header = TRUE, row.names = 1)

# calculate CPM values
cpm = cpm(counts)

# filter out genes with CPM < 1 in at least 2 samples
keep = filterByExpr(cpm, min.mean=1, N=2)

# subset count and CPM matrices to keep only selected genes
counts_filt = counts[keep,]
counts_filt_not = counts[!keep,]
cpm_filt = cpm[keep,]

write.csv(counts_filt, "../data/gene_expression/firebrowse_BRCA/BRCA_filtered_counts.csv")
