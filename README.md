# ClusterAnalysis_BRCA

### Python

## mutation_data_preparation (need to change name)
Can be used to check the correlation between the HRD scores (HRDsum, TAI, LST, LOH) and the mutations in HRR genes.
Currently it only checks if there is at least one mutation in any of the HRR genes.

## data_preparation
Prepares the count matrix, clinical data and HRD score file for differential expression analysis. The three dataframes are matched, so only cases that are in all three
dataframes are used. The count matrix and the clinical data are prepared in such a way that it matches the input for DeSeq2. The clinical data is the one from cBioportal which contains the mol. Subtypes.

## Clustering_testing
Not done yet

## checkout_match_data
Is now included in data_preparation but also compares with other RNAseq sources
Not documented yet

