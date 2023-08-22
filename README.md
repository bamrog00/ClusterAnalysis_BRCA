# ClusterAnalysis_BRCA

### Python

## mutation_data_preparation
Can be used to check the correlation between the HRD scores (HRDsum, TAI, LST, LOH) and the mutations in HRR genes.
Currently, it only checks if there is at least one mutation in any of the HRR genes.

## data_preparation
Prepares the count matrix, clinical data and HRD score file for differential expression analysis. The three dataframes are matched, so only cases that are in all three
data frames are used. The count matrix and the clinical data are prepared in such a way that it matches the input for DeSeq2. The clinical data is the one from cBioportal which contains the mol. Subtypes.
For BRCA, LUSC, LUAD and OV.

## Clustering_testing
Different clustering methods used to cluster the subtypes such as PCA and others.

## checkout_match_data
Is now included in data_preparation but also compares with other RNAseq sources
Was used to compare the RNAseq data from different sources. Not really used any more as the STAR count data from TCGA was used in the end for the pan-cancer analysis.

## General
This part was about the correlation of the subtypes of BRCA with HRD score and try to cluster the data.

