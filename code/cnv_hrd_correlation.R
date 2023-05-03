library(ggplot2)

path = '../../HRD_score/data/gene_cnv2/gene_cnv_norm/'


files = list.files(path = path, pattern = NULL, all.files = FALSE,
                   full.names = FALSE, recursive = FALSE,
                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

hrd_scores = data.frame(read.csv('../../HRD_score/data/HRD_scores_pan_cancer_annotated_typecorrect.csv'))

hrd_brca = data.frame(read.csv('../data/gene_expression/prepared_data/HRD_scores_BRCA.csv'))

hrdgenes<-read.gmt("../data/HRDgenes.gmt")


splitstring = strsplit(files, split='\\.')
case_id_gene_cnv = lapply(splitstring, function(x) paste0(x[2:(length(x)-3)], collapse= '.'))

case_id_hrd_brca = hrd_brca$case_id

common_ids = intersect(case_id_gene_cnv,case_id_hrd_brca)

brca_hrd_scores = hrd_scores[hrd_scores$Project.ID == 'TCGA-BRCA',]
brca_hrd_primary = brca_hrd_scores[brca_hrd_scores$Type == 'Primary',]

#Subt set files so only the 1077 are left, construct df, clustering
brca_files = lapply(brca_hrd_primary$case_id, function(x) paste('TCGA-BRCA',x,'gene_level_copy_number','v36','csv', sep = '.'))


common_genes_df = read.csv(paste(path,brca_files[1],sep = ''), sep = ',')
common_genes_df_hrd = subset(common_genes_df, common_genes_df$gene_name %in% hrdgenes$gene)

df = data.frame(matrix(ncol = length(brca_files), nrow = length(common_genes_df_hrd$gene_name)))
rownames(df) = common_genes_df_hrd$gene_name
colnames(df) = brca_hrd_primary$case_id

#df[,brca_hrd_primary$case_id[1]] = common_genes_df_hrd$copy_number


for (index in 1:length(brca_files)){
  cnv_df = read.csv(paste(path,brca_files[index],sep = ''), sep = ',')
  cnv_df_filtered = subset(cnv_df, cnv_df$gene_name %in% common_genes_df_hrd$gene_name)
  df[,brca_hrd_primary$case_id[index]] = cnv_df_filtered$copy_number
  if (index%%100 == 0){
    print(paste('Processed',as.character(index),sep=' '))
  }
}
write.csv(df,'../data/brca_norm_gene_cnv.csv', row.names = TRUE)

#df_raw = df + 2
#write.csv(df_raw,'../data/brca_gene_cnv.csv', row.names = TRUE)

df = read.csv('../data/brca_norm_gene_cnv.csv', sep = ',')
rownames(df) = df$X
df = df[,-1]

df_t = t(df)
df_t <- na.omit(df_t)

cases = rownames(df_t)
cases <- gsub("^X", "", cases)
cases <- gsub("\\.", "-", cases)
rownames(df_t) = cases

test = intersect(case_id_hrd_brca,cases)

brca_hrd_primary_no_na = subset(brca_hrd_primary, brca_hrd_primary$case_id %in% rownames(df_t))

df_t <- df_t[, intersect(colnames(df_t), genes_sig)]


pca = prcomp(df_t, scale. = TRUE)


plot(pca, type = "l", main = "Scree Plot")
biplot(pca, scale = 0)

plot(pca$x[,1], pca$x[,2], main = "Scores Plot")

pca_df <- data.frame(Case = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], HRD_sum = brca_hrd_primary_no_na$HRD_sum)

pca_df$HRD = ifelse(pca_df$HRD_sum >= 42, 'High', 'Low')

ggplot(pca_df, aes(x = PC1, y = PC2, color = HRD)) +
  geom_point() +
  labs(title = "PCA normalized CNV per HRD gene", x = "PC1", y = "PC2")



df_all = read.csv('../data/brca_gene_cnv.csv', sep = ',')
rownames(df_all) = df_all$X
df_all = df_all[,-1]

df_all_t = t(df_all)
df_all_t <- na.omit(df_all_t)

cases = rownames(df_all_t)
cases <- gsub("^X", "", cases)
cases <- gsub("\\.", "-", cases)
rownames(df_all_t) = cases

test = intersect(case_id_hrd_brca,cases)

brca_hrd_primary_no_na_all = subset(brca_hrd_primary, brca_hrd_primary$case_id %in% rownames(df_all_t))

pca_all = prcomp(df_all_t, scale. = TRUE)
pca_df_all <- data.frame(Case = rownames(pca_all$x), PC1 = pca_all$x[,1], PC2 = pca_all$x[,2], HRD_sum = brca_hrd_primary_no_na_all$HRD_sum)

pca_df_all$HRD = ifelse(pca_df_all$HRD_sum >= 42, 'High', 'Low')

ggplot(pca_df_all, aes(x = PC1, y = PC2, color = HRD)) +
  geom_point() +
  labs(title = "PCA normalized CNV all genes", x = "PC1", y = "PC2")
